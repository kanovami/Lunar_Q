# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.rc('font',family='Times New Roman',size=14)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

import os
import numpy as np
import emcee
import matplotlib.pyplot as plt
import scipy.constants as sc
from multiprocessing import Pool
import corner
from tidal_deformation import mrheo,mconst

os.environ["OMP_NUM_THREADS"] = "1"
max_n = 10000000

def mass(rho_mantle, rho_core, rho_crust, r_core, d_crust):
    r_mantle = r_planet-d_crust
    M = r_core**3*(rho_core-rho_mantle) + r_mantle**3*(rho_mantle-rho_crust) + r_planet**3*rho_crust
    return 4/3*sc.pi*M

def MoIF(rho_mantle, rho_core, rho_crust, r_core, d_crust):
    r_mantle = r_planet-d_crust
    moif = r_core**5*(rho_core-rho_mantle) + r_mantle**5*(rho_mantle-rho_crust) + r_planet**5*rho_crust
    moif /= r_core**3*(rho_core-rho_mantle) + r_mantle**3*(rho_mantle-rho_crust) + r_planet**3*rho_crust
    moif *= 2/5/r_planet**2
    return moif

def k2(mu, eta, alpha, zeta, Delta, trel, rho_mantle, rho_core, rho_crust, 
       r_core, d_crust, omg):
    
    etakv_mantle = trel/Delta
    
    # Allocate Fortran arrays
    # Planet interior structure - layerwise, from core to surface
    #-------------------------------------------------------------
    mconst.ire  = [0, 2, 0]               # rheological models
    mconst.rho  = [rho_core, rho_mantle, rho_crust]         # densities
    mconst.arad = [r_core, r_planet-d_crust, r_planet]# radii of interfaces
    mconst.eta  = [1e0, eta, 0]       # viscosities
    mconst.mu   = [1e-10, mu, mu]      # rigidities
    mconst.alpha= [None, alpha, None]      # parameter of the Andrade model
    mconst.zeta = [None, zeta, None]         # parameter of the Andrade model
    mconst.eta_kv_rel = [None, etakv_mantle, None] # parameter of the Sundberg-Cooper model
    mconst.cdef = [None, Delta, None]     # parameter of the Sundberg-Cooper model
    #-------------------------------------------------------------

    mconst.nl = len(mconst.arad)
    mconst.iellid = 1

    k2, h2 = mrheo.peltier(2, omg)
    k3, h3 = mrheo.peltier(3, omg)
    return k2,h2,k3

# Logarithm of a prior distribution
def log_prior(theta):
    mu, eta, alpha, zeta, Delta, trel, rho_mantle, rho_core, r_core = theta
    mu_cond = 60 < mu < 90
    eta_cond = 15 < eta < 30
    alpha_cond = 0.1 < alpha < 0.4
    zeta_cond = -5 < zeta < 5
    Delta_cond = -2 < Delta < 1
    trel_cond = -10 < trel < 0
    rhom_cond = 3 < rho_mantle < 4
    rhoc_cond = 4 < rho_core < 7
    Rc_cond = 0 < r_core < 450   
    if (mu_cond and eta_cond and alpha_cond and zeta_cond and Delta_cond 
        and trel_cond and rhom_cond and rhoc_cond and Rc_cond):
        return 0.0
    return -np.inf

def log_likelihood(theta, xobs, sig):
    '''
    Logarithm of probability density function
    '''  
    k_Love_m, h_Love_m, k_Love3_m = k2(
                theta[0]*1e9, 10**theta[1], theta[2], 10**theta[3],
                10**theta[4], 10**theta[5], theta[6]*1e3, theta[7]*1e3,
                2550, theta[8]*1e3, 40e3, omg_month)
    k_Love_a, h_Love_a, k_Love3_a = k2(
                theta[0]*1e9, 10**theta[1], theta[2], 10**theta[3],
                10**theta[4], 10**theta[5], theta[6]*1e3, theta[7]*1e3,
                2550, theta[8]*1e3, 40e3, omg_year)
    
    k2mod = k_Love_m.real
    k3mod = k_Love3_m.real
    h2mod = h_Love_m.real
    Qmmod = 1/np.sin(abs(np.angle(k_Love_m)))
    #Qamod = 1/np.sin(abs(np.angle(k_Love_a)))
    k2Q_a = -k_Love_a.imag
    MoIFmod = MoIF(theta[6]*1e3, theta[7]*1e3, 2550, theta[8]*1e3, 40e3)
    massmod = mass(theta[6]*1e3, theta[7]*1e3, 2550, theta[8]*1e3, 40e3)
    
    xmod = np.array([k2mod, Qmmod, k2Q_a, k3mod, h2mod, MoIFmod, massmod])
   
    a  = (xmod-xobs)
    b  = np.diag(1/sig**2)
    S  = np.matmul(np.matmul(a,b), np.transpose(a))
    return -0.5 * S

# Full probability given priors
def log_probability(theta, xobs, sig):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, xobs, sig)

#===============================
# Observed values
#===============================
# monthly tide
Q_month   = 38
sQ_month  = 4 #0.4
k2_month  = 0.02422
sk2_month = 0.00022
h2_month  = 0.0387  # Thor et al. (2021)
sh2_month = 0.0025
k3_month  = 0.0081  # unweighted mean, Konopliv+13, Lemoine+13
sk3_month = 0.0018
# annual tide
# Q_year    = 41
# sQ_year   = 0.1
k2Q_year  = 6.2e-4
sk2Q_year = 1.4e-4 #0.06e-4
# moment of inertia factor
MoIF_mean  = 0.392728  # Williams et al. (2014)
sMoIF_mean = 0.000012
# mass
M_planet  = 7.34630e22 # Williams et al. (2014)
sM_planet = 0.00088e22
#===============================

#===============================
# Other parameters
#===============================
omg_month = 2*sc.pi/(27.212*86400) # monthly
omg_year = 2*sc.pi/(365.260*86400) # annual
r_planet = 1737.151e3 # Wieczorek (2015)
#===============================

#--------------------------------
# Set values of hyperparameters
#--------------------------------
means = np.array([k2_month, Q_month, k2Q_year, k3_month, h2_month, MoIF_mean, M_planet])
sigms = np.array([sk2_month, sQ_month, sk2Q_year, sk3_month, sh2_month, sMoIF_mean, sM_planet])

# Number of variables
ndim = 9

# Number of walkers + initialisation
nwalkers = 32
initial  = [80, 21, 0.2, 0, -1, -5, 3.5, 5, 330]
theta = initial

p0 = initial + 1e-2 *np.random.rand(nwalkers, ndim)

# Autocorrelation time
index = 0
autocorr = np.empty(max_n)
old_tau = np.inf

with Pool() as pool:
    # Define the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=[means, sigms], pool=pool)

    # Production run
    for sample in sampler.sample(p0, iterations=max_n, progress=True):
    # Only check convergence every 1000 steps
    
        if sampler.iteration > 1000 and not sampler.iteration % 100:
            f = open("chain_core_REAL.dat", "a")
            #f.write("{0} ".format(sampler.iteration))
            for theta in sample[0]:
                for k in range(len(theta)):
                    f.write("{0:4f} ".format(theta[k]))
                    #f.write("{0:4f} {1:4f} {2:4f} {3:4f}".format(Qmmod, Qamod, k_Love_m.real, k_Love_a.real))
                f.write("\n")
            f.close()
    
        if sampler.iteration % 1000:
            continue

        # Compute the autocorrelation time so far
        tau = sampler.get_autocorr_time(tol=0)
        print(tau)
        autocorr[index] = np.mean(tau)
        index += 1

        # Check convergence
        converged = np.all(tau * 100 < sampler.iteration)
        converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
        if converged:
            break
        old_tau = tau

print(
    "Mean acceptance fraction: {0:.3f}".format(
        np.mean(sampler.acceptance_fraction)
    )
)

tau = np.mean(sampler.get_autocorr_time())
print(
    "Mean autocorrelation time: {0:.3f} steps".format(
        tau
    )
)

# Plot posterior probability density
samples = sampler.get_chain(discard=2*int(tau), thin=2000, flat=True) # get all samples
labels = [r"$\mu_{\rm{m}}$", r"$\log\;\eta_{\rm{m}}$", r"$\alpha$", r"$\zeta$", r"$\log\;\Delta$", 
          r"$\log\;t_{\rm{rel}}$", r"$\rho_{\rm{m}}$", r"$\rho_{\rm{c}}$", r"$R_{\rm{c}}$"]

fig = corner.corner(
    samples, labels=labels
)

plt.savefig("corner.png", dpi=300)

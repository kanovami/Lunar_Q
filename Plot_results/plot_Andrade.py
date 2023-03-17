# -*- coding: utf-8 -*-

import numpy as np
import corner
import scipy.constants as sc
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('font',family='Times New Roman',size=14)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'
from tidal_deformation import mrheo,mconst

#===============================================
# Functions calculating the geodetic parameters
#===============================================
# Total mass:
# -----------
def mass(rho_mantle, rho_core, rho_crust, r_core, d_crust):
    r_mantle = r_planet-d_crust
    M = r_core**3*(rho_core-rho_mantle) + r_mantle**3*(rho_mantle-rho_crust) + r_planet**3*rho_crust
    return 4/3*sc.pi*M

# Moment of inertia factor:
# -------------------------
def MoIF(rho_mantle, rho_core, rho_crust, r_core, d_crust):
    r_mantle = r_planet-d_crust
    moif = r_core**5*(rho_core-rho_mantle) + r_mantle**5*(rho_mantle-rho_crust) + r_planet**5*rho_crust
    moif /= r_core**3*(rho_core-rho_mantle) + r_mantle**3*(rho_mantle-rho_crust) + r_planet**3*rho_crust
    moif *= 2/5/r_planet**2
    return moif

# Tidal Love numbers:
# -------------------
def k2(mu, eta, alpha, zeta, rho_mantle, rho_core, rho_crust, 
       r_core, d_crust, omg):
    
    # Allocate Fortran arrays
    # Planet interior structure - layerwise, from core to surface
    #-------------------------------------------------------------
    mconst.ire  = [0, 1, 0]               # rheological models
    mconst.rho  = [rho_core, rho_mantle, rho_crust]         # densities
    mconst.arad = [r_core, r_planet-d_crust, r_planet]# radii of interfaces
    mconst.eta  = [1e0, eta, 0]       # viscosities
    mconst.mu   = [1e-10, mu, mu]      # rigidities
    mconst.alpha= [None, alpha, None]      # parameter of the Andrade model
    mconst.zeta = [None, zeta, None]         # parameter of the Andrade model
    mconst.eta_kv_rel = [None, None, None] # parameter of the Sundberg-Cooper model
    mconst.cdef = [None, None, None]     # parameter of the Sundberg-Cooper model
    #-------------------------------------------------------------

    mconst.nl = len(mconst.arad)
    mconst.iellid = 1

    k2, h2 = mrheo.peltier(2, omg)
    k3, h3 = mrheo.peltier(3, omg)
    return k2,h2,k3

#==========================================================
# Probability and likelihood calculation
#==========================================================
# Logarithm of probability density function:
# ------------------------------------------
def log_probability(theta, xobs, sig):
    '''
    Logarithm of probability density function
    '''  
    k_Love_m, h_Love_m, k_Love3_m = k2(
                theta[0]*1e9, 10**theta[1], theta[2], 10**theta[3],
                theta[4]*1e3, theta[5]*1e3,
                2550, theta[6]*1e3, 40e3, omg_month)
    k_Love_a, h_Love_a, k_Love3_a = k2(
                theta[0]*1e9, 10**theta[1], theta[2], 10**theta[3],
                theta[4]*1e3, theta[5]*1e3,
                2550, theta[6]*1e3, 40e3, omg_year)
    
    k2mod = k_Love_m.real
    k3mod = k_Love3_m.real
    h2mod = h_Love_m.real
    Qmmod = 1/np.sin(abs(np.angle(k_Love_m)))
    #Qamod = 1/np.sin(abs(np.angle(k_Love_a)))
    k2Q_a = -k_Love_a.imag
    MoIFmod = MoIF(theta[4]*1e3, theta[5]*1e3, 2550, theta[6]*1e3, 40e3)
    massmod = mass(theta[4]*1e3, theta[5]*1e3, 2550, theta[6]*1e3, 40e3)
    
    xmod = np.array([k2mod, Qmmod, k2Q_a, k3mod, h2mod, MoIFmod, massmod])
   
    a  = (xmod-xobs)
    b  = np.diag(1/sig**2)
    S  = np.matmul(np.matmul(a,b), np.transpose(a))

    return -0.5 * S, xmod, S

#===============================
# Observed values
#===============================
# monthly tide
Q_month   = 38      # Williams & Boggs (2015)
sQ_month  = 4
k2_month  = 0.02422 # Williams et al. (2014)
sk2_month = 0.00022
h2_month  = 0.0387  # Thor et al. (2021)
sh2_month = 0.0025
k3_month  = 0.0081  # unweighted mean, Konopliv+13, Lemoine+13
sk3_month = 0.0018
# annual tide
k2Q_year  = 6.2e-4  # Williams & Boggs (2015)
sk2Q_year = 1.4e-4
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

#--------------------------------
# Read from file
#--------------------------------
burnin = 1165

# Define file name
file_name = './Andrade/chain_core_REAL.dat'

df = np.loadtxt(file_name)
#print(np.shape(df))
length = len(df)

#==========================
# Find best-fitting sample
#==========================
best_vals = df[0]
best_vals_10 = np.zeros((10, len(df[0])))
best_fits_10 = np.full(10, -np.inf)
best_chi2_10 = np.zeros(10)
observables_10 = np.zeros((10, len(means)))
chi2 = []

for sample in df[burnin:,:]:
    logprob = log_probability(sample, means, sigms)
    chi2.append(logprob[2])
    if min(best_fits_10)<logprob[0]: 
        indx_min = np.argmin(best_fits_10)
        best_fits_10[indx_min] = logprob[0]
        best_chi2_10[indx_min] = logprob[2]
        best_vals_10[indx_min] = sample
        observables_10[indx_min] = logprob[1]
        
chi2 = np.array(chi2)
        
# print(best_vals_10)
# print(best_fits_10)
# print(best_chi2_10)

df[:,4] *= 1000
df[:,5] *= 1000

best_vals_10[:,4] *= 1000
best_vals_10[:,5] *= 1000

indx_min = np.argmin(best_chi2_10)
best_vals = best_vals_10[indx_min]

with open('Andrade-best_fits.csv', 'w', newline='') as f:
    writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
    writer.writerow(["chi2","mu","eta","alpha","zeta",
                     "rho_mantle","rho_core","R_core"])
    for i in range(10):
        writer.writerow(np.concatenate([[best_chi2_10[i]],best_vals_10[i]]))

#=============================
# Plot rheological parameters
#=============================
figure = corner.corner(df[burnin::2,:4],
    labels = [r"$\mu_{\rm{m}}\; \left[\rm{GPa}\right]$", 
              r"$\log\;\eta_{\rm{m}}\; \left[\rm{Pa\; s}\right]$", 
              r"$\alpha$", r"$\log\;\zeta$"], 
    smooth=True,
    color='indianred',
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={"fontsize": 12, "pad": 10},
    hist_kwargs={"density": True})

corner.overplot_lines(figure, best_vals[:4], color='black', ls="dotted")
corner.overplot_points(figure, best_vals[:4][None], marker="s", color="black")
corner.overplot_points(figure, best_vals_10[:,:4], marker="s", fillstyle="none", markeredgecolor="black")

plt.savefig('Andrade-RHEO.png', dpi=300, bbox_inches = 'tight')

#==============================
# Plot densities and core size
#==============================
figure = corner.corner(df[burnin::2,4:],
    labels = [r"$\rho_{\rm{m}}\; \left[\rm{kg\; m^{-3}}\right]$", 
              r"$\rho_{\rm{c}}\; \left[\rm{kg\; m^{-3}}\right]$", 
              r"$R_{\rm{c}}\; \left[\rm{km}\right]$"],
    smooth=True,
    color='indianred',
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={"fontsize": 12, "pad": 10},
    hist_kwargs={"density": True})

corner.overplot_lines(figure, best_vals[4:], color='black', ls='dotted')
corner.overplot_points(figure, best_vals[4:][None], marker="s", color="black")
corner.overplot_points(figure, best_vals_10[:,4:], marker="s", fillstyle="none", markeredgecolor="black")

plt.savefig('Andrade-GEO.png', dpi=300, bbox_inches = 'tight')

#=============================================
# Plot k2 and k2/Q as a function of frequency
#=============================================
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax.flatten()
plt.tight_layout()

indices = np.random.choice(length-burnin, 100)

samples  = df[indices]
chi2_100 = chi2[indices]

with open('Andrade-random100.csv', 'w', newline='') as f:
    writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
    writer.writerow(["chi2","mu","eta","alpha","zeta",
                     "rho_mantle","rho_core","R_core"])
    for i in range(100):
        writer.writerow(np.concatenate([[chi2_100[i]], samples[i]]))

# np.random.shuffle(df[burnin:,:])
# samples = df[burnin:burnin+100,:]

# Frequency limits
FREQ_MIN = -8
FREQ_MAX = -4

freqs = np.logspace(FREQ_MIN, FREQ_MAX, 100)

for theta in samples:
    k2s = []
    for omg in freqs:
        k_Love, h_Love, k_Love3 = k2(
            theta[0]*1e9, 10**theta[1], theta[2], 10**theta[3],
            theta[4], theta[5],
            2550, theta[6]*1e3, 40e3, omg)
        k2s.append(k_Love)
        
    k2s = np.array(k2s)

    ax[0].plot(np.log10(freqs), k2s.real,
                            lw=0.2, color='lightskyblue', ls='-')
    ax[1].plot(np.log10(freqs), -k2s.imag*10**3,
                            lw=0.2, color='lightskyblue', ls='-')
    
for theta in best_vals_10:
    #print(theta)
    k2s = []
    for omg in freqs:
        k_Love, h_Love, k_Love3 = k2(
            theta[0]*1e9, 10**theta[1], theta[2], 10**theta[3],
            theta[4], theta[5],
            2550, theta[6]*1e3, 40e3, omg)
        k2s.append(k_Love)
        
    k2s = np.array(k2s)

    ax[0].plot(np.log10(freqs), k2s.real,
                            lw=1, color='cornflowerblue', ls='-')
    ax[1].plot(np.log10(freqs), -k2s.imag*10**3,
                            lw=1, color='cornflowerblue', ls='-')
    
ax[0].set_xlabel(r"$\log\,\chi$ [rad/s]")
ax[0].set_ylabel(r"$\Re\;\left[\bar{k}_2\; (\,\chi)\right]$")
ax[0].hlines(k2_month, FREQ_MIN, FREQ_MAX, color='r')
ax[0].vlines(np.log10(omg_month), 0.022, 0.028, color='r')
ax[0].text(-5.5, 0.0223, r'1 month', color='r')
ax[0].set_xlim(FREQ_MIN, FREQ_MAX)
ax[0].set_ylim(0.022, 0.028)
ax[0].grid(c="gray", ls="dotted")

ax[1].set_xlabel(r"$\log\,\chi$ [rad/s]")
ax[1].set_ylabel(r"$-\Im\;\left[\bar{k}_2\; (\,\chi)\right]\; \left(\times10^{3}\right)$")
ax[1].hlines(k2_month/Q_month*10**3, FREQ_MIN, FREQ_MAX, color='r')
ax[1].hlines(k2Q_year*10**3, FREQ_MIN, FREQ_MAX, color='orange')
ax[1].vlines(np.log10(omg_month), 0, 0.002*10**3, color='r')
ax[1].vlines(np.log10(omg_year), 0, 0.002*10**3, color='orange')
ax[1].text(-6.6, 0.0001*10**3, r'1 yr', color='orange')
ax[1].text(-5.5, 0.0001*10**3, r'1 month', color='r')
ax[1].set_xlim(FREQ_MIN, FREQ_MAX)
ax[1].set_ylim(0, 0.002*10**3)
ax[1].grid(c="gray", ls="dotted")

fig.subplots_adjust(wspace=0.22)

plt.savefig('./Andrade-OVERVIEW.png', dpi=200, bbox_inches = 'tight')
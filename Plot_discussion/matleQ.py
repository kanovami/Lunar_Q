# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
mpl.rc('font',family='Times New Roman',size=14)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

import matplotlib.pyplot as plt
import scipy.constants as sc
from tidal_deformation import mrheo,mconst

def k2_calc(mu, eta, alpha, zeta, Delta, trel, rho_mantle, rho_core, rho_crust, 
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
    mconst.iellid = 0

    k2, h2, Qinv = mrheo.peltier(2, omg)
    return k2,h2,Qinv

#---------------------------------------------------------
# Input data
#---------------------------------------------------------
# Define file name
file_name = '../Sundberg-Cooper/chain_core_REAL.dat'

# Read from file
df = np.loadtxt(file_name)

burnin = 8896

df[:,6] *= 1000
df[:,7] *= 1000

# Williams and Boggs (2015):
freq_month = 2*sc.pi/(27.212*86400) # monthly
freq_year  = 2*sc.pi/(365.260*86400) # annual

r_planet = 1737.151e3 # Wieczorek (2015)

# Frequency limits
FREQ_MIN = np.log10(freq_year) #-8
FREQ_MAX = -4

freqs = np.logspace(FREQ_MIN, FREQ_MAX, 100)

fig, ax = plt.subplots(1, 1, figsize=(6.47, 4))
fig.tight_layout()

np.random.shuffle(df[burnin::2,:])
samples = df[burnin:burnin+10000,:]

Qmaxs = []
for theta in samples:
    Qinvs = []
    for omg in freqs:
        k_Love, h_Love, Qinv = k2_calc(
            theta[0]*1e9, 10**theta[1], theta[2], 10**theta[3],
            10**theta[4], 10**theta[5], theta[6], theta[7],
            2550, theta[8]*1e3, 40e3, omg)
        Qinvs.append(Qinv)

    Qmax=max(Qinvs)
    Qmaxs.append(Qmax)
    
Qmaxs = np.array(Qmaxs)
ax.scatter(samples[:,4], Qmaxs, s=2, c='dodgerblue')
ax.set_xlabel(r"$\log\;\Delta$")
ax.set_ylabel(r"seismic $Q^{-1}$ in the secondary peak")
ax.set_xlim(-2.1, -0.25)
ax.set_yticks(np.arange(0, 0.25, 0.05))
ax.grid(c="gray", ls="dotted")

plt.savefig('./Delta_Q.png', dpi=300, bbox_inches = 'tight')
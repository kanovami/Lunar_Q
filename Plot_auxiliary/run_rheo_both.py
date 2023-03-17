#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 11:22:17 2021

@author: Michaela Walterova
"""
import numpy as np
import matplotlib as mpl
mpl.rc('font',family='Times New Roman',size=26)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.constants as sc
from tidal_deformation import mrheo,mconst

# Figure parameters
fig, ax = plt.subplots(1, 2, figsize=(16.5, 6.5))
fig.tight_layout()
ax = ax.flatten()
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)

# Overwrite predefined constants by Scipy values
mconst.kappa = sc.G
mconst.pi    = sc.pi

NUM = 8

# Define colors
colors = []
cmap = plt.cm.get_cmap('Blues')
for c in np.linspace(0.05,1,NUM):
    colors.append(cmap(c))

# Williams and Boggs (2015):
freq_month = 2*sc.pi/(27.212*86400) # monthly
freq_year  = 2*sc.pi/(365.260*86400) # annual

# Frequency limits
FREQ_MIN = -15
FREQ_MAX = 0
freqs = np.logspace(FREQ_MIN, FREQ_MAX, 100)

# Interior structure - parameters
#---------------------------------
ETA_MELT = 16
ETA_MANTLE = 20

ire  = 2

rho_core = 5e3
rho_mantle = 3.3e3

r_core = 330e3
r_planet = 1737e3
d_melt = 500e3

cdef_mantle = 0.2 # relaxation strength
cdef_core = cdef_mantle
#---------------------------------

# Allocate Fortran arrays
# Planet interior structure - layerwise, from core to surface
#-------------------------------------------------------------
mconst.ire  = [0, 0, ire]               # rheological models
mconst.rho  = [rho_core, rho_mantle, rho_mantle]         # densities
mconst.arad = [r_core, r_core+d_melt, r_planet]# radii of interfaces
mconst.eta  = [1e0, 10**ETA_MELT, 10**ETA_MANTLE]       # viscosities
mconst.mu   = [1e-10, 1e10, 8e10]      # rigidities
mconst.alpha= [0.3, 0.3, 0.3]      # parameter of the Andrade model
mconst.zeta = [1, 1, 1]         # parameter of the Andrade model
mconst.cdef = [cdef_core, cdef_core, cdef_mantle]     # parameter of the Sundberg-Cooper model
#-------------------------------------------------------------

mconst.nl = len(mconst.arad)
mconst.iellid = 0

i = 0
for trel in np.logspace(0,-7,NUM):
  etakv_mantle = trel/cdef_mantle
  
  k2s  = []
  epss = []

  for freq in freqs:
    mconst.eta_kv_rel = [1e0, 1e0, etakv_mantle] # parameter of the Sundberg-Cooper model
    
    klove, hlove = mrheo.peltier(2,freq)
    k2s.append(abs(klove))
    epss.append(-np.angle(klove))
    
  color = colors[i]
    
  i+=1
    
  ax[0].plot(np.log10(freqs), np.log10(k2s*np.sin(epss)), 
             label=str(trel), lw=6, c=color)
  ax[0].set_xlabel(r"$\log\,\chi$ [rad/s]")
  ax[0].set_ylabel(r"$\log\;k_2/Q$")
  ax[0].set_xlim(FREQ_MIN, FREQ_MAX)
  ax[0].set_ylim(-5, 1)
  ax[0].set_xticks(np.arange(FREQ_MIN, FREQ_MAX+2.5, 2.5))
  
  ax[1].plot(np.log10(freqs), np.log10(np.sin(epss)), 
             label=str(trel), lw=6, c=color)
  ax[1].legend(fontsize=24, ncol=2)
  ax[1].set_xlabel(r"$\log\,\chi$ [rad/s]")
  ax[1].set_ylabel(r"$\log\;Q^{-1}$")
  ax[1].set_xlim(FREQ_MIN, FREQ_MAX)
  ax[1].set_ylim(-5, 1)
  
ax[0].grid(c="gray", ls="dotted")
ax[1].grid(c="gray", ls="dotted")

fig.subplots_adjust(wspace=0.2)

plt.savefig('Moon_trel_diff.png', dpi=300, bbox_inches = 'tight')
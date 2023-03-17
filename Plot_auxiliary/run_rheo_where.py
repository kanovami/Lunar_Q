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

# Overwrite predefined constants by Scipy values
mconst.kappa = sc.G
mconst.pi    = sc.pi

NUM = 8

# Define colors
colors = []
cmap = plt.cm.get_cmap('Blues')
for c in np.linspace(0.05,1,NUM):
    colors.append(cmap(c))

fig, ax = plt.subplots(1, 2, figsize=(16.5, 6.5))
fig.tight_layout()
ax = ax.flatten()
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)

# Williams and Boggs (2015):
freq_month = 2*sc.pi/(27.212*86400) # monthly
freq_year  = 2*sc.pi/(365.260*86400) # annual

# Frequency limits
FREQ_MIN = -10
FREQ_MAX = -4
freqs = np.logspace(FREQ_MIN, FREQ_MAX, 100)

# Interior structure - parameters
#---------------------------------
ETA_MANTLE = 22
MU_MANTLE  = 8e10
ire  = 2

rho_core = 5e3
rho_mantle = 3.3e3

alpha = 0.25
zeta  = 1

r_core = 330e3
r_planet = 1737e3

cdef_mantle = 0.05 # relaxation strength
cdef_core = cdef_mantle
#---------------------------------

# Allocate Fortran arrays
# Planet interior structure - layerwise, from core to surface
#-------------------------------------------------------------
mconst.ire  = [ire, ire]               # rheological models
mconst.rho  = [rho_mantle, rho_mantle]         # densities
mconst.arad = [r_core, r_planet]# radii of interfaces
mconst.eta  = [10**ETA_MANTLE, 10**ETA_MANTLE]       # viscosities
mconst.mu   = [MU_MANTLE, MU_MANTLE]      # rigidities
mconst.alpha= [alpha, alpha]      # parameter of the Andrade model
mconst.zeta = [zeta, zeta]         # parameter of the Andrade model
mconst.cdef = [cdef_core, cdef_mantle]     # parameter of the Sundberg-Cooper model
#-------------------------------------------------------------

mconst.nl = len(mconst.arad)
mconst.iellid = 0

i = 0
for trel in np.logspace(0,-7,NUM):
  etakv_mantle = trel/cdef_mantle
  k2s  = []
  epss = []
  
  mconst.eta_kv_rel = [etakv_mantle, etakv_mantle] # parameter of the Sundberg-Cooper model

  for freq in freqs:
    klove,hlove = mrheo.peltier(2,freq)
    k2s.append(abs(klove))
    epss.append(-np.angle(klove))
       
  color = colors[i]
  i+=1
    
  ax[0].plot(np.log10(freqs), np.log10(k2s*np.sin(epss)), label=str(trel), lw=6, c=color)

  ax[0].set_xlabel(r"$\log\,\chi$ [rad/s]")
  ax[0].set_ylabel(r"$\log\;k_2/Q$")
  ax[0].set_xlim(FREQ_MIN, FREQ_MAX)
  ax[0].set_ylim(-4.5, -2)
  ax[0].set_xticks(np.arange(FREQ_MIN, FREQ_MAX+1, 1))
  ax[0].legend(fontsize=20, ncol=2, loc=3)

  ax[1].plot(np.log10(freqs), np.log10(np.sin(epss)), label=str(trel), lw=6, c=color)
  ax[1].set_xlabel(r"$\log\,\chi$ [rad/s]")
  ax[1].set_ylabel(r"$\log\;Q^{-1}$")
  ax[1].set_xlim(FREQ_MIN, FREQ_MAX)
  ax[1].set_ylim(-2.5, -0.6)
  
  ax[1].hlines(np.log10(1./38), FREQ_MIN, FREQ_MAX, color='r')
  ax[1].hlines(np.log10(1./41), FREQ_MIN, FREQ_MAX, color='orange')
  ax[1].vlines(np.log10(freq_month), -4, -0.5, color='r')
  ax[1].vlines(np.log10(freq_year), -4, -0.5, color='orange')
  ax[1].text(-6.6, -2.4, r'1 yr', fontsize=24, color='orange')
  ax[1].text(-5.5, -2.4, r'1 month', fontsize=24, color='r')
  ax[1].set_yticks(np.arange(-2.5, 0, 0.5))
  
ax[0].grid(c="gray", ls="dotted")
ax[1].grid(c="gray", ls="dotted")

fig.subplots_adjust(wspace=0.23)

plt.savefig('Moon_trel_where.png', dpi=300, bbox_inches = 'tight')
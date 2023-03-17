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
ETA_MELT = 15
ETA_MANTLE = 20

rho_core = 5e3
rho_mantle = 3.3e3

r_core = 330e3
r_planet = 1737e3
d_melt = 500e3

cdef_mantle = 0.15 # relaxation strength
cdef_core = cdef_mantle

trel = 0.0001
etakv_mantle = trel/cdef_mantle

Dsc = 0e3 # depth at which the planet starts to follow the Sundberg-Cooper model
#---------------------------------

k2s_lay  = []
epss_lay = []
k2s_sc  = []
epss_sc = []
k2s_a  = []
epss_a = []
  
for freq in freqs:
   
    # LAYERED - ANDRADE
    # Allocate Fortran arrays
    # Planet interior structure - layerwise, from core to surface
    #-------------------------------------------------------------
    mconst.ire  = [0, 0, 1]               # rheological models
    mconst.rho  = [rho_core, rho_mantle, rho_mantle]         # densities
    mconst.arad = [r_core, r_core+d_melt, r_planet]# radii of interfaces
    mconst.eta  = [1e0, 10**ETA_MELT, 10**(ETA_MANTLE+0.2)]       # viscosities
    mconst.mu   = [1e-10, 1e10, 8e10]      # rigidities
    mconst.alpha= [0.3, 0.3, 0.3]      # parameter of the Andrade model
    mconst.zeta = [1, 1, 1]         # parameter of the Andrade model
    mconst.eta_kv_rel = [1e0, 1e0, etakv_mantle] # parameter of the Sundberg-Cooper model
    mconst.cdef = [cdef_core, cdef_core, cdef_mantle]     # parameter of the Sundberg-Cooper model
    #-------------------------------------------------------------
 
    mconst.nl = len(mconst.arad)
    mconst.iellid = 0

    klove, hlove = mrheo.peltier(2,freq)
    k2s_lay.append(abs(klove))
    epss_lay.append(-np.angle(klove))
    
    label_lay = "Andrade, layered"
    color_lay = "red"
    
    # HOMOGENEOUS, ANDRDAE
    # Allocate Fortran arrays
    # Planet interior structure - layerwise, from core to surface
    #-------------------------------------------------------------
    mconst.ire  = [1]               # rheological models
    mconst.rho  = [rho_mantle]         # densities
    mconst.arad = [r_planet]# radii of interfaces
    mconst.eta  = [10**ETA_MANTLE]       # viscosities
    mconst.mu   = [8e10]      # rigidities
    mconst.alpha= [0.3]      # parameter of the Andrade model
    mconst.zeta = [1]         # parameter of the Andrade model
    mconst.eta_kv_rel = [etakv_mantle] # parameter of the Sundberg-Cooper model
    mconst.cdef = [cdef_mantle]     # parameter of the Sundberg-Cooper model
    #-------------------------------------------------------------

    mconst.nl = len(mconst.arad)
    mconst.iellid = 0

    klove, hlove = mrheo.peltier(2,freq)
    k2s_a.append(abs(klove))
    epss_a.append(-np.angle(klove))
    
    label_a = "Andrade, homogeneous"
    color_a = "red"
    
    # HOMOGENEOUS, SUNDBERG-COOPER
    # Allocate Fortran arrays
    # Planet interior structure - layerwise, from core to surface
    #-------------------------------------------------------------
    mconst.ire  = [2, 1]               # rheological models
    mconst.rho  = [rho_mantle, rho_mantle]         # densities
    mconst.arad = [r_planet-Dsc, r_planet]# radii of interfaces
    mconst.eta  = [10**ETA_MANTLE, 10**ETA_MANTLE]       # viscosities
    mconst.mu   = [8e10, 8e10]      # rigidities
    mconst.alpha= [0.3, 0.3]      # parameter of the Andrade model
    mconst.zeta = [1, 1]         # parameter of the Andrade model
    mconst.eta_kv_rel = [etakv_mantle, etakv_mantle] # parameter of the Sundberg-Cooper model
    mconst.cdef = [cdef_mantle, cdef_mantle]     # parameter of the Sundberg-Cooper model
    #-------------------------------------------------------------

    mconst.nl = len(mconst.arad)
    mconst.iellid = 0

    klove, hlove = mrheo.peltier(2,freq)
    k2s_sc.append(abs(klove))
    epss_sc.append(-np.angle(klove))
    
    label_sc = "Sundberg-Cooper"
    color_sc = "dodgerblue"

ax[0].plot(np.log10(freqs), np.log10(k2s_a*np.sin(epss_a)), label=label_a, lw=4, ls='--', c=color_a)
ax[0].plot(np.log10(freqs), np.log10(k2s_lay*np.sin(epss_lay)), label=label_lay, lw=6, c=color_lay)    
ax[0].plot(np.log10(freqs), np.log10(k2s_sc*np.sin(epss_sc)), label=label_sc, lw=6, c=color_sc)

#ax[0].legend(fontsize=24)
ax[0].set_xlabel(r"$\log\,\chi$ [rad/s]")
ax[0].set_ylabel(r"$\log\;k_{2}/Q_{2}$")
ax[0].set_xlim(FREQ_MIN, FREQ_MAX)
ax[0].set_ylim(-5, 1)
ax[0].set_xticks(np.arange(FREQ_MIN, FREQ_MAX+2.5, 2.5))
  
ax[1].plot(np.log10(freqs), np.log10(np.sin(epss_a)), label=label_a, lw=4, ls='--', c=color_a)
ax[1].plot(np.log10(freqs), np.log10(np.sin(epss_lay)), label=label_lay, lw=6, c=color_lay)
ax[1].plot(np.log10(freqs), np.log10(np.sin(epss_sc)), label=label_sc, lw=6, c=color_sc)
ax[1].legend(fontsize=24, ncol=1, loc=8)
ax[1].set_xlabel(r"$\log\,\chi$ [rad/s]")
ax[1].set_ylabel(r"$\log\;Q_{2}^{-1}$")
ax[1].set_xlim(FREQ_MIN, FREQ_MAX)
ax[1].set_ylim(-5, 1)
  
ax[0].grid(c="gray", ls="dotted")
ax[1].grid(c="gray", ls="dotted")

fig.subplots_adjust(wspace=0.2)

plt.savefig('sc_vs_layered.png', dpi=300, bbox_inches = 'tight')
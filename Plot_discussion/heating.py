# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import matplotlib as mpl
mpl.rc('font',family='Times New Roman',size=18)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

N = 100
fig, ax = plt.subplots(2, 2, figsize=(14, 7))
fig.tight_layout()
ax = ax.flatten()

BOTTOM = 0.1
TOP = 0.95
WSPACE = 0.3
HSPACE = 0.3
RIGHT = 0.82

nlev  = 12
nlev2 = 30

pos_x = [0, np.pi/2, np.pi, 3*np.pi/2, 6.248278722]
lab_x = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]
pos_y = [0, np.pi/2, np.pi]
lab_y = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"]
pos_xx = [0.4363323130e-01, np.pi/4, np.pi/2, 3*np.pi/4, 3.097959422]
lab_xx = [r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$", r"$\frac{3\pi}{4}$", r"$\pi$"]

cmF = 'plasma'
cm = 'plasma'

#=======================================
# Subplot 1
#=======================================
filename = "./Tidal_heating/Surface_heat_flux/Moon_tidflux_SC.dat"
x, y, z = np.genfromtxt(filename, unpack=True)
filename = "./Tidal_heating/Surface_heat_flux/Moon_tidflux_MELT.dat"
x2, y2, z2 = np.genfromtxt(filename, unpack=True)

MIN = min(z.min(), z2.min())
MAX = max(z.max(), z2.max())

xi = np.linspace(x.min(), x.max(), N)
yi = np.linspace(y.min(), y.max(), N)
zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

ax[0].contourf(xi, yi, zi, np.arange(MIN-(MAX-MIN)/nlev, MAX+(MAX-MIN)/nlev, (MAX-MIN)/nlev), cmap=cmF)
ax[0].set_xlabel(r"$\varphi$")
ax[0].set_ylabel(r"$\vartheta$")
ax[0].set_xticks(pos_x)
ax[0].set_xticklabels(lab_x)
ax[0].set_yticks(pos_y)
ax[0].set_yticklabels(lab_y)
ax[0].set_title("Sundberg-Cooper mantle", fontweight='bold', pad=20)

#=======================================
# Subplot 2
#=======================================
xi = np.linspace(x2.min(), x2.max(), N)
yi = np.linspace(y2.min(), y2.max(), N)
zi = scipy.interpolate.griddata((x2, y2), z2, (xi[None,:], yi[:,None]), method='cubic')

cs1 = ax[1].contourf(xi, yi, zi, np.arange(MIN-(MAX-MIN)/nlev, MAX+(MAX-MIN)/nlev, (MAX-MIN)/nlev), cmap=cmF)
ax[1].set_xlabel(r"$\varphi$")
ax[1].set_ylabel(r"$\vartheta$")
ax[1].set_xticks(pos_x)
ax[1].set_xticklabels(lab_x)
ax[1].set_yticks(pos_y)
ax[1].set_yticklabels(lab_y)

cbar_ax = fig.add_axes([0.85, TOP-0.435*(TOP-BOTTOM), 0.02, 0.435*(TOP-BOTTOM)])
cbar = fig.colorbar(cs1, cax=cbar_ax)
cbar.set_label(r"$\Phi_{\rm{tide}}\; \left[{\rm mW/m}^2\right]$")
cbar.set_ticks(np.arange(0.011, 0.021, 0.002))

ax[1].set_title("Mantle with a low-viscosity layer", fontweight='bold', pad=20)

#=======================================
# Subplot 3
#=======================================
filename = "./Tidal_heating/Volumetric_heating/Moon_sc_rad.dat"
y, x, z = np.genfromtxt(filename, unpack=True)
filename = "./Tidal_heating/Volumetric_heating/Moon_melt_rad.dat"
y2, x2, z2 = np.genfromtxt(filename, unpack=True)

MIN = min(np.log10(z.min()), np.log10(z2.min()))
MAX = max(np.log10(z.max()), np.log10(z2.max()))

xi = np.linspace(x.min(), x.max(), N)
yi = np.linspace(y.min(), y.max(), N)
zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

ax[2].contourf(xi, yi, np.log10(zi), np.arange(MIN-(MAX-MIN)/nlev2, MAX+(MAX-MIN)/nlev2, (MAX-MIN)/nlev2), cmap=cm)
ax[2].set_xlabel(r"$\vartheta$")
ax[2].set_ylabel(r"$r/R$")
ax[2].set_xticks(pos_xx)
ax[2].set_xticklabels(lab_xx)

#=======================================
# Subplot 4
#=======================================
xi = np.linspace(x2.min(), x2.max(), N)
yi = np.linspace(y2.min(), y2.max(), N)
zi = scipy.interpolate.griddata((x2, y2), z2, (xi[None,:], yi[:,None]), method='cubic')

cs2 = ax[3].contourf(xi, yi, np.log10(zi), np.arange(MIN-(MAX-MIN)/nlev2, MAX+(MAX-MIN)/nlev2, (MAX-MIN)/nlev2), cmap=cm)
ax[3].set_xlabel(r"$\vartheta$")
ax[3].set_ylabel(r"$r/R$")
ax[3].set_xticks(pos_xx)
ax[3].set_xticklabels(lab_xx)

cbar_ax = fig.add_axes([0.85, BOTTOM, 0.02, 0.435*(TOP-BOTTOM)])
cbar = fig.colorbar(cs2, cax=cbar_ax)
cbar.set_label(r"$\log\;h_{\rm{tide}}\; \left[{\rm W/m}^3\right]$")
cbar.set_ticks(np.arange(-12, -7, 1))

fig.subplots_adjust(bottom=BOTTOM, wspace=WSPACE, hspace=HSPACE, right=RIGHT, top=TOP)

plt.savefig('./heating_map.png', dpi=300, bbox_inches = 'tight')
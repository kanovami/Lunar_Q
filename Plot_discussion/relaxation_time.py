# -*- coding: utf-8 -*-
# Written by Michaela Walterová, 2022

import numpy as np
import matplotlib as mpl
mpl.rc('font',family='Times New Roman',size=20)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'
import matplotlib.pyplot as plt
import scipy.constants as sc

tau_r = 10**1.15  # in s
T_r   = 1173.15   # in K
P_r   = 200e6     # in Pa
d_r   = 1e-5      # in m
V_act = 1e-5      # in m^3/mol
E_act = 259e3     # in J/mol

Rp = 1737.151e3
Rc = 441.85e3
rhom = 3.373e3
rhoc = 5.054e3
Rd = Rp-40e3

cmin = 0
cmax = 20
nlev = 10

def tau(d=d_r, P=P_r, T=T_r):
    m = 1.31
    return tau_r*(d/d_r)**m*np.exp(E_act/sc.R*(1/T-1/T_r))*np.exp(V_act/sc.R*(P/T-P_r/T_r))

def rev_tau(d=d_r, P=P_r, tau=10):
    m = 1.31
    return (E_act+P*V_act)/sc.R/(np.log(tau/tau_r)-m*np.log(d/d_r)+(E_act+P_r*V_act)/sc.R/T_r)

def cond_T(r, H2=9.5e-9):
    Ts = 250
    H1 = 160e-9
    H2 = H2 #9.5e-9
    k  = 3
    Fc = 0
    Tcon = (Ts - H2*r**2/6/k + (Fc*Rc**2-H2*Rc**3/3)/k*(1/r-1/Rp) + H1*Rp**2/6/k
            + Rd**3*(H1-H2)/3/k/Rp - Rd**2*(H1-H2)/2/k)
    return Tcon

def gravM(r):
    return 4/3*np.pi*sc.G*(rhoc*Rc**3 + rhom*(r**3-Rc**3))/r**2

def presM(r):
    z = Rp-r
    pres = (-4*np.pi*sc.G*rhom/3 * (rhoc*Rc**3*(-1/(Rp-z)+1/Rp) + rhom*((Rp-z)**2/2 + Rc**3/(Rp-z)
           -Rp**2/2 - Rc**3/Rp)))
    return pres

# for r in np.arange(Rc, Rp, 1e3):
#     print(r/1e3, gravM(r), presM(r)/1e9)

fig, ax = plt.subplots(1, 2, sharex=True, figsize=(14, 5))
fig.tight_layout()
ax.flatten()

grainsizes = np.logspace(-10, -2, 500)
radii = np.linspace(Rc, Rd, 500)
pressures  = presM(radii) #np.linspace(0, 5, 100)*1e9
temperatures = np.linspace(cond_T(Rd), cond_T(Rc), 500)
temps_cond = cond_T(radii)
temps_cond1 = cond_T(radii, 8e-9)
temps_cond2 = cond_T(radii, 11e-9)

#for r, p, t in zip(radii, pressures, temps_cond):
#    print((Rp-r)/1e3, p/1e9, np.log10(tau(d=1e-4, P=p, T=t)))
    
X, Y = np.meshgrid(radii, temperatures)

levels = np.linspace(cmin, cmax, nlev)

#print(rev_tau(d=1e-4, P=presM(X), tau=10**4.25))

cs = ax[0].contourf(X/1e3, Y, np.log10(tau(d=1e-4, P=presM(X), T=Y)), levels=levels, vmin=cmin, vmax=cmax, extend='both', cmap='plasma')
ax[0].plot(radii/1e3, rev_tau(d=1e-4, P=presM(radii), tau=10**4), lw=2, c='w')
ax[0].plot(radii/1e3, rev_tau(d=1e-4, P=presM(radii), tau=10**6), lw=2, c='w')
ax[0].plot(radii/1e3, temps_cond, lw=3, c='dodgerblue')
ax[0].plot(radii/1e3, temps_cond1, lw=3, c='dodgerblue', ls='--')
ax[0].plot(radii/1e3, temps_cond2, lw=3, c='dodgerblue', ls='--')
ax[0].set_ylim(cond_T(Rd), cond_T(Rc))
ax[0].set_xlabel(r"$r$ [km]")
ax[0].set_ylabel(r"$T$ [K]")
ax[0].set_xticks(np.arange(400, 1800, 300))
ax[0].text(0.95, 0.95, r'$d=10^{-4}$ m',
        verticalalignment='top', horizontalalignment='right',
        transform=ax[0].transAxes, bbox={'facecolor':'white', 'pad':10, 'alpha':1})

ax[1].contourf(X/1e3, Y, np.log10(tau(d=1e-2, P=presM(X), T=Y)), levels=levels, vmin=cmin, vmax=cmax, extend='both', cmap='plasma')
ax[1].plot(radii/1e3, rev_tau(d=1e-2, P=presM(radii), tau=10**4), lw=2, c='w')
ax[1].plot(radii/1e3, rev_tau(d=1e-2, P=presM(radii), tau=10**6), lw=2, c='w')
ax[1].plot(radii/1e3, temps_cond, lw=3, c='dodgerblue')
ax[1].plot(radii/1e3, temps_cond1, lw=3, c='dodgerblue', ls='--')
ax[1].plot(radii/1e3, temps_cond2, lw=3, c='dodgerblue', ls='--')
ax[1].set_ylim(cond_T(Rd), cond_T(Rc))
ax[1].set_xlabel(r"$r$ [km]")
ax[1].set_ylabel(r"$T$ [K]")
ax[1].text(0.95, 0.95, r'$d=10^{-2}$ m',
        verticalalignment='top', horizontalalignment='right',
        transform=ax[1].transAxes, bbox={'facecolor':'white', 'pad':10, 'alpha':1})

fig.subplots_adjust(right=0.80, bottom=0.15, top=0.90, wspace=0.30)
cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.75])
cbar = plt.colorbar(cs, cax=cbar_ax, ticks=np.arange(0, 24, 4))
cbar.set_label(r'$\log\;\tau$')

plt.savefig("tau_Jackson_etal_2014_wl.png", dpi=300, bbox_inches = 'tight')
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#main file, plotting figures, discussion on Moon basal layer
#coded Marie Behounkova marie.behounkova@matfyz.cuni.cz

import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d

#in house classes
import basal_melting as melting
import basal_temp as temp
import basal_viscosity as viscosity
import basal_porosity as por
import basal_process_results as process_results

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
   

#plot settings
mpl.rc('font',family='Times New Roman',size=24)
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

mpl.rc('xtick', labelsize=24) 
mpl.rc('ytick', labelsize=24) 

#----------------------------------------------------------------
#----------------------------------------------------------------

#plotting parameters
fs = 24
fs2 = 18
fs3 = 12
fs4 = 16
lw = 6
lw2 = 6
fact=1. #2.54
figsize = (16.5/2.*fact, 6*fact) #one plot
figsize_complex = (16.5*fact, 9*fact) #twoplots
figsize2 = (16.5*fact, 6*fact) #twoplots
figsize3 = (24, 8) #(16.5*fact*1.5, 6*fact*1.5) #twoplots


show_results = False


#------------------------------------------------------------------------
#Khan et al. (2014) a simple digitalization of temperature

khan_digi_depth2 = np.array([236,245,342,441,541,641,741,857])
khan_digi_min2 = np.array([1226,1360,1390,1418,1440,1468,1489,1519])
khan_digi_max2 = np.array([1230,1370,1410,1441,1461,1483,1513,1556])

khan_digi_depth2 = (khan_digi_depth2-236)*(1200)/(808-236)
khan_digi_min2 = (khan_digi_min2-1203)*(1500)/(1492-1203)
khan_digi_max2 = (khan_digi_max2-1203)*(1500)/(1492-1203)

scale_depth = 1200./(814-102.)
scale_temp = 1500./(425-77.)

khan_digi_depth = (np.array([102.,112.,115.,130.,180.,280.,330.,430.,480.,530.,630.,680.,730.,780.,830.,874])-102)*scale_depth
khan_digi_min =   (np.array([ 78.,251.,270.,271.,285.,317.,325.,353.,362.,378.,406.,416.,422.,432.,441.,460])-77.)*scale_temp #in celsius
khan_digi_max =   (np.array([ 78.,268.,280.,280.,308.,343.,353.,378.,387.,397.,423.,437.,451.,460.,467.,504])-77)*scale_temp #in celsius


#-------------------------------------------------------
class Moon(object):
        #object defining parameters of the Moon
        #----------------------------------------------------------------
        def __init__(self,*args, **kwargs):

            self.Rs = 1737e3 #radius in m
            self.Rc = 371.45e3 #radius of core in m
            self.k0 = 3. #thermal conductivity in m

            self.Ts = 250 #surface temperature in K
            self.qbot = 0e-3 #Wm^-2
                        

            #self.Rbasal = self.Rc+70e3 #basal layer radius
            self.Rbasal = 678.98e3 #500e3 #basal layer radius

            self.hcrust = 160e-9 #W/m^3 heating in crust

            self.hav = 8e-9 #average heating W/m^3
            self.Cenrich = 1. #default enrichment factor of basal layer


            self.alpha = 4e-5 #thermal expansion
            self.cp = 1250. #specific heat

            self.rho_crust = 2550.
            self.d_crust = 40e3

            #self.Rd = self.Rs-45e3 #crust radius            
            self.Rd = self.Rs-self.d_crust #crust radius            

            self.update()

        #----------------------------------------------------------------
        #updating structure of moon, including e.g. heating if 
        def update(self,*args,**kwargs):

            V_all = 4./3.*np.pi*(self.Rd**3-self.Rc**3)
            V_mantle = 4./3.*np.pi*(self.Rd**3-self.Rbasal**3)
            V_basal = 4./3.*np.pi*(self.Rbasal**3-self.Rc**3)

            self.hmantle = V_all/(self.Cenrich*V_basal+V_mantle)*self.hav
            self.hbase = self.hmantle*self.Cenrich

            self.r=np.array([self.Rc,self.Rbasal,self.Rd,self.Rs])
            self.k=np.array([self.k0,self.k0,self.k0])
            self.h=np.array([self.hbase,self.hmantle,self.hcrust])

            #ini of pressure
            #first density
            self.ini_pressure()


        #----------------------------------------------------------------
        def ini_pressure(self,*args,**kwargs):
            #computation of pressure and gravity acceleration as a function of as a function of radius
            G=6.6743e-11 #gravitational constant
            name = 'MoonKhan_JGR_2014_mean.nd' #name of the file with density
            data=np.genfromtxt(name)
            g_surf = 1.62
            rho = data[:,3]*1000.
            depth = data[:,0]
            frho = interp1d(depth, rho)
            nlayer= 100
            #only in mantle
            rarr =np.linspace(self.Rs,0.9*self.Rc,nlayer)
            press = np.zeros_like(rarr)
            grav = np.zeros_like(rarr)
            dr = (self.Rs-self.Rc)/nlayer
            press_aux = 0.
            grav[0] = g_surf
            for i in range(1,nlayer):
                raux = rarr[i]
                g_surf = g_surf*(raux)**2/(raux-dr)**2-G*4./3.*np.pi*(raux**3-(raux-dr)**3)*frho((self.Rs-raux)/1000.)/(raux-dr)**2
                press_aux += dr*frho((self.Rs-raux)/1000.)*g_surf
                press[i] = press_aux/1e9 #in GPa
                grav[i] = g_surf

            self.fpress = interp1d((self.Rs-rarr)/1000., press) #interpolation, for getting pressure in any point, input in depth in km
            self.fgrav = interp1d((self.Rs-rarr)/1000., grav) #interpolation, for getting the acceleration in any point, input in depth in km

        #----------------------------------------------------------------
        def adiabat_grad(self,z,T,*args,**kwargs): #for computation of adiabat
            return [self.alpha*self.fgrav(z/1000.)/self.cp*T] #z is depth in m

        #----------------------------------------------------------------
        def adiabat(self,tpot,*args,**kwargs): #numerical evaluation of adiabat through mantle, tpot potential temperature in K
            z0 = 0.
            zn = self.Rs-self.Rc
            solution = integrate.RK45(self.adiabat_grad, z0, [tpot], zn, rtol=0.001, atol=1e-06)

            depth = []
            temp = []
            while not (solution.status=='finished'):
                # get solution step state
                solution.step()
                depth.append(solution.t)
                temp.append(solution.y[0])

            return self.Rs-np.array(depth),np.array(temp)

        #----------------------------------------------------------------
        def pressure(self,r,*args,**kwargs): #for calling outside, IN: radius in m, out hydrostatic pressure in GPa
            return self.fpress((self.Rs-r)/1000.)

#----------------------------------------------------------------
#do i have all necessities?

import os
def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


necess=['MoonWeber_Science_2011.nd','MoonKhan_JGR_2014_mean.nd','MoonMatsumoto_GRL_2015.nd','MoonKhan_JGR_2014_std.nd','MoonMatsumoto_GRL_2015_std.nd']

if not (all([is_non_zero_file(nc) for nc in necess])):
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('I am missing files, I am sorry I cannot continue')
    print('please download *.nd files from https://doi.org/10.5281/zenodo.3372489 into the current directory')
    print("and try `python3 basal.py` again")
    print('files needed:')
    [print(nc) for nc in necess if not(is_non_zero_file(nc))]
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    exit()


#inicialization of class Moon containing Moon parameters
moon = Moon()

tpot=1600.
r_adiab,t_adiab = moon.adiabat(tpot)
print('----------------------------------------------')
print('adiabat temperature increase estimate')
print('surface temperature',tpot)
print('bottom temperature',t_adiab[-1])
print("estimated adiabatic gradient in Moon's mantle",t_adiab[-1]-tpot)
print('----------------------------------------------')

#inicialization of class for computing conductive temperature profile, numerically, matrix propagator
#automaticly used BC on temperature unless stated in kwargs mode_bot='flux' or mode_top='flux'
mode_temp = temp.Module_conduction_steady(moon.r,moon.k,moon.h,moon.Ts,moon.qbot,mode_bot='flux')
mode_temp.solve() #finding solution

#----------------------------------------------------------------
#radius for plotting
nrad = 101
rarr = np.linspace(moon.Rc,moon.Rs,nrad,endpoint=True)

#inicialization of classes computing solidus and liquidus based on papers
katz = melting.Katz2003()
herzberg = melting.Herzberg2000()
wyatt = melting.Wyatt1977()

#----------------------------------------------------------------
#----------------------------------------------------------------
#inversion results initialization
inv_results = process_results.Process_results(moon)
#----------------------------------------------------------------
#----------------------------------------------------------------


#----------------------------------------------------------------
#----------------------------------------------------------------
print('plotting porosity dependence')
#melting Kervazo style
fig, ax = plt.subplots(1, 2, figsize=figsize2)
fig.subplots_adjust(left=0.1,bottom=0.15,top=0.95,right=0.95,wspace=0.25)
#model Kervazo et al., initialization
kervazo2021 = por.Porosity()

phi_crit = 0.3

npor=100
por_max = 0.48 #max porosity for plotting

porosity = np.linspace(0.,por_max,npor) #porosity numpy array

#-------------------------------------
#viscosity scaled by viscosity of solidus as function of porosity (scaled by critical porosity)
eta_por = kervazo2021.eta_por(porosity)
ax[0].plot(porosity/phi_crit,np.log10(eta_por/eta_por[0]),label = 'Kervazo et al. (2021)',lw=lw,c='indianred')
inv_results.plot_log_rel_eta(ax[0]) #inversion resutls

#making the plot nicer and understendable
ax[0].set_xlim(left = 0.,right=por_max/phi_crit)
ax[0].set_ylim(bottom=-20,top=0.5)
ax[0].set_xlabel(r'relative porosity $\phi\,/\,\phi_c$',fontsize=fs)
ax[0].set_ylabel(r'$\log\;\eta(\phi)\,/\,\eta_{\rm{solid}}$',fontsize=fs)
ax[0].tick_params(axis="x", labelsize=fs)
ax[0].tick_params(axis="y", labelsize=fs)


#-------------------------------------
#shear modulus scaled by shear modulus at solidus as function of porosity (scaled by critical porosity)
mu_por = kervazo2021.mu_por(porosity)
ax[1].plot(porosity/phi_crit,np.log10(mu_por/mu_por[0]),label = 'Kervazo et al. (2021)',lw=lw,c='indianred')
inv_results.plot_log_rel_mu(ax[1]) #inversion results

#making the plot nicer and understendable
ax[1].set_xlim(left = 0.,right=por_max/phi_crit)
ax[1].set_ylim(bottom=None,top=None)
ax[1].set_xlabel(r'relative porosity $\phi\,/\,\phi_c$',fontsize=fs)
ax[1].set_ylabel(r'$\log\;\mu(\phi)\,/\,\mu_{\rm{solid}}$ [Pa]',fontsize=fs)
ax[1].legend(loc='best',fontsize=fs2)
ax[1].tick_params(axis="x", labelsize=fs)
ax[1].tick_params(axis="y", labelsize=fs)

#saving the plot
fig.savefig('Molten-melt-basal.png',dpi=300)
if (show_results): plt.show()
plt.close()


#----------------------------------------------------------------
#----------------------------------------------------------------
#viscosity, dry viscosity depending on T, wet viscosity depending on temperature
print('plotting viscosity depending on temperature')
fig, ax = plt.subplots(1, 2, figsize=(17, 6.))
fig.subplots_adjust(left=0.1,bottom=0.15,top=0.95,right=0.92,wspace=0.15)

#minmax for viscosity plotting
etamax=20.5
etamin=12.1

#aggregates
c_max = 0.16
c_min = 0.02

#experimental results
hirth1996 = viscosity.Hirth1996()
dygert2016 = viscosity.Dygert2016()
visco = viscosity.Viscosity()

#temperature field
ntemp=100 #number of temperatures
temp = np.linspace(1200,2400,ntemp) #in K
#setting other parameters
stress=1. #differential stress in MPa
press = moon.pressure(700e3)*1e6 #moon.pressure is providing pressure in GPa 

#------------------------------------------dry olivine
#viscosity of olivine an dilmenite
ax[0].plot(temp,np.log10(visco.visc_disl(temp,stress,press,hirth1996.par_disl_dry)),label = 'dry olivine, Kohlsted et al. (1996)',lw=lw,c='indianred')
ax[0].plot(temp,np.log10(visco.visc_disl(temp,stress,press,dygert2016.par_disl)),label = 'Dygert et al. (2016)',lw=lw,c='dodgerblue')

#aggregates
smin = np.log10(visco.isostress_disl(temp,stress,press,c_max,dygert2016.par_disl,hirth1996.par_disl_dry))
smax = np.log10(visco.isostress_disl(temp,stress,press,c_min,dygert2016.par_disl,hirth1996.par_disl_dry))
ax[0].fill_between(temp, smin,smax,alpha=0.5,label=r'isostress/Reuss model',color='powderblue',zorder=10)

smin = np.log10(visco.tullis1991(temp,stress,press,c_max,dygert2016.par_disl,hirth1996.par_disl_dry))
smax = np.log10(visco.tullis1991(temp,stress,press,c_min,dygert2016.par_disl,hirth1996.par_disl_dry))
ax[0].fill_between(temp, smin,smax,alpha=0.2,label=r'Tullis model',color='r',zorder=10)

#results
inv_results.plot_log_eta_layer(ax[0])

#solidus and liquidus comparison
ax[0].vlines(np.array([katz.solidus(moon.pressure(700e3),0)+273.15,katz.solidus(moon.pressure(moon.Rc),0)+273.15]),ymin=etamin,ymax=etamax,color='indianred',lw=lw/2,ls='--',label='peridotite, dry solidus, Katz et al. (2003)')
ax[0].vlines(np.array([wyatt.solidus(moon.pressure(700e3))+273.15,wyatt.solidus(moon.pressure(moon.Rc))+273.15]),ymin=etamin,ymax=etamax,color='dodgerblue',lw=lw/2,ls='--',label='ilmenite bearing, solidus, Wyatt (1977)')


#making the plot undersdandable and nicer
ax[0].set_xlim(left = 1200.,right=2400.)
ax[0].set_ylim(bottom=etamin,top=etamax)
ax[0].set_xlabel(r'$T$ [K]', fontsize=fs)
ax[0].set_ylabel(r'$\log\;\eta$ [Pa s]', fontsize=fs)
ax[0].tick_params(axis="x", labelsize=fs)
ax[0].tick_params(axis="y", labelsize=fs)
ax[0].set_title("Dry olivine", pad=20)


#------------------------------------------wet olivine

#olivine and ilmenite
line1, = ax[1].plot(temp,np.log10(visco.visc_disl(temp,stress,press,hirth1996.par_disl_wet)),label = 'olivine, Hirth and Kohlstedt (1996)',lw=lw,c='indianred')
line2, = ax[1].plot(temp,np.log10(visco.visc_disl(temp,stress,press,dygert2016.par_disl)),label = 'ilmenite, Dygert et al. (2016)',lw=lw,c='dodgerblue')

#aggregates
smin = np.log10(visco.isostress_disl(temp,stress,press,c_max,dygert2016.par_disl,hirth1996.par_disl_wet))
smax = np.log10(visco.isostress_disl(temp,stress,press,c_min,dygert2016.par_disl,hirth1996.par_disl_wet))
line3 = ax[1].fill_between(temp, smin,smax,alpha=0.5,label=r'isostress/Reuss model',color='powderblue',zorder=10)
smin = np.log10(visco.tullis1991(temp,stress,press,c_max,dygert2016.par_disl,hirth1996.par_disl_wet))
smax = np.log10(visco.tullis1991(temp,stress,press,c_min,dygert2016.par_disl,hirth1996.par_disl_wet))
line4 = ax[1].fill_between(temp, smin,smax,alpha=0.2,label=r'Tullis model',color='r',zorder=10)

#comparing to solidus and liquidus
line5 = ax[1].vlines(np.array([katz.solidus(moon.pressure(700e3),100)+273.15,katz.solidus(moon.pressure(moon.Rc),100)+273.15]),ymin=etamin,ymax=etamax,color='indianred',ls='--',lw=lw/2,label='peridotite, solidus, Katz et al. (2003)')
line6 = ax[1].vlines(np.array([wyatt.solidus(moon.pressure(700e3))+273.15,wyatt.solidus(moon.pressure(moon.Rc))+273.15]),ymin=etamin,ymax=etamax,color='dodgerblue',ls='--',lw=lw/2,label='ilmenite-bearing, solidus, Wyatt (1977)')

#plotting results
inv_results.plot_log_eta_layer(ax[1])

#making the plot nicer and understandable
ax[1].set_xlim(left = 1200.,right=2400.)
ax[1].set_ylim(bottom=etamin,top=etamax)
ax[1].set_xlabel(r'$T$ [K]')
ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs2)
ax[1].set_title("Wet olivine", pad=20)
ax[1].tick_params(axis="x", labelsize=fs)
ax[1].tick_params(axis="y", labelsize=fs)

fig.savefig('Molten-eta-basal.png',dpi=300, bbox_inches='tight')
if (show_results): plt.show()
plt.close()

#----------------------------------------------------------------
#----------------------------------------------------------------
print('plotting shear modulus')
#plotting three figures derived shear modulus
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.20,bottom=0.15,top=0.95,right=0.95,wspace=0.1)

#data from previous studies
datas=['MoonWeber_Science_2011.nd','MoonKhan_JGR_2014_mean.nd','MoonMatsumoto_GRL_2015.nd']
errors=[None,'MoonKhan_JGR_2014_std.nd','MoonMatsumoto_GRL_2015_std.nd']
colors =['green','indianred','dodgerblue']
describes = ['Weber et al. (2011)','Khan et al. (2014)','Matsumoto et al. (2015)']

eps=1e-8 #numerical epsilon if shear modulus of the core is set to be zero
for (dat_name,err_name,color,describe) in zip(datas,errors,colors,describes):
    
    data=np.genfromtxt(dat_name)
    vs = data[:,2]
    rho = data[:,3]*1000.
    mu=((1e3*vs)**2*rho)

    ax.plot(moon.Rs/1000.-data[:,0],np.log10(mu+eps),c=color,label=describe,lw=lw)

    if err_name:
        data_std=np.genfromtxt(err_name)
        vs_std = data_std[:,2]
        rho_std = data_std[:,3]*1000.
        if ('Khan' in dat_name):
            rho_std /=1000.
        mu_std = 2*(1000.*vs)*(1000.*vs_std)*rho+(1000*vs)**2*rho_std    
        ax.plot(moon.Rs/1000.-data[:,0],np.log10(mu+mu_std+eps),c=color,lw=lw2,ls='--')
        ax.plot(moon.Rs/1000.-data[:,0],np.log10(mu-mu_std+eps),c=color,lw=lw2,ls='--')

#plotting results
inv_results.plot_r_depend(ax,mode="log_mu")

#making the plot understendable and nicer
ax.set_xlabel(r'$r$ [km]',fontsize=fs)
ax.set_xlim(left = 330.,right=1000.)
ax.set_ylim(bottom=9.5,top=11)
ax.set_ylabel(r'$\log\;\mu$ [Pa]',fontsize=fs)
ax.legend(loc='lower right',fontsize=fs2)
ax.tick_params(axis="x", labelsize=fs)
ax.tick_params(axis="y", labelsize=fs)

fig.savefig('Molten-mu-basal.png',dpi=300)
if (show_results): plt.show()
plt.close()

#----------------------------------------------------------------
#----------------------------------------------------------------
#plotting temperature profiles 
print('plotting temperature profiles')
cenrich_max=2 #maximum enrichment in radiogenic elements in lower layer, the global production always remains the same
nums = 200 #number of models

#cm = plt.get_cmap('viridis')
#cm = plt.get_cmap('Blues_r')
cm = plt.get_cmap('pink')

fig = plt.figure(figsize=figsize_complex)
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.15,bottom=0.15,top=0.95,right=0.95,wspace=0.1)

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)

#trying nums models of different enrichement
cenrich=np.linspace(1.,cenrich_max,nums)

#average heating 1
moon.hav = 8e-9
for i,c in enumerate(cenrich):
    moon.Cenrich = c
    moon.update()
    mode_temp.h = moon.h
    mode_temp.solve() #finding solution
    tarr = mode_temp.values_r(mode_temp.tpoint,rarr)
    lines = ax.plot(rarr/1e3,tarr,c=cm(float(i)/nums),zorder=-10,alpha=0.5)

#average heating 2
moon.hav = 9.5e-9
for i,c in enumerate(cenrich):
    moon.Cenrich = c
    moon.update()
    mode_temp.h = moon.h
    mode_temp.solve() #finding solution
    tarr = mode_temp.values_r(mode_temp.tpoint,rarr)
    lines = ax.plot(rarr/1e3,tarr,c=cm(float(i)/nums),zorder=-20,alpha=0.5)

#average heating 3
moon.hav = 11e-9
for i,c in enumerate(cenrich):
    moon.Cenrich = c
    moon.update()
    mode_temp.h = moon.h
    mode_temp.solve() #finding solution
    tarr = mode_temp.values_r(mode_temp.tpoint,rarr)
    lines = ax.plot(rarr/1e3,tarr,c=cm(float(i)/nums),zorder=-30,alpha=0.5)


norm = mpl.colors.Normalize(vmin=1.,vmax=cenrich_max)
sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])

# Plot colorbar
cbar = plt.colorbar(sm, pad=0.01)
cbar.set_label(label=r'basal enrichment factor $f$')

#comparison with melting
press_arr = moon.pressure(rarr)

satur = katz.saturation(press_arr)
ax.plot(rarr/1000.,katz.solidus(press_arr,0.)+273.15,label='solidus, dry peridotite, Katz et al. (2003)',lw=lw,c='dodgerblue')
ax.plot(rarr/1000.,katz.solidus(press_arr,satur)+273.15,label='solidus, water saturated peridotite, Katz et al. (2003)',lw=lw,c='dodgerblue', ls='dotted')
ax.plot(rarr/1000.,katz.liquidus(press_arr,satur)+273.15,label='liquidus, peridotite, Katz et al. (2003)',lw=lw,c='dodgerblue',ls=(0, (5, 1)))
ax.plot(rarr/1000.,wyatt.solidus(press_arr)+273.15,label='solidus, ilmenite bearing, Wyatt (1977)',lw=lw,c='powderblue')
ax.plot(rarr/1000.,wyatt.liquidus(press_arr)+273.15,label='liquidus, ilmenite bearing, Wyatt (1977)',lw=lw,c='powderblue',ls=(0, (5, 1)))

#khan et al. 2014 model
ax.fill_between(moon.Rs/1000.-khan_digi_depth2, khan_digi_min2+273.15, khan_digi_max2+273.15,alpha=0.5,label=r'Khan et al. (2014)',color='gray',zorder=10)

#making the graph understendable and nicer
ax.set_xlabel(r'$r$ [km]')
ax.set_ylabel(r'$T$ [K]')
ax.set_xlim(left=370,right=moon.Rs/1e3)
ax.legend(loc='lower left',fontsize=22)
ax.tick_params(axis="x", labelsize=30)
ax.tick_params(axis="y", labelsize=30)
for t in cbar.ax.get_yticklabels():
      t.set_fontsize(30)


fig.savefig('Molten-temp-basal-diss.png',dpi=300, bbox_inches = 'tight')
    
if show_results: plt.show()
plt.close()



#----------------------------------------------------------------
#----------------------------------------------------------------
'''press = np.linspace(0,8,100) #in GPa
katz = melting.Katz2003()

fig = plt.figure()
ax = fig.add_subplot(111)
satur = katz.saturation(press)
ax.plot(press,katz.solidus(press,0.),label='dry peridotite, Katz et al. (2003)')
ax.plot(press,katz.solidus(press,30.),label=r'0.3\%wt peridotite, Katz et al. (2003)')
ax.plot(press,katz.solidus(press,satur),label='water saturated peridotite, Katz et al. (2003)')
#ax.plot(press,katz.solidus(press)-katz.water_corr(0.3),label='water saturated peridotite, Katz et al. (2003)')
ax.set_xlabel('pressure (GPa)')
ax.set_ylabel(r'solidus temperature ($^\circ$C)')
plt.show()

exit()'''

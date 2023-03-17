#!/usr/bin/env python3
#
#coded by Marie Behounkova marie.behounkova@matfyz.cuni.cz
#modules for expressing solidus and liquidus, different authors

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


#https://link.springer.com/article/10.1007/BF00375941
class Wyatt1977(object):
    #PRESSURE IN GPA
    #TEMPERATURE IN CELSIUS
        def __init__(self, *args, **kwargs):
            #tsolidus, parameterized in deg C
            self.A1 = 1300.            
            self.A2 = 1470

            self.P1 = 2. #GPa
            self.P2 = 4.7 #GPa

            #tsolidus, parameterized in deg C
            self.B2 = 1620.
            self.slope = 60.

        def solidus(self,press,*args,**kwargs):
            return np.where(np.logical_and(press>=self.P1,press<5), self.A2+(press-self.P2)*self.slope,float('NaN'))

        def liquidus(self,press,*args,**kwargs):
            return np.where(np.logical_and(press>=self.P1,press<5), self.B2+(press-self.P2)*self.slope,float('NaN'))


#https://doi.org/10.1029/2002GC000433
class Katz2003(object):
    #PRESSURE IN GPA
    #TEMPERATURE IN CELSIUS
        def __init__(self, *args, **kwargs):
            #tsolidus, parameterized in deg C
            self.A1 = 1085.7
            self.A2 = 132.9
            self.A3 = -5.1

            #correction for water
            self.K = 43.
            self.gamma = 0.75

            #saturation concentration
            self.xi1 = 12.
            self.xi2 = 1.
            self.lamb = 0.6

            #tliquidus (true)
            self.C1 = 1780.
            self.C2 = 45
            self.C3 = -2

        def solidus(self,press,conc,*args,**kwargs):
            dry = self.A1+self.A2*press+self.A3*press**2
            conc_sat = self.saturation(press)
            corr = self.water_corr(conc)
            corr_sat = self.water_corr(conc_sat)
            return np.where(conc<conc_sat,dry-corr,dry-corr_sat)

        def liquidus(self,press,conc,*args,**kwargs):
            return self.C1+self.C2*press+self.C3*press**2

        def saturation(self,press,*args,**kwargs):
            return self.xi1*np.power(press,self.lamb)+self.xi2*press

        def water_corr(self,conc,*args,**kwargs):
            return self.K*np.power(conc,self.gamma)


class Hirschman2000(object):
    #PRESSURE IN GPA
    #TEMPERATURE IN CELSIUS
        def __init__(self, *args, **kwargs):
            #tsolidus, parameterized in deg C
            self.a = -5.1404654
            self.b = 132.899012
            self.c = 1120.66061

        def solidus(self,press,*args,**kwargs):
            return self.a*press**2+self.b*press+self.c

class Herzberg2000(object):
    #PRESSURE IN GPA
    #TEMPERATURE IN CELSIUS
        def __init__(self, *args, **kwargs):
            #tsolidus, parameterized in deg C
            self.a = 1086
            self.b = -5.7
            self.c = 390.

        def solidus(self,press,*args,**kwargs):
            return self.a-self.b*press+390*np.log(press)

class Shear_mavko1980(object):
    def __init__(self, *args, **kwargs):
        self.mu_s = 70e9 #in Pa
        self.poiss = 0.25
        self.factor = (40.-24*self.poiss)/15.

    def shear(self, porosity,*args, **kwargs):
        return self.mu_s*(1+porosity*self.factor)**-1

    def shear_contrast(self, porosity,*args, **kwargs):
        return (1+porosity*self.factor)**(-1)

    def porosity(self, contrast,*args, **kwargs):
        return (contrast-1.)/self.factor


class Viscosity_exp(object):
    def __init__(self, *args, **kwargs):
        self.a = 26 #in dimensionless
        self.eta_s = 1e20 # in Pa s

    def viscosity(self, porosity,*args, **kwargs):
        return self.eta_s*np.exp(-self.a*porosity)

    def viscosity_contrast(self, porosity,*args, **kwargs):
        return np.exp(-self.a*porosity)

    def porosity(self, contrast,*args, **kwargs):
        return np.log(contrast)/(-self.a)


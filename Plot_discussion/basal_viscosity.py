#!/usr/bin/env python3
#
#coded by Marie Behounkova marie.behounkova@matfyz.cuni.cz
#module computing viscosities from different authors
#module computing viscosity of aggregate

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import constants


#--------------------------------------------------
#https://doi.org/10.1016/0012-821X(96)00154-9
class Hirth1996(object): #olivine
    #PRESSURE IN MPa
        def __init__(self, *args, **kwargs):
            #dislocation creep, dry
            self.A_disl_dry = 4.85e4 #MPa^-n
            self.Eact_disl_dry = 535e3 #J/mol
            self.Vact_disl_dry = 0.
            self.n_disl_dry = 3.5
            self.par_disl_dry = [self.A_disl_dry,self.Eact_disl_dry,self.Vact_disl_dry,self.n_disl_dry]


            #dislocation creep, wet
            self.A_disl_wet = 4.89e6 #MPa^-1
            self.Eact_disl_wet = 515e3 #J/mol
            self.Vact_disl_wet = 0.
            self.n_disl_wet = 3.5
            self.par_disl_wet = [self.A_disl_wet,self.Eact_disl_wet,self.Vact_disl_wet,self.n_disl_wet]


#--------------------------------------------------
#https://doi.org/10.1002/2015GL066546
class Dygert2016(object): #ilmenite
    #PRESSURE IN MPa
        def __init__(self, *args, **kwargs):
            #dislocation creep, 1.85 ± 1.68, E = 307 ± 26 kJ/mol, and n = 3.0 ± 0.3,
            self.lnA_disl = 1.85            
            self.lnA_disl_err = 1.68
            self.Eact_disl = 307e3 #J/mol
            self.Eact_disl_err = 26e3 #J/mol
            self.Vact_disl = 0.            
            self.n_disl = 3.
            self.n_disl_err = 0.3

            self.A_disl = np.exp(1.85)

            self.par_disl = [self.A_disl,self.Eact_disl,self.Vact_disl,self.n_disl]


#--------------------------------------------------
class Viscosity(object):
    def __init__(self, *args, **kwargs):
        self.R = constants.R
        self.fact = 1.

    def visc_disl(self,temp,sigma,press,par,*args,**kwargs): #pressure must be in MPa
        parA=par[0];parE=par[1];parV=par[2];parn=par[3]
        return self.fact/parA/np.power(sigma,parn-1)*np.exp((parE+parV*press)/self.R/temp)*1e6

    #harmonic mean
    def isostress_disl(self,temp,sigma,press,c1,par1,par2,*args,**kwargs):
        c2 = 1.-c1
        aux = c1*par1[0]*np.power(sigma,par1[3]-1)*np.exp(-(par1[1]+par1[2]*press)/self.R/temp)
        aux+= c2*par2[0]*np.power(sigma,par2[3]-1)*np.exp(-(par2[1]+par2[2]*press)/self.R/temp)
        return 1./aux*self.fact*1e6

    #https://doi.org/10.1029/90JB02491
    #geometric mean
    def tullis1991(self,temp,sigma,press,c1,par1,par2,*args,**kwargs):
        c2 = 1.-c1
        eta1=self.visc_disl(temp,sigma,press,par1)
        eta2=self.visc_disl(temp,sigma,press,par2)

        return eta1**c1*eta2**c2

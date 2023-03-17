#!/usr/bin/env python3
#
#coded by Marie Behounkova marie.behounkova@matfyz.cuni.cz
#module relating melt represented by the porosity and shear modulus and viscosity

from __future__ import print_function
import numpy as np
import scipy.special as spec #erf


#================================================================================================
#class computing viscosity and shear modulus depending on porosity
#based on Kervazo et al. 2021 parameterization
#https://doi.org/10.1051/0004-6361/202039433 
class Porosity(object):
    def __init__(self, *args, **kwargs):

        self.par_model = "Kervazo2021-default"

        #=======================================

        if (self.par_model == "Kervazo2021-default"):
            self.phi_c = 0.3
            self.eta_phicrit = self.phi_c
            self.mu_phicrit  = self.phi_c
            self.k_phicrit   = self.phi_c

            self.B_einstein= 2.5 #einstein parameter

            #table 3 in Kervazo et al. (2021)
            self.par_eta = [1.,25.7,1.17e-9,5.,1-0.569]
            self.par_mu  = [10.,2.10,7.08e-7,5.,1-0.597]
            self.par_k   = [1e9,2.62,.102,5.,1-0.712]
        else:
            print('unknown mode of the porosity')
            exit()

    #================================================================================================
    #-------------------------------------------------------
    def mavko(self,phi,A,B,*args,**kwargs): #see also Kervazo et al. (2021), eqs A.2 and A.3 and references within
        return A/(1+phi*B)

    #-------------------------------------------------------
    def kervazo_eq1(self,phi,par,*args,**kwargs):
        varl=par[0]
        delta=par[1]
        xi=par[2]
        gamma=par[3]
        phistar=par[4]
    
        theta=(1-phi)/(1.-phistar)
        f=self.kervazo_eq3(theta,xi,gamma,*args, **kwargs)
        return varl*(1+theta**delta)/(1-f)**(self.B_einstein*(1-phistar))

    #-------------------------------------------------------
    def kervazo_eq3(self,theta,xi,gamma,*args,**kwargs):
        arg = np.sqrt(np.pi)/(2.*(1-xi))*theta*(1+theta**gamma)
        return (1-xi)*spec.erf(arg)

    #================================================================================================
    #viscosity as a function of porosity
    def eta_por(self,phi,*args, **kwargs):
        return self.kervazo_eq1(phi,self.par_eta) #viscosity

    #================================================================================================
    #shear modulus as a function of porosity
    def mu_por(self,phi,*args, **kwargs):
        return self.kervazo_eq1(phi,self.par_mu) #shear modulus

    #================================================================================================
    #bulk modulus as a function of porosity
    def K_por(self,phi,*args, **kwargs):
        return self.kervazo_eq1(phi,self.par_k) #bulk modulus
#!/usr/bin/env python3
#
#coded by Marie Behounkova marie.behounkova@matfyz.cuni.cz
#computation of steady state temperature, spherical geometry, 1D conductive profile, piecewise constant properties
#using propagator method


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


#----------------------------------------------------------------
class Module_conduction_steady(object):
        #----------------------------------------------------------------
        #r,radii indexed from the bottom boundary to the top boundary, nlayer + 1 values, indexed from bottom to top
        #k,h conductivity and volumetric heating h for each individual layer, indexed from the bottom boundary to the top boundary, nlayer values
        #otop value of top boundary condition, in K for mode_top=='temp' and in W/m^2 for mode_top=='flux'
        #obot value of bottom boundary condition, in K for mode_bot=='temp' and in W/m^2 for mode_bot=='flux'
        #mode
        def __init__(self, r, k, h, otop,obot, mode_top='temp',mode_bot='temp',*args, **kwargs):
            self.nlayer = k.shape[0]
            self.nh     = h.shape[0]
            self.n      = r.shape[0]
            if ((self.n-self.nlayer)!=1) or (self.nh!=self.nlayer):
                print('non-acceptable input, check dimensions of k, r, h')
                exit()

            self.r = r
            self.h = h
            self.k = k

            self.otop = otop
            self.obot = obot

            self.vec_sol = np.zeros((self.nlayer,2))


            if mode_top=='temp':
                self.Pt = self.prop_T()
            elif mode_top=='flux':
                self.Pt = self.prop_q()
            else:
                print('unknown boundary condition top')
                exit()

            if mode_bot=='temp':
                self.Pb = self.prop_T()
            elif mode_bot=='flux':
                self.Pb = self.prop_q()
            else:
                print('unknown boundary condition bottom')
                exit()

        #----------------------------------------------------------------
        def solve(self,*args,**kwargs):
            self.find_coefs_one()
            for i in range(1,self.nlayer):
                self.vec_sol [i+1-1,:] = np.dot(self.matM_inv(i+1,self.r[i]),np.dot(self.matM(i,self.r[i]),self.vec_sol[i+0-1])+self.vecV(i,self.r[i])-self.vecV(i+1,self.r[i]))

        #----------------------------------------------------------------
        def find_coefs_one(self,*args,**kwargs):
            mataux = np.zeros((2,2))
            vecaux = np.zeros((2))

            Dtilde = np.dot(self.matM(self.nlayer,self.r[-1]),self.matD())
            Vtilde = np.dot(self.matM(self.nlayer,self.r[-1]),self.vecW())+self.vecV(self.nlayer,self.r[-1])

            mataux[0,:] = np.dot(self.Pt,Dtilde)
            vecaux[0]   = self.otop-np.dot(self.Pt,Vtilde)

            mataux[1,:] = np.dot(self.Pb,self.matM(1,self.r[0]))
            vecaux[1]   = self.obot-np.dot(self.Pb,self.vecV(1,self.r[0]))

            self.vec_sol [1-1,:] = np.linalg.solve(mataux,vecaux)

        #----------------------------------------------------------------
        def tpoint(self,raux,*args,**kwargs):
            i=self.layer_number(raux,*args,**kwargs)
            return (np.dot(self.matM(i,raux),self.vec_sol[i-1])+self.vecV(i,raux))[0]
        #----------------------------------------------------------------
        def qpoint(self,raux,*args,**kwargs):
            i=self.layer_number(raux,*args,**kwargs)
            return (np.dot(self.matM(i,raux),self.vec_sol[i-1]) +self.vecV(i,raux))[1]

        #----------------------------------------------------------------
        def values_r(self,func,rarr,*args,**kwargs):
            res = np.zeros_like(rarr)
            for i,raux in enumerate(rarr):
                res[i]=func(raux)
            return res

        #----------------------------------------------------------------
        def layer_number(self,raux,*args,**kwargs):
            return max(1,np.searchsorted(self.r, raux))

        #----------------------------------------------------------------
        def matM(self,i,raux,*args,**kwargs):
            kaux=self.k[i-1]
            return np.array([[-1./kaux/raux,1],[-1./raux**2,0]])

        #----------------------------------------------------------------
        def matM_inv(self,i,raux,*args,**kwargs):
            kaux = self.k[i-1]
            return np.array([[0,-raux**2],[1,-raux/kaux]])

        #----------------------------------------------------------------
        def vecV(self,i,raux,*args,**kwargs):
            haux = self.h[i-1]
            kaux = self.k[i-1]
            return np.array([-haux/kaux*raux**2/6.,haux*raux/3.])

        #----------------------------------------------------------------
        def matD(self,*args,**kwargs):
            mataux = np.identity(2)
            for i in range(1,self.nlayer):
                mataux = np.dot(np.dot(self.matM_inv(i,self.r[i]),self.matM(i+1,self.r[i])), mataux)
            return mataux

        #----------------------------------------------------------------
        def vecW(self,*args,**kwargs):
            vecaux = np.zeros(2)
            for i in range(1,self.nlayer):
                vecaux = np.dot(self.matM_inv(i+1,self.r[i]),self.vecV(i,self.r[i])-self.vecV(i+1,self.r[i])+np.dot(self.matM(i,self.r[i]),vecaux))
            return vecaux

        #----------------------------------------------------------------
        def prop_T(self,*args,**kwargs):
            return np.array([1,0])

        #----------------------------------------------------------------
        def prop_q(self,*args,**kwargs):
            return np.array([0,1])
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#main file, plotting figures, discussion on Moon basal layer
#coded Marie Behounkova marie.behounkova@matfyz.cuni.cz

import numpy as np
# import matplotlib as mpl
# mpl.rc('font',family='Times New Roman',size=24)
# mpl.rcParams['mathtext.fontset'] = 'custom'
# mpl.rcParams['mathtext.rm'] = 'Times New Roman'
# mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
# mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

# import matplotlib.pyplot as plt
# from matplotlib import ticker
# import melting
# import temp
# import viscosity
# import porosity as por
# from scipy import integrate



#-------------------------------------------------------
class Process_results(object):
        #object defining parameters of the Moon
        #----------------------------------------------------------------
        def __init__(self,body,*args, samples = '../Random_100/Molten-random100.csv',best='../Best_fits/Molten-best_fits.csv',**kwargs):
            self.data_best = np.genfromtxt(best,skip_header = 1,delimiter=',')
            self.data_samples = np.genfromtxt(samples,skip_header = 1,delimiter=',')

            self.order = {
                "chi2":0,
                "mu":1,
                "eta":2,
                "alpha":3,
                "zeta":4,
                "rho_mantle":5,
                "rho_core":6,
                "R_core":7,
                "rho_layer":8,
                "R_layer":9,
                "eta_layer":10,
                "mu_layer":11}
            self.body = body
            self.Rs = body.Rs

            self.c_sample = 'lightgrey' #'lightblue' #'lightblue'
            self.c_best = 'k' #'darkblue'

            self.zorder_sample = -10
            self.zorder_best = -5

            self.alpha_sample = 1.
            self.alpha_best = 1.

            self.alpha_sample = 1.
            self.alpha_best = 1.

            self.lw_sample = 1.
            self.lw_best = 1.


        #----------------------------------------------
        def plot_r_depend(self,ax,*args,mode='default',**kwargs):
            if (mode=='log_mu'):
                aux_func = self.mu
            elif (mode=='log_eta'):
                aux_func = self.log_eta
            elif (mode=='rho'):
                aux_func = self.rho

            n=10000
            r = np.linspace(0.,self.Rs*1e-3,n)
            aux = np.zeros_like(r)
            for sample in self.data_best:
                for i,auxr in enumerate(r):
                    aux[i]=aux_func(r[i],sample)
                if (mode=='log_mu'):
                    aux=np.log10(aux)
                ax.plot(r,aux,c=self.c_best,linewidth=self.lw_best,alpha=self.alpha_best,zorder = self.zorder_best)

            for sample in self.data_samples:
                for i,auxr in enumerate(r):
                    aux[i]=aux_func(r[i],sample)
                if (mode=='log_mu'):
                    aux=np.log10(aux)
                ax.plot(r,aux,c=self.c_sample,linewidth=self.lw_sample,alpha=self.alpha_sample,zorder = self.zorder_sample)
               

            

        #----------------------------------------------
        def mu(self,r,sample,*args,**kwargs):
            if r<sample[self.order["R_core"]]:
                return None
            elif r<sample[self.order["R_layer"]]:
                return sample[self.order["mu_layer"]]*1e9
            elif r<self.Rs:
                return sample[self.order["mu"]]*1e9
            else:
                return None

        #----------------------------------------------
        def log_eta(self,r,sample,*args,**kwargs):
            if r<sample[self.order["R_core"]]:
                return None
            elif r<sample[self.order["R_layer"]]:
                return sample[self.order["eta_layer"]]
            elif r<self.Rs:
                return sample[self.order["eta"]]
            else:
                return None
            
        #----------------------------------------------
        def rho(self,r,sample,*args,**kwargs):
            if r<sample[self.order["R_core"]]:
                return sample[self.order["rho_core"]]
            elif r<sample[self.order["R_layer"]]:
                return sample[self.order["rho_layer"]]
            elif r<self.Rs-self.body.d_crust:
                return sample[self.order["rho_mantle"]]
            elif r<self.Rs:
                return self.body.rho_crust
            else:
                return None

        #----------------------------------------------
        def plot_log_rel_mu(self,ax,*args,**kwargs):
            for sample in self.data_best:
               ax.axhline(np.log10(sample[self.order["mu_layer"]]/sample[self.order["mu"]]),c=self.c_best,lw=self.lw_best,zorder=self.zorder_best,alpha=self.alpha_best)
            for sample in self.data_samples:
               ax.axhline(np.log10(sample[self.order["mu_layer"]]/sample[self.order["mu"]]),c=self.c_sample,lw=self.lw_sample,zorder=self.zorder_sample,alpha=self.alpha_sample)

        #----------------------------------------------
        def plot_log_rel_eta(self,ax,*args,**kwargs):
            for sample in self.data_best:
               ax.axhline(sample[self.order["eta_layer"]]-sample[self.order["eta"]],c=self.c_best,lw=self.lw_best,zorder=self.zorder_best,alpha=self.alpha_best)
            for sample in self.data_samples:
               ax.axhline(sample[self.order["eta_layer"]]-sample[self.order["eta"]],c=self.c_sample,lw=self.lw_sample,zorder=self.zorder_sample,alpha=self.alpha_sample)

        #----------------------------------------------
        def plot_log_eta_layer(self,ax,*args,**kwargs):
            for sample in self.data_best:
               ax.axhline(sample[self.order["eta_layer"]],c=self.c_best,lw=self.lw_best,zorder=self.zorder_best,alpha=self.alpha_best)
            for sample in self.data_samples:
               ax.axhline(sample[self.order["eta_layer"]],c=self.c_sample,lw=self.lw_sample,zorder=self.zorder_sample,alpha=self.alpha_sample)
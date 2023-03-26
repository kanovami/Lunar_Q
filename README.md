# Lunar_Q

This repository contains the software used for the study "Is there a semi-molten layer at the base of the lunar mantle?" by Michaela Walterová, Marie Běhounková, and Michael Efroimsky. Below are the specifications of each directory, additional references, and brief tutorials.

## Best_fits

Overview of the ten best-fitting models for Models 1, 2, and 3, and for a reparameterised version of Model 2, where the grain size and grain-boundary viscosity are used instead of the relaxation time $\tau$.

## MCMC

Files used for running the MCMC inversion with Models 1, 2, and 3 as well as with the reparameterised Model 2. The output file is always called "chain_core_REAL.dat" and lists every *100th* accepted sample, starting with sample number 1001. Before running any of the files "k2_mcmc*.py", please use the existing Fortran files to create a Python module "tidal_deformation". This can be done with *Numpy* (https://numpy.org/) and its utility *f2py* (https://numpy.org/doc/stable/f2py/). A previously installed Fortran compiler is necessary.

In this example, we use the Fortran compiler *GFortran* (https://gcc.gnu.org/wiki/GFortran). A simple way to create the desired module is to first run the following command:

`gfortran -c nr.for`,

followed with

`python -m numpy.f2py -c nr.for mrheo.f90 -m tidal_deformation`.

Please ignore any warning arising from "nr.for" - this will be corrected soon.

The Python files also require previously installed *emcee* (https://emcee.readthedocs.io/en/stable/) and *corner* (https://corner.readthedocs.io/en/latest/) modules.

## Plot_auxiliary

Files used for plotting Figures 2, 3, and 4 from Subsection 5.2. As in the previous case, the Python files require the module "tidal_deformation", which can be created from the Fortran files.

## Plot_discussion

Files used for plotting Figures 13-19 from Section 6 (Discussion).

### Discussion of the Sundberg-Cooper model - Subsection 6.1

Figure 13 can be reproduced with "relaxation_time.py". Figure 14 was plotted with "mantleQ.py" and requires the module "tidal_deformation" as well as the output from the MCMC run with Model 2.

### Discussion of the model with a basal layer - Subsection 6.2

### Discussion of other sources of information - Subsection 6.3

Figure 19 was plotted with the file "heating.py", which requires the output from the tidal heating code. The code can be found in the folder "Tidal_heating" and comes in two versions: "Surface_heat_flux" and "Volumetric_heating". Both versions are written in *Fortran 90* and necessitate a Fortran compiler.

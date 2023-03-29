# Lunar_Q

This repository contains the software used for the study "Is there a semi-molten layer at the base of the lunar mantle?" by Michaela Walterová, Marie Běhounková, and Michael Efroimsky. Below are the specifications of each directory, additional references, and brief tutorials.

## Best_fits

Overview of the ten best-fitting samples for Models 1, 2, and 3, and for a reparameterised version of Model 2, where the grain size and grain-boundary viscosity are used instead of the relaxation time $\tau$.

## MCMC

Files used for running the MCMC inversion with Models 1, 2, and 3 as well as with the reparameterised Model 2. The output file is always called "chain_core_REAL.dat" and lists every *100th* accepted sample, starting with sample number 1001. Before running any of the files "k2_mcmc*.py", please use the existing Fortran files to create a Python module "tidal_deformation". This can be done with [*Numpy*](https://numpy.org/) and its utility [*f2py*](https://numpy.org/doc/stable/f2py/). A previously installed Fortran compiler is necessary.

In this example, we use the Fortran compiler [*GFortran*](https://gcc.gnu.org/wiki/GFortran). A simple way to create the desired module is to first run the following command:

`gfortran -c nr.for`,

followed with

`python -m numpy.f2py -c nr.for mrheo.f90 -m tidal_deformation`.

Please ignore any warning arising from "nr.for" - this will be corrected soon.

The Python files also require previously installed [*emcee*](https://emcee.readthedocs.io/en/stable/) and [*corner*](https://corner.readthedocs.io/en/latest/) modules.

## Plot_auxiliary

Files used for plotting Figures 2, 3, and 4 from Subsection 5.2. As in the previous case, the Python files require the module "tidal_deformation", which can be created from the Fortran files.

## Plot_discussion

Files used for plotting Figures 13-19 from Section 6 (Discussion).

### Discussion of the Sundberg-Cooper model - Subsection 6.1

Figure 13 can be reproduced with "relaxation_time.py". Figure 14 was plotted with "mantleQ.py" and requires the module "tidal_deformation" as well as the output from the MCMC run with Model 2.

### Discussion of the model with a basal layer - Subsection 6.2

### Discussion of other sources of information - Subsection 6.3

Figure 19 was plotted with the file "heating.py", which requires the output from the tidal heating code. The code can be found in the folder "Tidal_heating" and comes in two versions: "Surface_heat_flux" and "Volumetric_heating". Both versions are written in *Fortran 90* and necessitate a Fortran compiler.

The tidal heating code accepts input in the form indicated by the files "rheo_sc.in" and "rheo_melt.in". Each line of the input file specifies one interior layer, starting with the innermost one (which is, in our case, the core). The parameters of the layer are written in the following order:

1. rheological model (0 - Maxwell, 1 - Andrade, 2 - Sundberg-Cooper)
2. density (in kg/m^3)
3. outer radius (in m)
4. viscosity (in Pa s)
5. rigidity (in Pa)
6. $\alpha$ (for Andrade or Sundberg-Cooper rheology)
7. $\zeta$ (for Andrade or Sundberg-Cooper rheology)
8. relative viscosity of the dashpot in the Kelvin-Voigt element (for the Sundberg-Cooper rheology)
9. relaxation strength $\Delta$ of the secondary peak (for the Sundberg-Cooper rheology)

The last two parameters correspond to $\eta_{\rm{P}}/\eta_{\rm{S}}$ and $\delta J/J_{\rm{U}}$ from [Renaud and Henning (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aab784). The Sundberg-Cooper parameters used in our paper are defined as $\Delta=\delta J/J_{\rm{U}}$ and $t_{\rm{rel}}=\eta_{\rm{P}}\Delta/\eta_{\rm{S}}$.

The input file as well as other free parameters of the tidal heating code are specified in the module "mconst.f90". To compile the code, use all files ending with ".f90" and ".for" - the latter ("nr.for", in particular) contains subroutines from [Numerical Recipes in Fortran 77](http://numerical.recipes) by Press et al. (1992). File "main.f90" is the main body of the code and includes a subroutine for the calculation of tidal heating in a unit volume averaged over one orbital period. In "peltier.f90" is a subroutine calculating complex potential tidal Love number $k_2$ at a given frequency for the prescribed interior structure. Finally, file "glpq.f90" calculates Hansen coefficients (or, specifically, Kaula's eccentricity functions; see [Kaula, 1961](https://academic.oup.com/gji/article/5/2/104/669948)) using the Von Zeipel-Andoyer method (e.g., [Cherniack, 1972](https://ui.adsabs.harvard.edu/abs/1972SAOSR.346.....C/abstract)).

## Plot_results

Files used for plotting Figures 5-12 from Section 5 (Application to the Moon) and similar figures in the Supporting Information, given the output "chain_core_REAL.dat" of the MCMC runs. The same files also generate the lists of ten best-fitting samples and 100 random samples for Models 1-3 and the reparameterised version of Model 2. The module "tidal_deformation", mentioned in one of the previous sections, is required.

## Random_100

Overview of the 100 randomly-chosen samples for Models 1, 2, and 3, and the reparameterised version of Model 2 that are depicted in Figures 7, 9, 15-17, and S1.

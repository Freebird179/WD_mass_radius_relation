# WD_mass_radius_relation
Module for calculating the radius and mass of a white dwarf from effective temperature (Teff) and surface gravity (logg) by interpolating the Argentenian (http://evolgroup.fcaglp.unlp.edu.ar/TRACKS/newtables.html) (Althaus et al. 2013) or Montreal (https://www.astro.umontreal.ca/~bergeron/CoolingModels/) (Bedard et al. 2020) WD evolutionary models to use as a mass-radius relation.  

1. Required packages:Python 3.6 or higher, numpy, scipy, astropy, pandas
 
2. Download the folder WD_mass_radius_relation. Go to the downloaded path in the Terminal and enter python or ipython.

3. Type the following in the terminal :
import WD_MR_relation as MR

4. To calculate the radius and mass from input Teff (in K) and log(g), use the functions MR.Radius_from_teff_logg and MR.Mass_from_teff_logg respectively.

5. Read the documentation for more details.

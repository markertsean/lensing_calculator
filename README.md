This code is to analyze FITS images output by MAGNMS (Mass mAp GeNerator for the Multidark Suite)

5 input files are required (in order of reading):
  
  lensUserParams.dat - Contains variables used by the code, many have default values
  
  foxH2012.dat       - FoxH function table to interpolate over for Einasto fitting
  
  foxH2123.dat       - " "
  
  haloList.dat       - File containing all the FITS images to analyze, DO NOT MIX DIFFERENT HALOS IN THIS FILE
  
  param_build        - Glamer file containing variables used by the program


The flow of the code is as follows:

1 - Read user input file, and set variables that will not be changing over course of program

2 - Read the FoxH tables, as they will also be unchanging over the course of the program

3 - Loop over haloList.dat, combining data from param_build and haloList.dat to build paramfile, which is used by GLAMER

4 - GLAMER sets internal parameters

5 - Read the input FITS file header to retrieve information about the halo and image

6 - Set all remaining parameters to govern the program

7 - Generate PixelMaps used by GLAMER

8 - Construct the Lens using GLAMER

9 - Generate grid onto lens

10 - Extract lensing quantities onto our PixelMaps, from the grid

11 - OUTPUT: 2D maps of reduced shear quantites, binned radially and from orientation of halo

12 - Generate random source locations, remove any within R_vir ( roughly strong lensing regime )

13 - OUTPUT: Source distances and shear values averaged over nearby pixels

14 - Loop over subsets, calculating fits, and omitting subsets for jacknife errors

15 - Generate and populate RTS and dist bins

16 - Fit using "rolling ball" chi-squared fitting of Einasto, NFW, and truncated NFW profile

17 - Calculate the jacknife errors based on subset results

18 - OUTPUT: Density fits and associated errors

The outputs for the 2D maps and sources will be unique to each FITS in haloList. The Density file will contain
information from all the images in haloList, which is why the FITS images should all be the same halo.

#ifndef LENS_FITTER
#define LENS_FITTER

#include <lensing_classes.h>

void fitDensProfile( 
		    lensProfile  &densProfile,  //Density profile we will write results to 
		    haloInfo            &halo,  //Information about the halo we are trying to fit
		    userInfo                u,  //Input information from the user
		    double              *gArr,  //Radially averaged rts measurements, to fit
		    double              *dArr,  //Radial distances from above bins
		    double           *gErrArr,  //Errors in the radiallt averaged rts bins
		    double          *sourceSC,  //Source critical surface densities
		    double        *sourceDist); //Projected radial distances of sources to lens centers

void generateNFWRTS( 
		    double           *gArr,  //Radially averaged RTS array function will return
		    lensProfile      &lens,  //Input density profile
		    haloInfo         &halo,  //Information on the halo
                    userInfo             u,  //Info from user
		    double       *sourceSc,  //Critical surface density of a source
		    double     *sourceDist); //Projected radial distance of sources to lens centers

double    SDNFW(
		double     r,  //Distance to evaluate SD of NFW profile at
		double   r_s,  //Scale radius of profile
		double rho_o); //Initial density of profile

double SDAvgNFW(
		double     r,  //Distance to evaluate SD of NFW profile at
		double   r_s,  //Scale radius of profile
		double rho_o); //Initial density of profile


#endif

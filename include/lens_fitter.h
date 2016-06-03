#ifndef LENS_FITTER
#define LENS_FITTER

#include <lensing_classes.h>


void foo();
void fitDensProfile(
                    densProfile   &  profile ,  // Density profile we are outputting
                    haloInfo      &     halo ,  // Info about parent halo
                    userInfo               u ,  // Info from the user
                    double       *      gArr ,  // RTS binned array we "observed"
                    double       *      dArr ,  // Distance binned array
                    double       *   gErrArr ,  // Error array in RTS
                    double       *sourceSigC ,  // Crit surface densities of sources
                    double       *sourceDist ); // 2D distance from source to lens
/*
void generateNFWRTS(
		    double           *gArr,  //Radially averaged RTS array function will return
		    lensProfile      &lens,  //Input density profile
		    haloInfo         &halo,  //Information on the halo
                    userInfo             u,  //Info from user
		    double       *sourceSc,  //Critical surface density of a source
		    double     *sourceDist); //Projected radial distance of sources to lens centers

void generateEinRTS(
		    double           *gArr,  //Radially averaged RTS array function will return
		    lensProfile      &lens,  //Input density profile
		    haloInfo         &halo,  //Information on the halo
                    userInfo             u,  //Info from user
		    double       *sourceSc,  //Critical surface density of a source
		    double     *sourceDist); //Projected radial distance of sources to lens centers

double    SDNFW(
                double               r ,  //Input radius to calc SD at
                lensProfile inpProfile ); //Input NFW profile

double    SDAvgNFW(
                double               r ,  //Input radius to calc SD at
                lensProfile inpProfile ); //Input NFW profile


double    SDNFWFull(
		double     r,  //Distance to evaluate SD of NFW profile at
		double   r_s,  //Scale radius of profile
		double rho_o); //Initial density of profile

double SDAvgNFWFull(
		double     r,  //Distance to evaluate SD of NFW profile at
		double   r_s,  //Scale radius of profile
		double rho_o); //Initial density of profile

double foxH2012(
                double         z       ,  // Z from fox H function
                double     alpha       ,  // Shape parameter Ein profile
                double tolerance = 1e-2); // Tolence level for convergence

double foxH2123(
                double         z       ,  // Z from fox H function
                double     alpha       ,  // Shape parameter Ein profile
                double tolerance = 1e-2); // Tolence level for convergence
*/

#endif

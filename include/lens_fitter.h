#ifndef LENS_FITTER
#define LENS_FITTER

#include <lensing_classes.h>


void fitDensProfile(
                          densProfile   &  profile ,  // Density profile we are outputting
                    const haloInfo      &     halo ,  // Info about parent halo
//                    const userInfo               u ,  // Info from the user
                    const userInfo               u ,  // Info from the user
                    const double       *      gArr ,  // RTS binned array we "observed"
                    const double       *      dArr ,  // Distance binned array
                    const double       *   gErrArr ,  // Error array in RTS
                    const COSMOLOGY          cosmo );

void generateNFWRTS(
                          double          *gArr ,  // RTS array to output
                    const densProfile     &lens ,  // Input density profile to generate profile for
                    const double         N_bins ,  // Actual information from the halo
                    const double          *dist ,  // Projected distances between source and lens
                    const double           SigC ); // Critical surface density of sources


double    SDNFW( const      double               r ,  // Input radius to calc SD at
                 const densProfile inpProfile ); // Input NFW profile

double    SDAvgNFW( const double               r ,  //Input radius to calc SD at
                    const densProfile inpProfile ); //Input NFW profile


double    SDNFWFull( const double     r ,  //Distance to evaluate SD of NFW profile at
                     const double   r_s ,  //Scale radius of profile
                     const double rho_o ); //Initial density of profile

double SDAvgNFWFull( const double     r ,  //Distance to evaluate SD of NFW profile at
                     const double   r_s ,  //Scale radius of profile
                     const double rho_o ); //Initial density of profile

/*

void generateEinRTS(
		    double           *gArr,  //Radially averaged RTS array function will return
		    lensProfile      &lens,  //Input density profile
		    haloInfo         &halo,  //Information on the halo
                    userInfo             u,  //Info from user
		    double       *sourceSc,  //Critical surface density of a source
		    double     *sourceDist); //Projected radial distance of sources to lens centers


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

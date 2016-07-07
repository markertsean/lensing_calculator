#ifndef LENS_FITTER
#define LENS_FITTER

#include <lensing_classes.h>
#include <gmpxx.h>

#define ln_ln10   2.302585092994045684017991454684364207601101488628
#define ln_ln10s "2.302585092994045684017991454684364207601101488628"
#define ln_e      2.718281828459045235360287471352662497757247093699
#define ln_es    "2.718281828459045235360287471352662497757247093699"


void rollBall(        densProfile   &ball ,  // Ball to roll
                      double        &chi2 ,  // Chi2 value
               const  double        *gArr ,  // g values we are fitting
               const  double        *dArr ,  // distance array corresponding to g values
               const  double     *gErrArr ,  // Errors associated with g
               const  double       sigmaC ,  // Critical surface density to use, of sources
               const  userInfo          u );


void rollingFitDensProfile(
                          densProfile   &  profile ,  // Density profile we are outputting
                    const haloInfo      &     halo ,  // Info about parent halo
                    userInfo               u ,  // Info from the user
//                    const userInfo               u ,  // Info from the user
                    const double       *      gArr ,  // RTS binned array we "observed"
                    const double       *      dArr ,  // Distance binned array
                    const double       *   gErrArr ,
                    const COSMOLOGY          cosmo ); // Error array in RTS


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



void generateEinRTS(
		    double           *gArr,  //Radially averaged RTS array function will return
		    densProfile      &lens,  //Input density profile
        userInfo             u,  //Info from user
		    double        sourceSc,  //Critical surface density of a source
		    double     *sourceDist); //Projected radial distance of sources to lens centers


double foxH2012(
                double         z       ,  // Z from fox H function
                double     alpha       ,  // Shape parameter Ein profile
                double tolerance = 1e-2); // Tolence level for convergence

double foxH2123(
                double         z       ,  // Z from fox H function
                double     alpha       ,  // Shape parameter Ein profile
                double tolerance = 1e-2); // Tolence level for convergence

mpf_class diGamma( mpf_class z );

// Roughly 1e2 - 1e4 times slower than exp
mpf_class exp( mpf_class inpVal );


// Roughly 1e5 times slower than log
mpf_class ln( mpf_class inpVal );


#endif

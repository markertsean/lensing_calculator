/*

File name from commandline, set paramfile variable to filename, allows looping over file for file


Find & loop over Box files


Number of sources
Source error distribution


Seperate section for triaxiality and g+
Uniform distribution of sources, need alignment relative to halo orientation

*/


#include <ctime>
#include <slsimlib.h>
#include <iostream>
#include <cstring>
#include <simpleTree.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
#include <cosmo.h>
#include "grid_maintenance.h"
#include "gridmap.h"

#include <CCfits/CCfits>

//My files
#include "astro_constants.h"
#include "lensing_classes.h"
#include "my_utilities.h"
#include "input_functions.h"
#include "lens_fitter.h"
#include "pixelmap_functions.h"
#include "precise_math.h"

#include <gmp.h>
#include <gmpxx.h>

// Name of the code log file
std::string logFileName = "";


int main(int arg,char **argv){

//int prec = 100;
int prec = 500;
mpf_set_default_prec( prec );

/*
int foo = 50;
std::cout.precision(foo);
std::scientific;
mpf_class piV(ln_pis);

mpf_class bar( -151.1 );

std::cout << std::setw( foo ) << bar             << std::endl;
std::cout << std::setw( foo ) << bar.get_d()     << std::endl;
std::cout << std::setw( foo ) << spouges( bar ) << std::endl;
exit(0);
//*/

  // Initializes the log file, generates logfiles directory
  //  and a file name based on current time
  initLogFile();

  // Default values that will not be changing
  long    seed    = time(NULL); //-1827674;

  logMessage( std::string("Seed = ") + std::to_string( (long long) seed ) );

  srand(seed); // Sets random seed

  //////////////////////////////////
  ////////////READ IN///////////////
  //////////////////////////////////


  haloInfo    myHalo;  // Contains all the info on our halo
  userInfo userInput;  // Contains all the info from the user


  std::string   userFile = "lensUserParams.dat" ; // User specified input
  std::string  paramFile = "paramfile"          ; // Glamer input


  // File for glamer read in, contains information used by the library
  // Must contain variables:
  //    deflection_off          0
  //    field_off               1
  //    lensing_off             0
  //    main_halo_on            0
  //    pixelmaps_input_file    ???.FITS
  //    pixelmaps_on            1
  //    z_source                #

  std::cout << "Reading file: "              << std::endl;
  std::cout << "              " << paramFile << std::endl;

  InputParams      params     ( paramFile   );
  logMessage( std::string("Read paramfile: ") + paramFile );

  std::cout << "              Done."  << std::endl << std::endl;



  // Determine our FITS file name, will read the header for parameters
  std::string                                    fitsFileName  ;
  params.get(          "pixelmaps_input_file"  , fitsFileName );
  logMessage( std::string("Using FITS file: ") + fitsFileName );





  std::cout << "Reading file: "              << std::endl;
  std::cout << "              " << userFile  << std::endl;

  // Reads info from user file
  readInpFile( userInput, userFile );

  std::cout << "              Done."  << std::endl << std::endl;



  std::cout << "Reading file: "                  << std::endl;
  std::cout << "              " << fitsFileName  << std::endl;

  // Read FITS image, populates halo and user inputs
  readFitsHeader( fitsFileName, myHalo, userInput );

  std::cout << "              Done."  << std::endl << std::endl;


  ///////////////////////////////////////////////////
  /////////////INITIALIZE NEEDED PARAMETERS//////////
  ///////////////////////////////////////////////////


  COSMOLOGY cosmo;

  // Read in from lensUserParams, sets the cosmology for glamer
  if ( userInput.getCosmology() == "PLANCK" ) {

    cosmo = Planck1yr;
    logMessage( std::string( "Using Planck cosmology") );
                std::cout << "Using Planck cosmology"   << std::endl << std::endl;

  } else
  if ( userInput.getCosmology() == "WMAP"   ) {

    cosmo = WMAP5yr;
    logMessage( std::string( "Using WMAP 5yr cosmology") );
                std::cout << "Using WMAP 5yr cosmology" << std::endl << std::endl;

  } else {
                std::cout << "Unrecognized cosmology: " << userInput.getCosmology() << std::endl;
    logMessage( std::string( "Unrecognized cosmology: ") + userInput.getCosmology() );
    logMessage( "Aborting." );
    exit(1);
  }

  // Source redshift moves from glamer parameters to userInfo class
  {
    double srcZ;
    params.get("z_source",srcZ);
    userInput.setSourceZ(srcZ);
  }

  double   angRange =   userInput.getAngFOV  ()  * M_PI/180  ;  // Angular field of view, in degrees
  double  realWidth =   userInput.getPhysFOV ()              ;  // Field of view in Mpc
  omp_set_num_threads ( userInput.getNthreads()             );  // For parallelization, default 1
  double center[]   = { 0, 0 };                                 // Center of grid, 0,0 aligns grid right

  userInput.setNgridPoints( userInput.getNpixH() * 2       );   // Number of gridpoints, in reality want it more refined


  logMessage( std::string("angRange        = " ) + std::to_string((long double)  angRange    ));
  logMessage( std::string("realWidth       = " ) + std::to_string((long double)  realWidth   ));
  logMessage( std::string("omp_num_threads = " ) + std::to_string((long double)  userInput.getNthreads()));
  logMessage( std::string("counter         = " ) + std::to_string((long double) center[0]    )
                                                 + std::to_string((long double) center[1]    ));



userInput.setNpixH( 9 );
userInput.setNpixV( 9 );
userInput.setNpix( 9*9);
  std::cout << "Constructing PixelMaps..." << std::endl;

  // PixelMaps we need to keep
  PixelMap  g_tanMap( center, userInput.getNpixH(), userInput.getNpixV(), angRange / userInput.getNpixH() );
  PixelMap  g_secMap( center, userInput.getNpixH(), userInput.getNpixV(), angRange / userInput.getNpixH() );
  PixelMap   distMap( center, userInput.getNpixH(), userInput.getNpixV(), angRange / userInput.getNpixH() );

  {

  // PixelMaps we don't need to keep
  PixelMap  kappaMap( center, userInput.getNpixH(), userInput.getNpixV(), angRange / userInput.getNpixH() );
  PixelMap gamma1Map( center, userInput.getNpixH(), userInput.getNpixV(), angRange / userInput.getNpixH() );
  PixelMap gamma2Map( center, userInput.getNpixH(), userInput.getNpixV(), angRange / userInput.getNpixH() );
  PixelMap invMagMap( center, userInput.getNpixH(), userInput.getNpixV(), angRange / userInput.getNpixH() );

  std::cout << "               Done." << std::endl;

  logMessage( std::string("PixelMaps allocated") );



  ///////////////////////////////////////////////////////////
  ///////////////////////CONSTRUCT LENS//////////////////////
  ///////////////////////////////////////////////////////////


  std::cout <<                           std::endl << std::endl;
  std::cout << "Constructing lens..." << std::endl << std::endl;

  // Generates lens from paramfile, glamer side
  Lens myLens( params, &seed, cosmo );

  std::cout <<                                        std::endl;
  std::cout << "Lens constructed."    << std::endl << std::endl;

  logMessage( std::string("Lens constructed") );


  //////////////////////////////////////////////////////
  //////////////////CONSTRUCT GRID//////////////////////
  //////////////////////////////////////////////////////

  std::cout << "Constructing grid..." << std::endl;

  // Construct grid using physical units of Mpch
  GridMap myGrid( &myLens                    ,
                  userInput.getNgridPoints() ,
                  center                     ,
                  angRange                   ,
                  angRange                   *
                  userInput.getNpixV()       /
                  userInput.getNpixH()       );

  std::cout << "Grid constructed."    << std::endl << std::endl;
  logMessage( std::string("Grid constructed") );



  std::cout << "Generating PixelMaps from grid..." << std::endl;

  calcLensMaps( myGrid,
                         kappaMap ,
                        gamma1Map ,
                        gamma2Map ,
                        invMagMap ,
                         g_tanMap ,
                         g_secMap ,
                          distMap ,
             userInput.getNpixH() ,
             userInput.getNpixV() ,
                         angRange ,
                           center );

  }
  logMessage( std::string("Lens, Grid deallocated") );

  std::cout << "PixelMaps generated" << std::endl << std::endl;


//This section needs work
  ////////////////////////////////////////////////////////////
  ///////////////////Generate sources/////////////////////////
  ////////////////////////////////////////////////////////////


  // Generates random source positions, distances, and errors for the shape measurements
  // Random positions are stored in a 1D array, with pixelmap index
  // Distance is also stored in a 1D array

  std::cout << "Generating sources..." << std::endl;

  double  srcErrArr[ userInput.getNsrc() ]; // Error
  double  srcSCRArr[ userInput.getNsrc() ]; // Sigma Crit
  double  srcDArr  [ userInput.getNsrc() ]; // Redshift

  logMessage( std::string("Allocated src arrays of size: ") + std::to_string((long long) userInput.getNsrc()) );

// Need read error, Z distibution
  for (int i = 0; i < userInput.getNsrc(); ++i ){
    srcErrArr[i] = 0.3;                           // Errors all a temporary 0.3
    srcSCRArr[i] = 1.0;//cosmo.SigmaCrit( myHalo.getZ(), sourceInfo.getZ() ); // Sigma Crit, depends on source Z
  }


  std::cout << "  generating indexes..." << std::endl;


  int indexes[ userInput.getNsrc() ];
  getRandomSourcesIndexes( indexes, userInput );


  logMessage( std::string("Sources placed") );


  std::cout << "  calculating distances..." << std::endl;

  // Determine distances of the sources from center of cluster
  distArrCalc( srcDArr, indexes, &distMap, userInput.getPhysFOV() / userInput.getAngFOV() * 180 / M_PI, userInput.getNsrc() );

  logMessage( std::string("Source distances from center found") );

  std::cout << "Done." << std::endl << std::endl;


  ////////////////////////////////////////////////////////////
  ///////////////////Determine radial average/////////////////
  ///////////////////Bin different source vals////////////////
  ////////////////////////////////////////////////////////////


  // Array to bin distances, RTS, and their errors
  double    distArr[ userInput.getNbins() ],
            gTanArr[ userInput.getNbins() ],
            gErrArr[ userInput.getNbins() ];

  radialDistAverage( distArr, srcDArr, userInput, center );

  logMessage( std::string("Source distances averaged") );

  radialShearAverage( gTanArr, gErrArr, indexes, g_tanMap, srcErrArr, srcDArr, userInput, center );

  logMessage( std::string("Shear values averaged") );

  std::cout << " Sources averaged" << std::endl;;

double m_to_use= 1e13;
double R_to_use=  2.0;
double C_to_use=  7.0;

densProfile myProfile;
myHalo.setRmax(     R_to_use );
myHalo.setC(        C_to_use );
myHalo.setM(        m_to_use );
myProfile.setR_max( R_to_use );
myProfile.setM_enc( m_to_use );
myProfile.setC(     C_to_use );

for ( int i = 0; i < userInput.getNbins(); ++i ){
  distArr[i] = i * userInput.getPhysFOV() / 2 / userInput.getNbins() +1e-3;

  double  srcZ;


  double    SD =    SDNFW( distArr[i], myProfile ); //At radius
  double avgSD = SDAvgNFW( distArr[i], myProfile ); //Average
  double SigCr = cosmo.SigmaCrit( myHalo.getZ(), userInput.getSourceZ() );


  gTanArr[i] = ( avgSD - SD ) / ( SigCr - SD );
  gErrArr[i] = 0.3;

printf("%7.3f %14.5e %14.5e\n",distArr[i], gTanArr[i], gErrArr[i]);
}
printf("\n");
{

densProfile einProfile(0.41);
einProfile.setR_max( R_to_use );
einProfile.setM_enc( m_to_use );
einProfile.setC(     C_to_use );

std::cout << einProfile.getR_s() << std::endl;
generateEinRTS( gTanArr, einProfile, userInput, cosmo.SigmaCrit( myHalo.getZ(), userInput.getSourceZ() ), distArr );
}
exit(0);
  //////////////////////////////////////////////////////////
  ////////////////////////FIT PROFILE///////////////////////
  //////////////////////////////////////////////////////////

  // Attempts to fit the density using the radial averages of distance and RTS


  densProfile nfwProfile, einProfile( 0.2 ); // 0.2 sets profile as Einasto with alpha = 0.2

  nfwProfile.setR_max( myHalo.getRmax() );
  einProfile.setR_max( myHalo.getRmax() );


  rollingFitDensProfile( nfwProfile, myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );
//  fitDensProfile( nfwProfile, myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );

printf("%7.5f %14.4e\n",nfwProfile.getC(), nfwProfile.getM_enc() );


//printf("%12.3e %5.3lf\n",nfwProfile.getM_enc(),nfwProfile.getC());


/*
printf("\n\n\n");
  fitDensProfile( einProfile, lensInfo, userParams, gTanArr, distArr, gErrArr,
                  tempSCRArr, sourceDArr );
printf("%12.3e %5.3lf\n",einProfile.getM_enc(),einProfile.getC());
//*/

  exit(0);
  return 0;
}

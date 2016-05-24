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
//#include "my_utilities.h"
#include "input_functions.h"
//#include "lens_fitter.h"
//#include "pixelmap_functions.h"


// Name of the code log file
std::string logFileName = "";


int main(int arg,char **argv){

  // Initializes the log file, generates logfiles directory
  //  and a file name based on current time
  initLogFile();

  long    seed    = -1827674;


  logMessage( std::string("Seed = ") + std::to_string( (long long) seed ) );


  //////////////////////////////////
  ////////////READ IN///////////////
  //////////////////////////////////


  haloInfo    myHalo;  // Contains all the info on our halo
  userInfo userInput;  // Contains all the info from the user


  std::string   userFile = "lensUserParams.dat" ; // User specified input
  std::string  paramfile = "paramfile"          ; // Glamer input


//Edit these
  // File for glamer read in, contains information used by the library
  // Must contain variables:
  //  field_off
  //  z_source
  //  main_halo_on
  //  pixelmaps_input_file
  //  pixelmaps_on

  InputParams      params     ( paramfile   );
  logMessage( std::string("Read paramfile: ") + paramfile );




  // Determine our FITS file name, will read the header for parameters
  std::string                                    fitsFileName  ;
  params.get(          "pixelmaps_input_file"  , fitsFileName );
  logMessage( std::string("Using FITS file: ") + fitsFileName );



  // Reads info from user file
  readInpFile( userInput, userFile );

  // Read FITS image, populates halo and user inputs
  readFitsHeader( fitsFileName, myHalo, userInput );


  COSMOLOGY cosmo;

  // Read in from lensUserParams
  if ( userInput.getCosmology() == "PLANCK" ) {
    cosmo = Planck1yr;
  } else
  if ( userInput.getCosmology() == "WMAP"   ) {
    cosmo = WMAP5yr;
  } else {
                std::cout << "Unrecognized cosmology: " << userInput.getCosmology() << std::endl;
    logMessage( std::string( "Unrecognized cosmology: ") + userInput.getCosmology() );
    logMessage( "Aborting." );
    exit(1);
  }



exit(0);
//  std::cout << "Using parameter file: " << paramfile << std::endl;
//  std::cout << "Using user lens file: " <<  userFile << std::endl << std::endl;


//  ReadInpFile( userParams, userFile );
  double center[] =    {0,0}; // Center of grid


  // Read in values of cosmology, lens, and source properties from paramfile
  //setCosmoParameters( params, planck ); // Comment if initialized Planck1yr
//  setHaloParameters ( params, planck,   lensInfo);
//  setHaloParameters ( params, planck, sourceInfo, "source");
//  std::cout << std::endl;


//  lensInfo.setRmax( userParams.R_max );

size_t Ninit=1024;
double range=0.3*M_PI/180.;
std::cout<<std::endl<<std::endl;
std::cout<<"Constructing grid"<<std::endl;
Lens lens(params,&seed);
std::cout<<"Lens constructed"<<std::endl;
Grid grid(&lens,Ninit,center,range);
std::cout<<"Grid constructed"<<std::endl;
exit(0);

  ///////////////////////////////////////////////////
  /////////////INITIALIZE NEEDED PARAMETERS//////////
  ///////////////////////////////////////////////////

  /*
      Sets some parameters, and initializes pixelmaps for later use
  */
/*
  srand(seed);                                           // Sets random seed

  omp_set_num_threads( userParams.num_threads );         // For parallelization, default 1

  lensType = userParams.model   + userParams.nbody * 2 + // Either model 1 or nbody 2 or mass map 3
             userParams.massMap * 3;

  if ( lensType!=1 ||  // Can only have 1 lens type
       lensType!=2 ||
       lensType!=3 ){

    std::cout << "Error: LensType = " << lensType << std::endl;
    exit(0);

  }


  double Sigma_crit = sourceInfo.getSigmaCrit();                   // Critital surface density
  double   lensDist =   lensInfo.getAngDist();                     // Ang diam distance to lens
  double  realWidth =   lensInfo.getRealFOV( userParams.angFOV );  // Field of view in Mpc


  PixelMap  kappaMap( center, userParams.N_pixels, userParams.angFOV / userParams.N_pixels );
  PixelMap gamma1Map( center, userParams.N_pixels, userParams.angFOV / userParams.N_pixels );
  PixelMap gamma2Map( center, userParams.N_pixels, userParams.angFOV / userParams.N_pixels );
  PixelMap invMagMap( center, userParams.N_pixels, userParams.angFOV / userParams.N_pixels );
  PixelMap  g_tanMap( center, userParams.N_pixels, userParams.angFOV / userParams.N_pixels );
  PixelMap  g_aziMap( center, userParams.N_pixels, userParams.angFOV / userParams.N_pixels );
  PixelMap   distMap( center, userParams.N_pixels, userParams.angFOV / userParams.N_pixels );
*/



  ///////////////////////////////////////////////////////////
  ///////////////////////CONSTRUCT LENS//////////////////////
  ///////////////////////////////////////////////////////////

  /*

      Constructs a lens based on input from paramfile
      if a model or nbody particle run, need to generate new halo and replace

  */


/*
  //Generates lens from paramfile
  Lens myLens( params, &seed );
//myLens.printMultiLens();

  std::cout <<                    std::endl << std::endl;
  std::cout << "Constructing lens..."       << std::endl;


  //1 model, 2 nbody input, 3 massMap
  if ( lensType==1 ){
    //mass, rmax, redshift, rs, axis ratio, position angle, number of stars
    LensHaloNFW myNFWLens(  lensInfo.getM()     ,
                            lensInfo.getRmax()  ,
                            lensInfo.getZ()     ,
                            lensInfo.getRscale(),
                                             1.0,
                                             0.0,
                                             0  ); //Spherical, normal NFW

    myLens.insertMainHalo( &myNFWLens);
  }

  else
  if ( lensType==2 ){
    //For rotating particles, could be useful for triaxiality
    Point_2d rotation_vector;
    rotation_vector *= 0;
    LensHaloParticles pHalo(  userParams.readFile    ,

                                      lensInfo.getZ(),

                              userParams.N_partSmooth,
                                               planck,
                                      rotation_vector);

    myLens.replaceMainHalos( &pHalo);
  }

  else
  if ( lensType==3 ){
    //Nothing to do if it's SD?
  }

  else{
  std::cout << " Error: lensType="<<lensType<< std::endl << std::endl;
  exit(0);
  }

  std::cout << "Lens constructed."          << std::endl << std::endl;
*/



  /*
  //////////////////////////////////////////////////////
  //////////////////CONSTRUCT GRID//////////////////////
  //////////////////////////////////////////////////////

  Generates a grid, and fills the pixel maps in
  */


/*
  std::cout << "Constructing grid..."       << std::endl;

  Grid myGrid( &myLens, userParams.N_pixels, center, userParams.angFOV);

  std::cout << "Grid constructed."          << std::endl << std::endl;

  calcLensMaps( myGrid,
                         kappaMap,
                        gamma1Map,
                        gamma2Map,
                        invMagMap,
                         g_tanMap,
                         g_aziMap,
                          distMap,
              userParams.N_pixels,
                        realWidth,
                           center);
*/



  /*
  ////////////////////////////////////////////////////////////
  ///////////////////Generate source positions////////////////
  ////////////////////////////////////////////////////////////

  Generates random source positions, distances, and errors for the shape measurements
  Random positions are stored in a 1D array, with pixelmap index
  Distance is also stored in a 1D array
NEED TO CHANGE SOURCE ERRORS
  */



/*
  int indexes[ userParams.N_sources ];
  getRandomSourcesIndexes( indexes, userParams);

////////////////////////////////
////////////////////////////////
////////////////////////////////
//False errors and redshifts, for now
  double  tempErrArr[ userParams.N_sources ];
  double  tempSCRArr[ userParams.N_sources ];
  double  sourceDArr[ userParams.N_sources ];


  for (int i = 0; i < userParams.N_sources; ++i ){
    tempErrArr[i] = 0.3;                                                    //Errors all a temporary 0.3
    tempSCRArr[i] = planck.SigmaCrit( lensInfo.getZ(), sourceInfo.getZ() ); //Just fills an array with Sigma Crit
  }

  //Determine distances of the sources from center of cluster
  distArrCalc( sourceDArr, indexes, userParams, realWidth/userParams.N_pixels, center );
*/




  /*
  ////////////////////////////////////////////////////////////
  ///////////////////Determine radial average/////////////////
  ///////////////////Bin different source vals////////////////
  ////////////////////////////////////////////////////////////

  Generates radial averages of the distances and RTS values of the sources
  */




/*
  //Array to bin distances, RTS, and their errors
  double    distArr[ userParams.N_bins ], gTanArr[ userParams.N_bins ];
  double distErrArr[ userParams.N_bins ], gErrArr[ userParams.N_bins ];

  //Determine radial averages and bin distances and RTS
  radialSourceAverage( distArr, distErrArr, indexes,  distMap, tempErrArr, userParams, center );
  radialSourceAverage( gTanArr,    gErrArr, indexes, g_tanMap, tempErrArr, userParams, center );
*/




  /*
  //////////////////////////////////////////////////////////
  ////////////////////////FIT PROFILE///////////////////////
  //////////////////////////////////////////////////////////

  Attempts to fit the density using the radial averages of distance and RTS
  */



/*
  lensProfile nfwProfile, einProfile( 0.2 ); // 0.2 sets profile as Einasto with alpha = 0.2

  nfwProfile.setR_max( lensInfo.getRmax() );
  einProfile.setR_max( lensInfo.getRmax() );

  fitDensProfile( nfwProfile, lensInfo, userParams, gTanArr, distArr, gErrArr, tempSCRArr, sourceDArr );

printf("%12.3e %5.3lf\n",nfwProfile.getM_enc(),nfwProfile.getC());
*/

/*
printf("\n\n\n");
  fitDensProfile( einProfile, lensInfo, userParams, gTanArr, distArr, gErrArr,
                  tempSCRArr, sourceDArr );
printf("%12.3e %5.3lf\n",einProfile.getM_enc(),einProfile.getC());
//*/

  exit(0);
  return 0;
}

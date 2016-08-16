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
#include "output_functions.h"


// Name of the code log file
std::string logFileName = "";

// Tables for Einasto interpolation
einTable einKappa    ;
einTable einKappaAvg ;

int main(int arg,char **argv){


  // Initializes the log file, generates logfiles directory
  //  and a file name based on current time
  initLogFile();

  // Default values that will not be changing
  long    seed    = time(NULL); //-1827674;

  logMessage( std::string("Seed = ") + std::to_string( (long long) seed ) );

  srand(seed); // Sets random seed


  //////////////////////////////////
  ///////READ IN USERINFO///////////
  //////////////////////////////////



  // User input, most info constant over halos so read in once

  userInfo userInput;                                                 // Contains all the info from the user
  std::string   userFile = "lensUserParams.dat" ;                     // User specified input


  std::cout << "Reading file: "              << std::endl;
  std::cout << "              " << userFile  << std::endl;

  readInpFile( userInput, userFile );                                 // Reads info from user file

  std::cout << "              Done."  << std::endl << std::endl;




  // Read in from lensUserParams, sets the cosmology for glamer

  COSMOLOGY cosmo;

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



  // Sets constant parameters from user input file

  omp_set_num_threads ( userInput.getNthreads()             );  // For parallelization, default 1
  double center[]   = { 0, 0 };                                 // Center of grid, 0,0 aligns grid right


  logMessage( std::string("omp_num_threads = " ) + std::to_string((long double)  userInput.getNthreads()));
  logMessage( std::string("counter         = " ) + std::to_string((long double) center[0]    )
                                                 + std::to_string((long double) center[1]    ));



  // Generates Fox H tables for interpolating over the Einasto profiles

  std::cout << "Reading foxH tables: "                  << std::endl;

  einKappa    = readFoxH( userInput, 1 );
  einKappaAvg = readFoxH( userInput, 2 );

  std::cout << "              Done."  << std::endl << std::endl;





  // Stuff for setting halo for paramfile
  // Loops over FITS files for GLAMER read in

  int halo_index = 0;

  std::string halo_file_name = "";
  std::string paramFileStart = "pixelmaps_input_file    ";




  //////////////////////////////////////
  //////Loop over haloList.dat//////////
  //////////////////////////////////////



  do{



    //////////////////////////////////////
    /////////Generate PARAMFILE///////////
    //////////////////////////////////////


    halo_file_name = getHaloFile( halo_index );           // Get the halo fits file, and write it to paramfile

    if ( halo_file_name == ""  ||
         halo_file_name == " " )  break;                  // Will leave the loop upon EOF

    ++halo_index;

    generateParamfile( paramFileStart + halo_file_name ); // Generates the paramfile for GLAMER from our halo list





    //////////////////////////////////////
    /////////Reads    PARAMFILE///////////
    //////////////////////////////////////


    std::string  paramFile = "paramfile"          ;           // Glamer input

    std::cout << "Reading file: "              << std::endl;
    std::cout << "              " << paramFile << std::endl;

    InputParams      params     ( paramFile   );
    logMessage( std::string("Read paramfile: ") + paramFile );

    std::cout << "              Done."  << std::endl << std::endl;



    // Determine our FITS file name, will read the header for parameters

    std::string                                    fitsFileName  ;
    params.get(          "pixelmaps_input_file"  , fitsFileName );
    logMessage( std::string("Using FITS file: ") + fitsFileName );




    // Read the header of the FITS to determine halo info

    haloInfo    myHalo;

    std::cout << "Reading file: "                  << std::endl;
    std::cout << "              " << fitsFileName  << std::endl;

    // Read FITS image, populates halo and user inputs
    readFitsHeader( fitsFileName, myHalo, userInput );

    std::cout << "              Done."  << std::endl << std::endl;




    // Set parameters based on FITS image

    double   angRange =   userInput.getAngFOV  ()  * M_PI/180  ;  // Angular field of view, in degrees
    double  realWidth =   userInput.getPhysFOV ()              ;  // Field of view in Mpc

    logMessage( std::string("angRange        = " ) + std::to_string((long double)  angRange    ));
    logMessage( std::string("realWidth       = " ) + std::to_string((long double)  realWidth   ));

    userInput.setNgridPoints( userInput.getNpixH() * 2       );   // Number of gridpoints, in reality want it more refined


    // Source redshift moves from glamer parameters to userInfo class
    {
      double srcZ;
      params.get("z_source",srcZ);
      userInput.setSourceZ(srcZ);
    }


    ///////////////////////////////////////////////////
    /////////////INITIALIZE NEEDED MAPS////////////////
    ///////////////////////////////////////////////////



/*
userInput.setNpixH( 9 );
userInput.setNpixV( 9 );
userInput.setNpix( 9*9);
//*/
    std::cout << "Constructing PixelMaps..." << std::endl;

    // PixelMaps we need to keep

    PixelMap  g_totMap( center, userInput.getNpixH(), userInput.getNpixV(), angRange / userInput.getNpixH() );
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
                             g_totMap ,
                              distMap ,
                 userInput.getNpixH() ,
                 userInput.getNpixV() ,
                             angRange ,
                               center );


      std::cout << "PixelMaps generated" << std::endl << std::endl;
    }

    logMessage( std::string("Lens, Grid deallocated") );



    // If the orientation of the halo is known, generate 2D maps of the Reduced shear components

    if ( myHalo.getPhi() != -1 ){

      if ( writeAngRTS( myHalo, userInput , g_tanMap, g_secMap, distMap ) == 2){
      std::cout << "2D maps found" << std::endl << std::endl;
      logMessage(  "2D maps found" );
      }

    } else {
      std::cout << "Halo orientation uknown, ignoring 2D map" << std::endl << std::endl;
      logMessage(  "Halo orientation uknown, ignoring 2D map" );
    }



    ////////////////////////////////////////////////////////////
    ///////////////////Generate sources/////////////////////////
    ////////////////////////////////////////////////////////////


    // Generates random source positions, distances, and errors for the shape measurements
    // Random positions are stored in a 1D array, with pixelmap index
    // Distance is also stored in a 1D array

    std::cout << "Generating sources..." << std::endl;

    double  srcErrArr[ userInput.getNsrc() ]; // Error
    double  srcDArr  [ userInput.getNsrc() ]; // Redshift

    logMessage( std::string("Allocated src arrays of size: ") + std::to_string((long long) userInput.getNsrc()) );

    for (int i = 0; i < userInput.getNsrc(); ++i ){
      srcErrArr[i] = userInput.getShapeNoise();                           // Errors all a temporary 0.3
    }



    // Place in random index locations

    std::cout << "  generating indexes..." << std::endl;

    int indexes[ userInput.getNsrc() ];
    getRandomSourcesIndexes( indexes, userInput );

    logMessage( std::string("Sources placed") );



    // Generate distances based on index

    std::cout << "  calculating distances..." << std::endl;

    // Determine distances of the sources from center of cluster
    distArrCalc( srcDArr, indexes, &distMap, userInput.getPhysFOV() / userInput.getAngFOV() * 180 / M_PI, userInput.getNsrc() );

    logMessage( std::string("Source distances from center found") );

    std::cout << "Done." << std::endl << std::endl;



    ////////////////////////////////////////////////////////////
    ///////////////////To calculate errors//////////////////////
    //////////////Loop over and jacknife sample/////////////////
    ////////////////////////////////////////////////////////////


    densProfile nfwFits[ userInput.getNsrc() + 1 ];
    densProfile nfTFits[ userInput.getNsrc() + 1 ];
    densProfile einFits[ userInput.getNsrc() + 1 ];


                std::cout <<"Generating fits... "  << std::endl;
    logMessage( std::string("Generating fits... ") );


    for ( int  omitIndex = -1; omitIndex < userInput.getNsrc(); ++ omitIndex ) {

      nfwFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      nfTFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      einFits[ omitIndex + 1 ].setR_max( myHalo.getRmax() );
      einFits[ omitIndex + 1 ].setType( 2 );
      nfTFits[ omitIndex + 1 ].setType( 0 );



      ////////////////////////////////////////////////////////////
      ///////////////////Determine radial average/////////////////
      ///////////////////Bin different source vals////////////////
      ////////////////////////////////////////////////////////////



      // Here loop over sources, indicating which index to omit (start at -1)
      // Pass to dist and shear calculators, to indicate ommited index
      if ( omitIndex > -1 ){
                  std::cout <<"  Omitting source " << omitIndex << std::endl;
      logMessage( std::string("  Omitting source ") +
                  std::to_string( (long long) omitIndex ) );
      }

      // Array to bin distances, RTS, and their errors
      double    distArr[ userInput.getNbins() ],
                gTanArr[ userInput.getNbins() ],
                gErrArr[ userInput.getNbins() ];

      radialDistAverage( distArr, srcDArr, userInput, center );

      logMessage( std::string("  Source distances averaged") );

      radialShearAverage( gTanArr, gErrArr, indexes, g_totMap, srcErrArr, srcDArr, userInput, center );

      logMessage( std::string("  Shear values averaged") );

      std::cout << "  Sources averaged" << std::endl;

/*
densProfile testProfile;
testProfile.setM_enc( 1.0e14 );
testProfile.setR_max( 5.0    );
testProfile.setC    ( 5.0    );

nfwFits[ omitIndex + 1 ].setR_max( testProfile.getR_max() );
nfTFits[ omitIndex + 1 ].setR_max( testProfile.getR_max() );
einFits[ omitIndex + 1 ].setR_max( testProfile.getR_max() );
einFits[ omitIndex + 1 ].setType( 2 );
nfTFits[ omitIndex + 1 ].setType( 0 );

for ( int i = 0; i < userInput.getNbins(); ++i ){

  distArr[i] = 1e-3 + i * userInput.getPhysFOV() / userInput.getNbins() / 2.0;
  gErrArr[i] = 0.3;

}
gErrArr[0] = 0;
generateNFWRTS( gTanArr, testProfile, userInput.getNbins(), distArr, cosmo.SigmaCrit( myHalo.getZ(), userInput.getSourceZ() ) );
//*/

      //////////////////////////////////////////////////////////
      ////////////////////////FIT PROFILE///////////////////////
      //////////////////////////////////////////////////////////

      // Attempts to fit the density using the radial averages of distance and RTS

                  std::cout <<"  Calculating NFW fit..." << std::endl;
      logMessage( std::string("  Calculating NFW fit...") );
      rollingFitDensProfile( nfwFits[ omitIndex + 1], myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );


                  std::cout <<"  Calculating NFW trunc fit..." << std::endl;
      logMessage( std::string("  Calculating NFW trunc fit...") );
      rollingFitDensProfile( nfTFits[ omitIndex + 1], myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );


                  std::cout <<"  Calculating EIN fit..." << std::endl;
      logMessage( std::string("  Calculating EIN fit...") );
      rollingFitDensProfile( einFits[ omitIndex + 1], myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );

      std::cout << std::endl ;

    }


    // Calculate the jacknife errors

    double nfwErr[3]; // 0 is C, 1 is log(M), 2 is alpha
    double nfTErr[3];
    double einErr[3];

    jacknife( nfwFits, userInput.getNsrc(), nfwErr );
    jacknife( nfTFits, userInput.getNsrc(), nfTErr );
    jacknife( einFits, userInput.getNsrc(), einErr );


                std::cout <<"Done.              " << std::endl;
    logMessage( std::string("Fitting complete"   ));




    writeProfileFits( userInput, myHalo, einFits[0], nfwFits[0], nfTFits[0], einErr, nfwErr, nfTErr, halo_index );


  } while ( true );

  std::cout << "Reached edge of halo list, exiting" << std::endl;

  logMessage( std::string("Read ") +
              std::to_string( (long long) halo_index ) +
              std::string(" halos"));

/*
Need to check tan/azi uses in the program

Need to work on error introduction
//*/

  exit(0);
  return 0;
}

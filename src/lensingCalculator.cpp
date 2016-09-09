/**
This code is to analyze FITS images output by MAGNMS (Mass mas GeNerator for the Multidark Suite)

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

**/


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

    userInput.setNsrc( -1 ); // Reset each time, as we will change during the program


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

    // On first iteration, check if output exists. If it does, abort
    if ( halo_index == 1 &&
         checkOutputExists( userInput, myHalo ) ) {

                    std::cout << "Output exists for this halo, aborting..." << std::endl;

        logMessage( std::string( "Output exists for this halo, aborting...") );
        exit(0);

    }


    ///////////////////////////////////////////////////
    /////////////INITIALIZE NEEDED MAPS////////////////
    ///////////////////////////////////////////////////



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
                               center ,
             userInput.getImageMass() /
             userInput.getNpix     () /
             cosmo.SigmaCrit(
             myHalo.getZ() ,
             userInput.getSourceZ() ) );


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

    double *srcErrArr ;
    double *srcDArr   ;
    int    *indexes   ;


    double *srcErrArrTemp;
    double *srcDArrTemp  ;
    int    *indexesTemp  ;

    srcErrArrTemp = new double[ userInput.getNsrc() ];
    srcDArrTemp   = new double[ userInput.getNsrc() ];
    indexesTemp   = new int   [ userInput.getNsrc() ];


    logMessage( std::string("Allocated src arrays of size: ") + std::to_string((long long) userInput.getNsrc()) );

    for (int i = 0; i < userInput.getNsrc(); ++i ){
      srcErrArrTemp[i] = userInput.getShapeNoise();                           // Errors all a temporary 0.3
    }



    // Place in random index locations

    std::cout << "  generating indexes..." << std::endl;

    getRandomSourcesIndexes( indexesTemp, userInput );

    logMessage( std::string("Sources placed") );



    // Generate distances based on index

    std::cout << "  calculating distances..." << std::endl;


    // Comb through sources and only keep those outside our halo radius
    {

      int new_Nsrc = 0;

      // Determine distances of the sources from center of cluster, if within Rmax flag as -1
      new_Nsrc = distArrCalc( srcDArrTemp, indexesTemp, &distMap, userInput.getPhysFOV() / userInput.getAngFOV() * 180 / M_PI, userInput.getNsrc(), myHalo.getRmax() );

      logMessage( std::string("Source distances from center found") );


      srcErrArr = new double [ new_Nsrc ];
      srcDArr   = new double [ new_Nsrc ];
      indexes   = new int    [ new_Nsrc ];

      int validCounter = 0;

      for ( int i = 0; i < userInput.getNsrc(); ++i ){
        if ( srcDArrTemp[i] != -1 ){
          srcErrArr[ validCounter ] = srcErrArrTemp[i];
          srcDArr  [ validCounter ] = srcDArrTemp  [i];
          indexes  [ validCounter ] = indexesTemp  [i];
          ++validCounter;
        }
      }
      userInput.setNsrc( new_Nsrc );

    }

    delete[] srcErrArrTemp ;
    delete[] srcDArrTemp   ;
    delete[] indexesTemp   ;


    std::cout << "Done." << std::endl << std::endl;

    // Write the sources to file
    writeDataPoints( userInput, myHalo, srcDArr, indexes, g_totMap );


    ////////////////////////////////////////////////////////////
    ///////////////////To calculate errors//////////////////////
    //////////////Loop over and jacknife sample/////////////////
    ////////////////////////////////////////////////////////////

    int N_jackbins = userInput.getJacknifeBins() * userInput.getJacknifeBins() ;


                std::cout <<"Using " <<     N_jackbins    << " subsets to calculate errors" << std::endl;
    logMessage( std::string   (            "N_jackbins = ") +
                std::to_string( (long long) N_jackbins    ) );


    densProfile nfwFits[ N_jackbins + 1 ];
    densProfile nfTFits[ N_jackbins + 1 ];
    densProfile einFits[ N_jackbins + 1 ];


                std::cout <<"Generating fits... "  << std::endl;
    logMessage( std::string("Generating fits... ") );


    for ( int  omitIndex = -1; omitIndex < N_jackbins ; ++ omitIndex ) {

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
                  std::cout <<"  Omitting subset " << omitIndex << std::endl;
      logMessage( std::string("  Omitting subset ") +
                  std::to_string( (long long) omitIndex ) );
      }

      // Array to bin distances, RTS, and their errors
      double    distArr[ userInput.getNbins() ],
                gTanArr[ userInput.getNbins() ],
                gErrArr[ userInput.getNbins() ];


      radialDistAverage( distArr, srcDArr, userInput, indexes, center, omitIndex );

      logMessage( std::string("  Source distances averaged") );

      radialShearAverage( gTanArr, gErrArr, indexes, g_totMap, srcErrArr, srcDArr, userInput, center, omitIndex );

      logMessage( std::string("  Shear values averaged") );

      std::cout << "  Sources averaged" << std::endl;


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

    jacknife( nfwFits, N_jackbins , nfwErr );
    jacknife( nfTFits, N_jackbins , nfTErr );
    jacknife( einFits, N_jackbins , einErr );


                std::cout <<"Done.              " << std::endl;
    logMessage( std::string("Fitting complete"   ));




    writeProfileFits( userInput, myHalo, einFits[0], nfwFits[0], nfTFits[0], einErr, nfwErr, nfTErr, halo_index );



    delete[] srcErrArr ;
    delete[] srcDArr   ;
    delete[] indexes   ;


  } while ( true );

  std::cout << "Reached edge of halo list, exiting" << std::endl;

  logMessage( std::string("Read ") +
              std::to_string( (long long) halo_index ) +
              std::string(" halos"));


  exit(0);
  return 0;
}

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
/*
einTable foo;
foo = readFoxH( 1 );

for ( int j = 1; j < foo.getA_bins()-1; ++j ){
for ( int i = 0; i < foo.getX_bins()  ; ++i ){

int s1 = (int) (fabs(( foo.getVal( j + 1, i ) - foo.getVal( j     , i ) )) / ( ( foo.getVal( j + 1, i ) - foo.getVal( j     , i ) ) ));
int s2 = (int) (fabs(( foo.getVal( j    , i ) - foo.getVal( j - 1 , i ) )) / ( ( foo.getVal( j    , i ) - foo.getVal( j - 1 , i ) ) ));
//printf("%2i %2i, ",s1,s2);
// NaNs
//if ( foo.getVal( j, i ) != foo.getVal( j, i ) ) {
if ( s1 < 0 && s2 > 0 ) {

  double y0 = foo.getVal( j-1, i );
  double y1 = foo.getVal( j+1, i );

  double x0 = ( j - 1 ) * ( foo.getA_max() - foo.getA_min() ) / foo.getA_bins() + foo.getA_min() ;
  double x  = ( j     ) * ( foo.getA_max() - foo.getA_min() ) / foo.getA_bins() + foo.getA_min() ;
  double x1 = ( j + 1 ) * ( foo.getA_max() - foo.getA_min() ) / foo.getA_bins() + foo.getA_min() ;

printf("%4i %4i %16.12f %16.12f %16.12f\n", j, i, foo.getVal(j-1,i), foo.getVal(j,i), foo.getVal(j+1,i));

  foo.setVal( j, i,

                    y0 + ( y1 - y0 ) * ( x - x0 ) / ( x1 - x0 )

  );

printf("%16.12f\n",foo.getVal(j,i));

}

}
}
exit(0);
FILE *pFile;


pFile = fopen("foxH2012.dat","w");
fprintf(pFile ,"%15.10f %15.10f %3i\n",foo.getX_min(),foo.getX_max(),foo.getX_bins());
fprintf(pFile ,"%15.10f %15.10f %3i\n",foo.getA_min(),foo.getA_max(),foo.getA_bins());


for ( int i = 0; i < foo.getA_bins() ; ++i ){
for ( int j = 0; j < foo.getX_bins() ; ++j ){

fprintf(pFile ,"%15.10f ", foo.getVal( i, j ));

}
fprintf(pFile ,"\n");
}

fclose(pFile);
//std::cout<< interpolateEinRTS( 4, 0.3911, foo )<< std::endl;
//std::cout<<  interpolateEinRTS( , )  <<std::endl;
exit(0);//*/
/*
mpf_set_default_prec(500);

// Generate interpolation table

FILE *pFile;
FILE *pFile2;
pFile = fopen("foxH2012.dat","w");
pFile2= fopen("foxH2123.dat","w");

int  N_bins = 203;
int  a_bins = 203;

double minA =1.3e-1;
double minX =  1e-3;
double maxA =  1   ;
double maxX = 10   ;
double xArr[N_bins];
double aArr[a_bins];
double kArr[N_bins];
double kAvg[N_bins];

fprintf(pFile ,"%15.10f %15.10f %3i\n",minX,maxX,N_bins);
fprintf(pFile ,"%15.10f %15.10f %3i\n",minA,maxA,a_bins);
fprintf(pFile2,"%15.10f %15.10f %3i\n",minX,maxX,N_bins);
fprintf(pFile2,"%15.10f %15.10f %3i\n",minA,maxA,a_bins);


for ( int i = 0; i < N_bins; ++i ){
  xArr[i] = (maxX-minX) / N_bins * i + minX;
}

for ( int i = 0; i < a_bins; ++i ){
  aArr[i] = (maxA-minA) / a_bins * i + minA;
}

// Generate tables

for ( int i = 0; i < a_bins; ++i ){


foxH2012( kArr, kAvg, xArr, N_bins, aArr[i] );


for ( int j = 0; j < N_bins; ++j ){

fprintf(pFile ,"%15.10f ",std::log10(kArr[j]));
fprintf(pFile2,"%15.10f ",std::log10(kAvg[j]));

}
fprintf(pFile ,"\n");
fprintf(pFile2,"\n");

}
fclose( pFile  ) ;
fclose( pFile2 ) ;

exit(0);//*/
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

densProfile myProfile(0.41984538);
//densProfile myProfile;
myProfile.setR_max( myHalo.getRmax() );
myProfile.setC( 5.0 );
myProfile.setM_enc( 1e14 );
//printf("%7.2f %7.2f %14.4e\n",myProfile.getR_max(),myProfile.getC(),myProfile.getM_enc());

for ( int i = 0; i < userInput.getNbins(); ++i ){
  distArr[i] = i * userInput.getPhysFOV() / 2 / userInput.getNbins() +1e-3;

/*
  double    SD =    SDNFW( distArr[i], myProfile ); //At radius
  double avgSD = SDAvgNFW( distArr[i], myProfile ); //Average
  double SigCr = cosmo.SigmaCrit( myHalo.getZ(), userInput.getSourceZ() );


  gTanArr[i] = ( avgSD - SD ) / ( SigCr - SD );
//*/
  gErrArr[i] = 0.3;

//printf("%7.3f %14.5e %14.5e %14.4e\n",distArr[i], gTanArr[i], gErrArr[i], SigCr);
}

generateEinRTS( gTanArr, myProfile, userInput, distArr, cosmo.SigmaCrit( myHalo.getZ(), userInput.getSourceZ() ) );
for ( int i = 0; i < userInput.getNbins(); ++i ){
  printf("%7.3f %14.5e %14.5e\n",distArr[i], gTanArr[i], gErrArr[i]);
}
//*/
printf("\n");
gErrArr[0] = 0;
gErrArr[userInput.getNbins() - 1] = 0;
gErrArr[userInput.getNbins() - 2] = 0;
exit(0);

  //////////////////////////////////////////////////////////
  ////////////////////////FIT PROFILE///////////////////////
  //////////////////////////////////////////////////////////

  // Attempts to fit the density using the radial averages of distance and RTS


  densProfile nfwProfile, einProfile( 0.415843 ); // 0.2 sets profile as Einasto with alpha = 0.2

  nfwProfile.setR_max( myHalo.getRmax() );
  einProfile.setR_max( myHalo.getRmax() );

double gtanNFW[ userInput.getNbins() ];
double gtanEIN[ userInput.getNbins() ];

/*
  rollingFitDensProfile( nfwProfile, myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );
printf("%7.5f %14.4e\n",nfwProfile.getC(), nfwProfile.getM_enc() );
//*/
  rollingFitDensProfile( einProfile, myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );
//  fitDensProfile( nfwProfile, myHalo, userInput, gTanArr, distArr, gErrArr, cosmo );
printf("%7.5f %14.4e %7.5f\n",einProfile.getC(), einProfile.getM_enc(), einProfile.getAlpha() );
printf("               %7.5f %14.4e %7.5f\n",myProfile.getC(), myProfile.getM_enc(),myProfile.getAlpha() );



generateEinRTS( gtanEIN, einProfile, userInput, distArr, cosmo.SigmaCrit( myHalo.getZ(), userInput.getSourceZ() ) );

for ( int i = 0; i < userInput.getNbins(); ++i ){

  double    SD =    SDNFW( distArr[i], nfwProfile ); //At radius
  double avgSD = SDAvgNFW( distArr[i], nfwProfile ); //Average
  double SigCr = cosmo.SigmaCrit( myHalo.getZ(), userInput.getSourceZ() );


printf("%14.4e %14.4e %14.4e\n",gTanArr[i],( avgSD - SD ) / ( SigCr - SD ), gtanEIN[i]);
}
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

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
//My files
#include "astro_constants.h"
#include "lensing_classes.h"
#include "my_utilities.h"
#include "input_functions.h"
#include "lens_fitter.h"
#include "pixelmap_functions.h"



int main(int arg,char **argv){

  long    seed    = -1827674;
  double center[] =    {0,0}; // Center of grid
  int  lensType   =        0; // Will serve as flag, 1 model 2 nbody input

  /*
  ////////////READ IN///////////////
  */

  //Read input of what this source code will be doing, nbins outputfiles etc.
  userInfo userParams;
  std::string userFile  = "lensUserParams.dat";
  std::string paramfile = "paramfile";
  //COSMOLOGY    planck; //Comment if initialized Planck1yr
  COSMOLOGY planck(Planck1yr); //Comment if initialized user read in
  haloInfo   lensInfo, sourceInfo;

  std::cout << "Using parameter file: " << paramfile << std::endl << std::endl;
  std::cout << "Using user lens file: " <<  userFile << std::endl << std::endl;

  ReadInpFile( userParams, userFile );
  InputParams params(paramfile);

  //Read in values of cosmology, lens, and source properties from paramfile
  //setCosmoParameters( params, planck ); //Comment if initialized Planck1yr
  setHaloParameters ( params, planck,   lensInfo);
  setHaloParameters ( params, planck, sourceInfo, "source");
  std::cout << std::endl;

  /*
  /////////////INITIALIZE NEEDED PARAMETERS//////////
  */

  omp_set_num_threads( userParams.num_threads );    //For parallelization, default 1
  lensType = userParams.model + userParams.nbody*2; //Either model 1 or nbody 2
  srand(seed);


  double Sigma_crit = sourceInfo.getSigmaCrit();
  double   lensDist =   lensInfo.getAngDist();
  double  realWidth =   lensInfo.getRealFOV( userParams.angFOV );
  PixelMap  alphaMap(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap alpha1Map(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap alpha2Map(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap  kappaMap(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap  gammaMap(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap gamma1Map(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap gamma2Map(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap invMagMap(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap  g_tanMap(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap  g_aziMap(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);
  PixelMap   distMap(center,userParams.N_pixels,userParams.angFOV/userParams.N_pixels);

  /*
  ///////////////////////CONSTRUCT LENS//////////////////////
  */


//  Lens myLens( &seed, sourceInfo.getZ() );
  Lens myLens(params,&seed);
  myLens.printMultiLens();

/*
  std::cout << "Constructing lens..."       << std::endl;
  //1 model, 2 nbody input. Should almost exclusively use nbody
  if ( lensType==1 ){
    //mass, rmax, redshift, rs, axis ratio, position angle, number of stars
    LensHaloNFW myNFWLens( lensInfo.getM(), lensInfo.getRmax(), lensInfo.getZ(),
                   lensInfo.getRscale(), 1.0, 0.0, 0 ); //Spherical, normal NFW

    myLens.insertMainHalo(&myNFWLens);
//myLens.RevertSourcePlane();
  }
  else
  if ( lensType==2 ){
    //For rotating particles, could be useful for triaxiality
    Point_2d rotation_vector;
    rotation_vector *= 0;
    LensHaloParticles pHalo(userParams.readFile    ,          lensInfo.getZ(),
                            userParams.N_partSmooth, planck, rotation_vector);

//    myLens.replaceMainHalos(&pHalo);
    myLens.insertMainHalo(&pHalo);
    myLens.RevertSourcePlane();
    myLens.printMultiLens();
  }

//exit(0);*/
Point_2d rotation_vector;
rotation_vector *= 0;
LensHaloParticles pHalo(userParams.readFile, lensInfo.getZ(), userParams.N_partSmooth, planck, rotation_vector);

myLens.replaceMainHalos(&pHalo);



  std::cout << "Lens constructed."          << std::endl << std::endl;

  std::cout << "Constructing grid..."       << std::endl;
  Grid myGrid( &myLens, userParams.N_pixels, center, userParams.angFOV);


  std::cout << "Grid constructed."          << std::endl << std::endl;
//exit(0);
  calcLensMaps( myGrid,  alphaMap, alpha1Map, alpha2Map, kappaMap, gammaMap,
                        gamma1Map, gamma2Map, invMagMap, g_tanMap, g_aziMap,
                          distMap,  userParams.N_pixels,realWidth, center);
exit(0);
  /*
  ///////////////////Generate source positions////////////////
  */

  int indexes[userParams.N_sources];
  getRandomSourcesIndexes( indexes, userParams);


  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////
//False errors and redshifts, for now
  double  tempErrArr[userParams.N_sources];
  double  tempSCRArr[userParams.N_sources];
  double  sourceDArr[userParams.N_sources];
  for (int i=0;i<userParams.N_sources;++i){
    tempErrArr[i] = 0.3;
    tempSCRArr[i] = planck.SigmaCrit(lensInfo.getZ(), sourceInfo.getZ());
  }

  //Determine distances of the sources from center of cluster
  distArrCalc( sourceDArr, indexes, userParams,
               realWidth/userParams.N_pixels, center );


  /*
  ///////////////////Determine radial average/////////////////
  ///////////////////Bin different source vals////////////////
  */


  double    distArr[userParams.N_bins], gTanArr[userParams.N_bins];
  double distErrArr[userParams.N_bins], gErrArr[userParams.N_bins];

  radialSourceAverage( distArr, distErrArr, indexes,  distMap,
                    tempErrArr, userParams, center );
  radialSourceAverage( gTanArr,    gErrArr, indexes, g_tanMap,
                    tempErrArr, userParams, center );

printf("GLAMER\n");
for (int i=0;i<userParams.N_bins;++i)
printf("%5.2lf %12.3e\n",distArr[i],gTanArr[i]);
printf("\n");
  /*
    Have radial averages of sources,
    need to compare vs predictions,
    should probably have error, can use dummy

  //*/


lensProfile densProfile;
densProfile.setR_max( lensInfo.getRmax() );
densProfile.setC    (  5.0  );
densProfile.setM_enc( 1e15  );
double gTanArr2[userParams.N_sources];

//for (int i=0;i<userParams.N_sources;++i)
//printf("%12.3e\n",sourceDArr[i]);

//Generate mock NFW for testing purposes
generateNFWRTS( gTanArr2, densProfile, lensInfo, userParams, tempSCRArr, sourceDArr);

printf("\nAnalytic NFW\n");
for (int i=0;i<userParams.N_bins;++i)
printf("%5.2lf %12.3e\n",distArr[i],gTanArr2[i]);

/*
lensProfile densProfile2(2);
densProfile2.setR_max( lensInfo.getRmax() );
densProfile2.setC    (  3.0  );
densProfile2.setM_enc( 1e14  );
densProfile2.setAlpha( 0.12  );

generateEinRTS(gTanArr, densProfile, lensInfo, userParams, tempSCRArr, sourceDArr);
printf("\nAnalytic Ein\n");
for (int i=0;i<userParams.N_bins;++i)
printf("%5.2lf %12.3e\n",distArr[i],gTanArr[i]);
*/

exit(0);
  lensProfile nfwProfile, einProfile(2); // 2 sets profile as Einasto
  nfwProfile.setR_max( lensInfo.getRmax() );
  einProfile.setR_max( lensInfo.getRmax() );

  fitDensProfile( nfwProfile, lensInfo, userParams, gTanArr, distArr, gErrArr,
                  tempSCRArr, sourceDArr );
printf("%12.3e %5.3lf\n",nfwProfile.getM_enc(),nfwProfile.getC());

/*
printf("\n\n\n");
  fitDensProfile( einProfile, lensInfo, userParams, gTanArr, distArr, gErrArr,
                  tempSCRArr, sourceDArr );
printf("%12.3e %5.3lf\n",einProfile.getM_enc(),einProfile.getC());
//*/

  exit(0);
  return 0;
}



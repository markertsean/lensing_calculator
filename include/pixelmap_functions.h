#ifndef PIXELMAP_FUNCTIONS
#define PIXELMAP_FUNCTIONS

#include <lensing_classes.h>
#include <image_processing.h>

void getRandomSourcesIndexes( 
			     int *indexes,  //Array of pixelmap indexes
			     userInfo   u); //Input user information

void radialSourceAverage( 
			 double     *avgArr,  //Output averaged array
			 double     *errArr,  //Radial errors for value we are averaging
			 int       *indexes,  //Index of the source pixels for inpMap
			 PixelMap    inpMap,  //The input map we are radially averaging
			 double     *inpErr,  //The errors in the quantities we are averging
			 userInfo         u,  //Input user information
			 double   center[2]); //Center pixels

void printPixelMap(     
		   PixelMap   &inpMap,  //input pixel map
		   int       N_pixels); //Number of pixels on a side, assumes square map

void findRadialAverage(
		       PixelMap      &inpMap,   //input map to average
		       int          N_pixels,   //Number of pixels on a side, assumes square map
		       int            N_bins,   //Number of bins to average to
		       double   *avgInpArray,   //Output array to store average in
		       double  *distInpArray,   //Outpput array to store distances in
		       double       realSize,   //Real size on the plane of the FOV
		       double      center[2]);  //Center coordinates


void calcLensMaps(
		  Grid        &inpGrid,  //GLAMER grid to calc values on
		  PixelMap   &alphaMap,  //Various input maps, only need g_tan/azi out
		  PixelMap  &alpha1Map,
		  PixelMap  &alpha2Map,
		  PixelMap   &kappaMap,
		  PixelMap   &gammaMap,
		  PixelMap  &gamma1Map,
		  PixelMap  &gamma2Map,
		  PixelMap  &invMagMap,
		  PixelMap   &g_tanMap,
		  PixelMap   &g_aziMap,
		  PixelMap    &distMap,
		  int         N_pixels,  //Number of pixels on a side, assumes square map
		  double      realSize,  //Real width on the 2D sky plane
		  double     center[2]); //Center location of halo


void calcMapsFromKappa(
		       PixelMap  &kappaMap, //Various input maps, some overwritten
		       PixelMap    &phiMap,
		       PixelMap   &massMap,
		       PixelMap   &distMap,
		       PixelMap  &sigmaMap,
		       double    center[2],  //Center location of halo
		       userInfo    inpInfo,  //User information
		       double    sigmaCrit,  //Critical surface density to use for whole map
		       double     realSize); //Real size on the 2D sky plane


void sigmaMapCalc( 
		  PixelMap   &inpMap, 
		  PixelMap     &kMap, 
		  int       N_pixels, 
		  double  Sigma_crit);


void   phiMapCalc( 
		  PixelMap   &phiMap, 
		  PixelMap &sigmaMap, 
		  int       N_pixels, 
		  double    realSize);


void   magMapCalc( 
		  PixelMap   &magMap, 
		  PixelMap     &kMap, 
		  PixelMap     &gMap, 
		  int       N_pixels);


void  massMapCalc( 
		  PixelMap  &massMap, 
		  PixelMap   &phiMap, 
		  PixelMap  &distMap, 
		  int       N_pixels);


void  distMapCalc( 
		  PixelMap  &distMap, 
		  int       N_pixels, 
		  double     inpSize,
		  double   center[2]);


void  distArrCalc(
		  double *sourceDistArr,
		  int          *indexes, 
		  userInfo            u,
                  double          scale,
		  double      center[2]);


#endif

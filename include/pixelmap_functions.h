#ifndef PIXELMAP_FUNCTIONS
#define PIXELMAP_FUNCTIONS

#include <lensing_classes.h>
#include <gridmap.h>
#include <image_processing.h>


void radialDistAverage( double       *avgArr , // Array to overwrite
                        double    *distances , // Array of source distances
                        userInfo           u , // User parameters
                        double     center[2] );

void radialShearAverage( double      *avgArr ,  // Array to overwrite
                         double      *errArr ,  // Error array to overwrite
                         int        *indexes ,  // Indexes of sources
                         PixelMap     inpMap ,  // Array to sample from
                         double      *errors ,  // Errors to use in weighting
                         double        *dist ,  // Distances of the halos
                         userInfo          u ,
                         double    center[2] );



void  distArrCalc(  double *sourceDistArr ,
                    int          *indexes ,
                    PixelMap     *distMap ,
                    int         N_srouces );



void  distMapCalc(  PixelMap  &distMap ,
                    int      N_pixelsH ,
                    int      N_pixelsV ,
                    double     inpSize ,
                    double   center[2] );



void calcLensMaps(  GridMap     &inpGrid ,  //GLAMER grid to calc values on
                    PixelMap   &kappaMap ,
                    PixelMap  &gamma1Map ,
                    PixelMap  &gamma2Map ,
                    PixelMap  &invMagMap ,
                    PixelMap   &g_tanMap ,
                    PixelMap   &g_aziMap ,
                    PixelMap    &distMap ,
                    int       N_pixels_h ,  //Number of pixels on a side
                    int       N_pixels_v ,  //Number of pixels on a side
                    double      realSize ,  //Real width on the 2D sky plane
                    double     center[2] ); //Center location of halo


void printPixelMap( PixelMap   &inpMap ,  //input pixel map
                    int       N_pixels ); //Number of pixels on a side, assumes square map

void printPixelMap( PixelMap   &inpMap   ,  //input pixel map
                    int       N_pixels_h ,
                    int       N_pixels_v ); //Number of pixels on a side


void getRandomSourcesIndexes( int      *indexes ,  //Array of pixelmap indexes
                              userInfo        u ); //Input user information



#endif

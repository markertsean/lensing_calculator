#ifndef PIXELMAP_FUNCTIONS
#define PIXELMAP_FUNCTIONS

#include <lensing_classes.h>
#include <gridmap.h>
#include <image_processing.h>
#include "my_utilities.h"


void radialDistAverage( double       *avgArr , // Array to overwrite
                        double    *distances , // Array of source distances
                        userInfo           u , // User parameters
                        double     center[2] ,
                        int ignoreIndex = -1 );// Index to ignore, for errors

void radialShearAverage( double      *avgArr ,  // Array to overwrite
                         double      *errArr ,  // Error array to overwrite
                         int        *indexes ,  // Indexes of sources
                         PixelMap     inpMap ,  // Array to sample from
                         double      *errors ,  // Errors to use in weighting
                         double        *dist ,  // Distances of the halos
                         userInfo          u ,
                         double    center[2] ,
                         int ignoreIndex = 1 ); // Index to ignore, for errors


void distArrCalc( double *sourceDistArr ,  // Array to overwrite
                  int          *indexes ,  // Source locations
                  PixelMap     *distMap ,  // Map of distances
                  double          scale ,  // Mpc/rad conversion
                  int         N_sources ); // Number of source



void  distMapCalc(  PixelMap  &distMap ,
                    int      N_pixelsH ,
                    int      N_pixelsV ,
                    double     inpSize ,
                    double   center[2] );



double calcLensMaps(  GridMap     &inpGrid ,  //GLAMER grid to calc values on
                      PixelMap   &kappaMap ,
                      PixelMap  &gamma1Map ,
                      PixelMap  &gamma2Map ,
                      PixelMap  &invMagMap ,
                      PixelMap   &g_tanMap ,
                      PixelMap   &g_aziMap ,
                      PixelMap   &g_totMap ,
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

// Box-Muller transformation to provide gaussian distribution
double gaussErr( userInfo     u ,
                 int       Ngal ); // Number of galaxies


void jacknife( densProfile  *profile ,
               int         N_samples ,
               double       * errArr );


#endif

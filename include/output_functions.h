#ifndef OUTPUT_FUNCTIONS
#define OUTPUT_FUNCTIONS

#include <gridmap.h>
#include <image_processing.h>

#include "lensing_classes.h"



int  writeAngRTS( haloInfo   & h ,
                  userInfo     u ,
                  PixelMap  gTan ,
                  PixelMap  gSec ,
                  PixelMap  dMap );


void generateParamfile( std::string haloName );

std::string getHaloFile( int index );


void writeProfileFits( userInfo        u ,   // User input
                       haloInfo        h ,   // Info on our halo
                       densProfile   ein ,   // Einasto   density profile
                       densProfile   nfw ,   // NFW Full  density profile
                       densProfile   nfT ,   // NFW trunc density profile
                       double    *einErr ,   // Einasto   errors
                       double    *nfwErr ,   // NFW Full  errors
                       double    *nfTErr ,   // NFW trunc errors
                       int       haloNum );  // How many times we've written, first time we need to write halo info

// Writes source data points for stacking
void writeDataPoints(  userInfo        u ,
                       haloInfo        h ,
                       double      *dArr ,
                       int      *indexes ,
                       PixelMap    &gMap );

bool checkOutputExists( userInfo       u , // If output files exist before first run, abort
                        haloInfo       h );


#endif // OUTPUT_FUNCTIONS

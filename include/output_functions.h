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
                       densProfile   ein ,   // Einasto  density profile
                       densProfile   nfw ,   // NFW Full density profile
                       double    *einErr ,   // Einasto  errors
                       double    *nfwErr ,   // NFW Full errors
                       int       haloNum );  // How many times we've written, first time we need to write halo info



#endif // OUTPUT_FUNCTIONS

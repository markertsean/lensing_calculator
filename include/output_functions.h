#ifndef OUTPUT_FUNCTIONS
#define OUTPUT_FUNCTIONS

#include <gridmap.h>
#include <image_processing.h>

#include "lensing_classes.h"



int  writeAngRTS( haloInfo     h ,
                  userInfo     u ,
                  PixelMap  gTan ,
                  PixelMap  gSec ,
                  PixelMap  dMap );


void generateParamfile( std::string haloName );

std::string getHaloFile( int index );


#endif // OUTPUT_FUNCTIONS

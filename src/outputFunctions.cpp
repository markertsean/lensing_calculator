#include "output_functions.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <vector>


void checkDir( std::string dirName ){
    struct stat sb;
    if ( stat( dirName.c_str(), &sb ) !=0 ){
      char str[100];
      sprintf( str, "mkdir %s", dirName.c_str());
      system( str );
      logMessage( std::string("Wrote directory: ") + dirName );
    }
}



void writeAngRTS( haloInfo     h ,
                  userInfo     u ,
                  PixelMap  gTan ,
                  PixelMap  gSec ,
                  PixelMap  dMap ){


  // If output directory does not exist, create one
  checkDir( u.getOutputPath() );

//h.setTheta(0.5*M_PI);
printf("Theta: %10.6f   cos: %10.6f   sin: %10.6f\n",h.getTheta()/M_PI, cos(h.getTheta()), sin(h.getTheta()) );
//h.setPhi  (0.0*M_PI);
printf("Phi  : %10.6f   tan: %10.6f   sin: %10.6f\n",h.getPhi  ()/M_PI, tan(h.getPhi  ()), sin(h.getPhi  ()));

  double x_0 = sin( h.getTheta() ) * sin( h.getPhi() );
  double y_0 = cos( h.getPhi  () ) ;

// Open files, write info on halo and image
// gTan table, gSec table, g table
  double zInc    = acos( cos( h.getTheta() ) *  // Inclination out of the page, 0 directly out
                         sin( h.getPhi  () ) );

  double alpha_0 = atan2( x_0, y_0 );

  if ( alpha_0 < 0 ) alpha_0 += 2*M_PI;

  double gTan_binned[ u.getNbins_A2D() ][ u.getNbins_R2D() ]; // Our bins
  double gSec_binned[ u.getNbins_A2D() ][ u.getNbins_R2D() ];
  double gTot_binned[ u.getNbins_A2D() ][ u.getNbins_R2D() ];
  double NTot_binned[ u.getNbins_A2D() ][ u.getNbins_R2D() ];


  for ( int i = 0; i < u.getNbins_A2D(); ++ i ){ // Initialize to 0
  for ( int j = 0; j < u.getNbins_R2D(); ++ j ){
    gTan_binned[i][j] = 0;
    gSec_binned[i][j] = 0;
    gTot_binned[i][j] = 0;
    NTot_binned[i][j] = 0;
  }
  }


  // Boundaries for bins

  double rMin = 0;
  double rMax = u.getPhysFOV() / 2.;
  double aMin = 0;
  double aMax = 2 * M_PI;

  double dScale = u.getPhysFOV() / u.getAngFOV() * 180 / M_PI;

printf("Z_inc = %10.6f\n",zInc/M_PI);
printf("alpha = %10.6f\n",alpha_0/M_PI);

  // Bin the data
  for ( int i = 0; i < u.getNpix(); ++i ){

    double x =    i % u.getNpixH()   + 0.5 - u.getNpixH() / 2.0;
    double y = -( i / u.getNpixH() ) - 0.5 + u.getNpixV() / 2.0;

    double alpha = atan2( x, y ) - alpha_0;

    while ( alpha < 0 ){
      alpha += 2*M_PI;
    }


    int rBin = std::min(
               std::max(

                     (int) ( dScale * dMap[ i ] / ( rMax-rMin ) * u.getNbins_R2D() )

                                                  , 0                )
                                                  , u.getNbins_R2D() );

    int aBin = std::min(
               std::max(

                          (int) ( alpha        / ( aMax-aMin ) * u.getNbins_A2D() )

                                                  , 0                )
                                                  , u.getNbins_A2D() );
    gTan_binned[aBin][rBin] += gTan[i];
    gSec_binned[aBin][rBin] += gSec[i];
    gTot_binned[aBin][rBin] += sqrt( gTan[i] * gTan[i] + gSec[i] * gSec[i] );
    NTot_binned[aBin][rBin] += 1;
  }

  for ( int i = 0; i < u.getNbins_A2D(); ++ i ){ // Take the average
  for ( int j = 0; j < u.getNbins_R2D(); ++ j ){
    gTan_binned[i][j] /= NTot_binned[i][j];
    gSec_binned[i][j] /= NTot_binned[i][j];
    gTot_binned[i][j] /= NTot_binned[i][j];
  }
  }

  for ( int i = u.getNbins_R2D()-1; i > -1; -- i ){ // Take the average
  for ( int j = 0; j < u.getNbins_A2D(); ++ j ){
    gTot_binned[j][i] /= NTot_binned[j][i];
printf("%10.6f",log10(gTot_binned[j][i]));
  }
printf("\n");
  }


}



/*
  double posArr[2]= { 0, 0 }; // Pixel position
  double phi      =   0;      // Position angle


  for (int i=0;i<N_pixels_v;++i){
    posArr[1] = (-i - 0.5 + N_pixels_v/2.0)-center[1];

  for (int j=0;j<N_pixels_h;++j){
    posArr[0] = ( j + 0.5 - N_pixels_h/2.0)-center[0];

    int k = j+i*N_pixels_h;
    phi = atan2(posArr[1],posArr[0]);
    //gamma1  <0 |   >0 -
    //gamma2  <0 \   >0 /
    //g_tan = -g1*cos - g2*sin
    //g_sec =  g1*sin - g2*cos

//*/




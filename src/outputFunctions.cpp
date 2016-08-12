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

bool checkFile( char dirName[] ){
    struct stat sb;
    if ( stat( dirName, &sb ) !=0 ){
      return false;
    }
      return true;
}


// Reads the halo fits file from a file, to generate file for GLAMER
// If done with file, returns ""
std::string getHaloFile ( int index ){

  std::string   inputLine;
  std::ifstream halo_file("haloList.dat" );

  int counter = 0;

  if ( halo_file.is_open() ){
    while( getline( halo_file, inputLine ) ){
      if ( counter == index ) break;
      if ( halo_file.eof()  ) break;
      ++counter;
    }
  } else{
    std::cout << "Code requires files \"haloList.dat\" in executing directory" << std::endl;
    logMessage("Couldn't open haloList");
    exit(1);
  }

  halo_file.close();

  if ( counter != index ) inputLine = "";

  return inputLine;
}


// Generates paramfile for GLAMER from param_build and haloList.dat
// param_build contains bulk of file, copied directly
// haloList.dat contains list of halo fits files to use
void generateParamfile( std::string haloName ){

  std::string   inputLine;
  std::ifstream param_build("param_build" );
  std::ofstream param_file ("paramfile");

  if ( param_build.is_open() && param_file.is_open() ){

    while( getline( param_build, inputLine ) ){
      param_file << inputLine << std::endl;
    }

  } else{
    std::cout << "Code requires files \"param_build\" in executing directory" << std::endl;
    logMessage("Couldn't open param_build");
    exit(1);
  }

      param_file << haloName << std::endl;

  param_build.close();
  param_file.close();

}



int  writeAngRTS( haloInfo     h ,
                  userInfo     u ,
                  PixelMap  gTan ,
                  PixelMap  gSec ,
                  PixelMap  dMap ){


  // If output directory does not exist, create one
  checkDir( u.getOutputPath() );


  double integ = u.getIntegLength(); // Need for file name

  if ( integ == -1 ){                // If sphere, just mark as 0 integ length
    integ = 0.0;
  }

  char     fileName[100];
  sprintf( fileName, "%sHalo_%010li_%06.1f_AngDistMap.dat", u.getOutputPath().c_str(), h.getID(), integ);


  double x_0 = sin( h.getTheta() ) * sin( h.getPhi() );
  double y_0 = cos( h.getPhi  () ) ;

  double zInc    = acos( cos( h.getTheta() ) *  // Inclination out of the page, 0 directly out
                         sin( h.getPhi  () ) );

  double alpha_0 = atan2( x_0, y_0 );           // Orientation of halo on xy plane, 0 directly up

  if ( alpha_0 < 0 ) alpha_0 += 2*M_PI;

  h.setAlpha( alpha_0 );                        // Store orientation as halo info
  h.setGamma( zInc    );

  if ( checkFile( fileName ) )                  // If file exists, don't bother writing a new one
    return 2;


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

  // Bin the data
  for ( int i = 0; i < u.getNpix(); ++i ){

    double x =    i % u.getNpixH()   + 0.5 - u.getNpixH() / 2.0; // Positions of the pixels
    double y = -( i / u.getNpixH() ) - 0.5 + u.getNpixV() / 2.0;

    double alpha = atan2( x, y ) - alpha_0;                      // Angle relative to halo orientation

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


  sprintf( fileName, "%sHalo_%010li_%06.1f_AngDistMap.dat", u.getOutputPath().c_str(), h.getID(), integ);

  FILE *pFile;

  pFile = fopen( fileName, "w" );

  fprintf( pFile, "M       %14.6e\n", h.getM    () );
  fprintf( pFile, "C       %14.6e\n", h.getPhi  () );
  fprintf( pFile, "R_max   %14.6e\n", h.getPhi  () );


  fprintf( pFile, "phi     %14.6e\n", h.getPhi  () );
  fprintf( pFile, "theta   %14.6e\n", h.getTheta() );

  fprintf( pFile, "alpha_0  %14.6e\n", alpha_0     );
  fprintf( pFile, "z_inc    %14.6e\n", zInc        );

  fprintf( pFile, "a %10.6f %10.6f %4i\n", aMin, aMax, u.getNbins_A2D() );
  fprintf( pFile, "r %10.6f %10.6f %4i\n", rMin, rMax, u.getNbins_R2D() );

  for ( int i = 0; i < u.getNbins_R2D(); ++ i ){ // Prints total shear to file
  for ( int j = 0; j < u.getNbins_A2D(); ++ j ){
    fprintf( pFile, "%14.6e",gTot_binned[j][i]);
  }
    fprintf( pFile, "\n");
  }
    fprintf( pFile, "\n");


  for ( int i = 0; i < u.getNbins_R2D(); ++ i ){ // Prints tangential shear to file
  for ( int j = 0; j < u.getNbins_A2D(); ++ j ){
    fprintf( pFile, "%14.6e",gTan_binned[j][i]);
  }
    fprintf( pFile, "\n");
  }
    fprintf( pFile, "\n");


  for ( int i = 0; i < u.getNbins_R2D(); ++ i ){ // Prints secantial shear to file
  for ( int j = 0; j < u.getNbins_A2D(); ++ j ){
    fprintf( pFile, "%14.6e",gSec_binned[j][i]);
  }
    fprintf( pFile, "\n");
  }

    fprintf( pFile, "\n");

  fclose( pFile );

  std::cout << "Wrote 2D map file: " << fileName << std::endl << std::endl;
  logMessage(  "Wrote 2D map file" );

  return 1;

}





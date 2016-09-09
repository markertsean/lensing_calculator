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

// Returns true if file exists
bool checkFile( char dirName[] ){
    struct stat sb;
    if ( stat( dirName, &sb ) !=0 ){
      return false;
    }
      return true;
}


bool checkOutputExists( userInfo       u , // If output files exist before first run, abort
                        haloInfo       h ){

  char     fileName[100];
  sprintf( fileName, "%sHalo_%010li_densFits.dat", u.getOutputPath().c_str(), h.getID() );

  return checkFile( fileName );

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



int  writeAngRTS( haloInfo   & h ,
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
    if ( gTan[i] == gTan[i] &&
         gSec[i] == gSec[i] ){

      gTan_binned[aBin][rBin] += gTan[i];
      gSec_binned[aBin][rBin] += gSec[i];
      gTot_binned[aBin][rBin] += sqrt( gTan[i] * gTan[i] + gSec[i] * gSec[i] );
      NTot_binned[aBin][rBin] += 1;

    }

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


  fprintf( pFile, "ID      %10li\n" , h.getID   () );
  fprintf( pFile, "M       %14.6e\n", h.getM    () );
  fprintf( pFile, "C       %14.6e\n", h.getC    () );
  fprintf( pFile, "R_max   %14.6e\n", h.getRmax () );

  fprintf( pFile, "integ   %14.6e\n", integ        );
  fprintf( pFile, "M_img   %14.6e\n", u.getImageMass() );

  fprintf( pFile , "b/a    %10.6f\n", h.getBA   () );
  fprintf( pFile , "c/a    %10.6f\n", h.getCA   () );

  fprintf( pFile, "phi     %14.6e\n", h.getPhi  () );
  fprintf( pFile, "theta   %14.6e\n", h.getTheta() );

  fprintf( pFile, "alpha   %14.6e\n", alpha_0     );
  fprintf( pFile, "gamma   %14.6e\n", zInc        );

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




void writeProfileFits( userInfo        u ,   // User input
                       haloInfo        h ,   // Info on our halo
                       densProfile   ein ,   // Einasto   density profile
                       densProfile   nfw ,   // NFW Full  density profile
                       densProfile   nfT ,   // NFW trunc density profile
                       double    *einErr ,   // Einasto   errors
                       double    *nfwErr ,   // NFW Full  errors
                       double    *nfTErr ,   // NFW trunc errors
                       int       haloNum ){  // How many times we've written, first time we need to write halo info


  checkDir( u.getOutputPath() );

  double integ = u.getIntegLength(); // Need for file name

  if ( integ == -1 ){                // If sphere, just mark as 0 integ length
    integ = 0.0;
  }

  char     fileName[100];
  sprintf( fileName, "%sHalo_%010li_densFits.dat", u.getOutputPath().c_str(), h.getID() );



  if ( haloNum == 1 ){ // First time through, write the general halo info
    FILE *pFile;

    pFile = fopen( fileName, "w" );

    fprintf( pFile , "ID          %10li\n" , h.getID   () );
    fprintf( pFile , "M           %14.6e\n", h.getM    () );
    fprintf( pFile , "C           %10.6f\n", h.getC    () );
    fprintf( pFile , "R_max       %10.6f\n", h.getRmax () );
    fprintf( pFile , "Z           %10.6f\n", h.getZ    () );
    fprintf( pFile , "b/a         %10.6f\n", h.getBA   () );
    fprintf( pFile , "c/a         %10.6f\n", h.getCA   () );
    fprintf( pFile , "phi         %10.6f\n", h.getPhi  () );
    fprintf( pFile , "theta       %10.6f\n", h.getTheta() );
    fprintf( pFile , "alpha       %10.6f\n", h.getAlpha() );
    fprintf( pFile , "gamma       %10.6f\n", h.getGamma() );

    fprintf( pFile , "FOV         %10.6f\n", u.getPhysFOV() );
    fprintf( pFile , "N_pixH      %10i\n"  , u.getNpixH  () );
    fprintf( pFile , "N_pixV      %10i\n"  , u.getNpixV  () );

    fprintf( pFile , "N_src       %10i\n"  , u.getNsrc      () );
    fprintf( pFile , "Z_src       %10.6f\n", u.getSourceZ   () );
    fprintf( pFile , "sigma_shape %10.6f\n", u.getShapeNoise() );

    fclose(  pFile );
  }

  FILE *pFile;

  pFile = fopen( fileName, "a+" );



  fprintf( pFile , "IntegLength %10.6f\n"        , u.getIntegLength() );
  fprintf( pFile , "ImageMass   %14.6e\n"        , u.getImageMass  () );
  fprintf( pFile , "NFW_Full %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" , log10( nfw.getM_enc() ), nfw.getC(),           -1.0, nfwErr[1], nfwErr[0],      -1.0);
  fprintf( pFile , "NFW_Trnc %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" , log10( nfT.getM_enc() ), nfT.getC(),           -1.0, nfTErr[1], nfTErr[0],      -1.0);
  fprintf( pFile , "Ein      %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n" , log10( ein.getM_enc() ), ein.getC(), ein.getAlpha(), einErr[1], einErr[0], einErr[2]);


  fclose( pFile );

  std::cout << "Appended file: " << fileName << std::endl << std::endl;

}



// Writes source data points for stacking
void writeDataPoints(  userInfo        u ,
                       haloInfo        h ,
                       double      *dArr ,
                       int      *indexes ,
                       PixelMap    &gMap ) {


  checkDir( u.getOutputPath() );

  double integ = u.getIntegLength(); // Need for file name

  if ( integ == -1 ){                // If sphere, just mark as 0 integ length
    integ = 0.0;
  }

  char     fileName[100];
  sprintf( fileName, "%sHalo_%010li_%06.1f_Sources.dat", u.getOutputPath().c_str(), h.getID(), integ);


  FILE *pFile;

  pFile = fopen( fileName, "w" );

  fprintf( pFile , "ID          %10li\n" , h.getID   () );
  fprintf( pFile , "M           %14.6e\n", h.getM    () );
  fprintf( pFile , "C           %10.6f\n", h.getC    () );
  fprintf( pFile , "R_max       %10.6f\n", h.getRmax () );
  fprintf( pFile , "Z           %10.6f\n", h.getZ    () );
  fprintf( pFile , "b/a         %10.6f\n", h.getBA   () );
  fprintf( pFile , "c/a         %10.6f\n", h.getCA   () );
  fprintf( pFile , "phi         %10.6f\n", h.getPhi  () );
  fprintf( pFile , "theta       %10.6f\n", h.getTheta() );
  fprintf( pFile , "alpha       %10.6f\n", h.getAlpha() );
  fprintf( pFile , "gamma       %10.6f\n", h.getGamma() );

  fprintf( pFile , "IntegLength %10.6f\n", u.getIntegLength() );
  fprintf( pFile , "IntegMass   %14.8e\n", u.getImageMass() );
  fprintf( pFile , "FOV         %10.6f\n", u.getPhysFOV() );
  fprintf( pFile , "N_pixH      %10i\n"  , u.getNpixH  () );
  fprintf( pFile , "N_pixV      %10i\n"  , u.getNpixV  () );

  fprintf( pFile , "N_src       %10i\n"  , u.getNsrc      () );
  fprintf( pFile , "Z_src       %10.6f\n", u.getSourceZ   () );
  fprintf( pFile , "sigma_shape %10.6f\n", u.getShapeNoise() );




  int posArr[] = { 0, 0 };

  // For each source, find which pixel it's in,
  //  and use nearby pixels to generate an average
  //  shear value to use
  for(int i = 0; i < u.getNsrc(); ++i ){


    // 1 is y, or row
    // 0 is x, or columns
    // Pixel numbers, from top corner
    posArr[0] = indexes[i] % u.getNpixH();
    posArr[1] = indexes[i] / u.getNpixH();

    // Locate indexes for source averaging
    int startXIndex = std::max( (int) ( posArr[0] - u.getSourceRadius() ),            0 );
    int   endXIndex = std::min( (int) ( posArr[0] + u.getSourceRadius() ), u.getNpixH() );

    int startYIndex = std::max( (int) ( posArr[1] - u.getSourceRadius() ), 0 );
    int   endYIndex = std::min( (int) ( posArr[1] + u.getSourceRadius() ), u.getNpixV() );

    double avgVal = 0;
    int    n_srcs = 0;

    // Average over the nearby pixels for a given source
    for ( int j = startXIndex; j <= endXIndex; ++j ){
    for ( int k = startYIndex; k <= endYIndex; ++k ){

      avgVal += gMap.getValue( k * u.getNpixH() + j );
      n_srcs += 1;

    }
    }

  fprintf( pFile , "%10.6f %14.8e\n", dArr[i], avgVal / n_srcs );

  }



  fclose(  pFile );

  std::cout << "Wrote file: " << fileName << std::endl;


}


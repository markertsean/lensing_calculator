#include <slsimlib.h>
#include <image_processing.h>
#include <gridmap.h>
#include <grid_maintenance.h>
#include <astro_constants.h>
#include <lensing_classes.h>
#include <pixelmap_functions.h>



//
//  Radially average pixelmap values for given source positions
//
void radialShearAverage( double      *avgArr ,  // Array to overwrite
                         double      *errArr ,  // Error array to overwrite
                         int        *indexes ,  // Indexes of sources
                         PixelMap     inpMap ,  // Array to sample from
                         double      *errors ,  // Errors to use in weighting
                         double        *dist ,  // Distances of the halos in Mpc
                         userInfo          u ,
                         double    center[2] ,
                         int     ignoreIndex ){

  int    iBin;
  int posArr[2]={0,0};

  int N_countsArr[ u.getNbins() ];

  // Make sure average is 0 to start
  for(int i = 0; i < u.getNbins(); ++i ){
         avgArr[i] = 0;
         errArr[i] = 0;
    N_countsArr[i] = 0;
  }

  // For each source, find which bin it's in,
  //  and use nearby pixels to generate an average
  //  shear value to use
  for(int i = 0; i < u.getNsrc(); ++i ){

    // Ignore for error analysis
    if ( i != ignoreIndex ){

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

        avgVal += inpMap.getValue( k * u.getNpixH() + j );
        n_srcs += 1;

      }
      }


      // Distance of radians from center, converted to "bin" units
      iBin = round( dist[i] / u.getPhysFOV() * 2 * u.getNbins() );


           avgArr[iBin] += avgVal / n_srcs / ( errors[i]*errors[i] ); // Weighted average sum
           errArr[iBin] += 1               / ( errors[i]*errors[i] ); // inverse of the sum of variance
      N_countsArr[iBin] += 1;
    }
  }


  // Calculate weighted average of each bin
  for(int i = 0; i < u.getNbins(); ++i){

    if (N_countsArr[i]==0){ // If nothing in bin, mark as -1
      avgArr[i] = -1.0;
      errArr[i] =  0.0;
    }

    else{

      errArr[i] = sqrt( 1   /   errArr[i]             );
      avgArr[i] = avgArr[i] * ( errArr[i] * errArr[i] ) ;//* gaussErr( u, N_countsArr[i] );

    }
  }

}



//
//  Radially average pixelmap values for given source positions
//
void radialDistAverage( double       *avgArr ,  // Array to overwrite
                        double    *distances ,  // Distance array of sources
                        userInfo           u ,
                        double     center[2] ,
                        int      ignoreIndex ){

  int    iBin;
  double dist;
  double posArr[2]={0,0};

  int N_countsArr[ u.getNbins() ];

  // Make sure average is 0 to start
  for(int i = 0; i < u.getNbins(); ++i ){
         avgArr[i] = 0;
    N_countsArr[i] = 0;
  }

  // For each source, find which bin it's in
  for(int i = 0; i < u.getNsrc(); ++i ){

    // For errors
    if ( i != ignoreIndex ) {
      // Distance of pixels from center, converted to "bin" units
      iBin = round( distances[i] / u.getPhysFOV() * 2 * u.getNbins() );

           avgArr[iBin] += distances[i];
      N_countsArr[iBin] += 1;
    }
  }

  // Perform the averaging
  for(int i = 0; i < u.getNbins(); ++i){

    if (N_countsArr[i]==0){ // If nothing in bin, mark as -1
      avgArr[i] = -1.0;
    }

    else{

      avgArr[i] = avgArr[i] / N_countsArr[i];

    }
  }

}



//
// Get array of distances for sources
//
void distArrCalc( double *sourceDistArr ,  // Array to overwrite
                  int          *indexes ,  // Source locations
                  PixelMap     *distMap ,  // Map of distances
                  double          scale ,  // Mpc/rad conversion
                  int         N_sources ){ // Number of source

  for(int i = 0 ; i < N_sources ; ++i ){
    sourceDistArr[i] = (*distMap).getValue( indexes[i] ) * scale;
  }
}



//
//  Get random indexes to use as source positions, sources not allowed within some pixels
//
void getRandomSourcesIndexes( int     *indexes ,
                              userInfo       u ){

  //Will use 2d points, convert to one index later
  int xPos[ u.getNsrc() ], yPos[ u.getNsrc() ], counter;
  //Initialize so we can compare, not use same index twice or too close
  for(int i = 0; i < u.getNsrc(); ++i ){
       xPos[i] = -5;
       yPos[i] = -5;
    indexes[i] = -1;
  }
  bool accept;

  // Fill in each source
  // Assume point is good when pick point
  // Compare to all previous points
  // If any points too close, dont accept, repeat
  // Closeness can be specified by user

  for(int i = 0 ; i < u.getNsrc() ; ++i ){
    counter = 0;

    do{
      accept=true;

      // Generate random source pixels, within edge buffer
      xPos[i] =floor( rand()/(float)RAND_MAX * ( u.getNpixH() - 2*u.getEdgePix() )) + u.getEdgePix();
      yPos[i] =floor( rand()/(float)RAND_MAX * ( u.getNpixV() - 2*u.getEdgePix() )) + u.getEdgePix();

      // Compare against previous points
      for(int j=i-1;j>=0;--j){
        if(  sqrt( (xPos[i]-xPos[j])*(xPos[i]-xPos[j])
                 + (yPos[i]-yPos[j])*(yPos[i]-yPos[j]) ) < u.getMinNeighborDist() )
        accept=false;
      }
      // Abort if have too hard a time fitting sources into grid
      ++counter;
      if ( counter > 1e6 ){
        printf("Error: Unable to fit %4i sources in %4i x %4i pixelmap\n", u.getNsrc()  ,
                                                                           u.getNpixH() ,
                                                                           u.getNpixV() );
        printf(" Nearest neighbor tolerance: %5.2lf, Edgedist: %4i\n",     u.getMinNeighborDist(),
                                                                           u.getEdgePix());
        logMessage( std::string("Unable to place sources") );
        logMessage( std::string("           Nsrc  = ") + std::to_string( u.getNsrc()      ) +
                    std::string("           NpixH = ") + std::to_string( u.getNpixH()     ) +
                    std::string("           NpixV = ") + std::to_string( u.getNpixV()     ) +
                    std::string("  Min Seperation = ") + std::to_string( u.getMinNeighborDist() ) +
                    std::string("  Min Edge Dist  = ") + std::to_string( u.getEdgePix() ) );
        logMessage( std::string("Aborting."));
        exit(1);
      }
    }while(accept==false);

    indexes[i] = xPos[i] + yPos[i] * u.getNpixH();

  }
}




//
// Determines distance as function from center of pixelmap, in ang units
//
void distMapCalc( PixelMap  &distMap ,  // pixelmap to output
                  int      N_pixelsH ,  // Number of pixels on a side
                  int      N_pixelsV ,  // Number of pixels on a side
                  double     inpSize ,  // Angular size of x side
                  double   center[2] ){ // Center values

  double posArr[2]={0,0};
  double step     = inpSize/N_pixelsH;
//  #pragma omp parallel for private(posArr)
  for (int i=0;i<N_pixelsV;++i){
    posArr[1] = -i - 0.5 + N_pixelsV/2.0;
  for (int j=0;j<N_pixelsH;++j){
    posArr[0] =  j + 0.5 - N_pixelsH/2.0;
    distMap[j+i*N_pixelsH]= step * \
              sqrt( pow( posArr[0]-center[0] ,2) + pow( posArr[1]-center[1] ,2) );
  }
  }

}



//
//  Calculate the lensing quantities from the grid
//  Overwrites all the input PixelMaps
//
void calcLensMaps(  GridMap     &inpGrid ,  //GLAMER grid to calc values on
                    PixelMap   &kappaMap ,
                    PixelMap  &gamma1Map ,
                    PixelMap  &gamma2Map ,
                    PixelMap  &invMagMap ,
                    PixelMap   &g_tanMap ,
                    PixelMap   &g_secMap ,
                    PixelMap    &distMap ,
                    int       N_pixels_h ,  // Number of pixels on a side
                    int       N_pixels_v ,  // Number of pixels on a side
                    double      realSize ,  // Angular width in horizontal direction
                    double     center[2] ){ // Center location of halo




   kappaMap = inpGrid.writePixelMapUniform( center, N_pixels_h, N_pixels_v, KAPPA );  logMessage( std::string("Kappa  map populated") );
  gamma1Map = inpGrid.writePixelMapUniform( center, N_pixels_h, N_pixels_v, GAMMA1);  logMessage( std::string("Gamma1 map populated") );
  gamma2Map = inpGrid.writePixelMapUniform( center, N_pixels_h, N_pixels_v, GAMMA2);  logMessage( std::string("Gamma2 map populated") );
  invMagMap = inpGrid.writePixelMapUniform( center, N_pixels_h, N_pixels_v, INVMAG);  logMessage( std::string("invMag map populated") );


  distMapCalc( distMap, N_pixels_h, N_pixels_v, realSize, center);

  logMessage( std::string("Dist   map populated") );

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

    double a = 2*phi;

    g_tanMap[k] = (-gamma1Map[k]*cos(a)-gamma2Map[k]*sin(a)) / (1-kappaMap[k]);
    g_secMap[k] = (-gamma1Map[k]*sin(a)+gamma2Map[k]*cos(a)) / (1-kappaMap[k]);
  }
  }


  logMessage( std::string("gtan   map populated") );
  logMessage( std::string("gsec   map populated") );

/*
  printPixelMap(   distMap, N_pixels_h, N_pixels_v );
  printPixelMap(  kappaMap, N_pixels_h, N_pixels_v );
  printPixelMap( gamma1Map, N_pixels_h, N_pixels_v );
  printPixelMap( gamma2Map, N_pixels_h, N_pixels_v );
  printPixelMap(  g_tanMap, N_pixels_h, N_pixels_v );
  printPixelMap(  g_secMap, N_pixels_h, N_pixels_v );
*/
}



void printPixelMap( PixelMap   &inpMap ,
                    int       N_pixels ){
  printPixelMap( inpMap, N_pixels, N_pixels);
}

void printPixelMap( PixelMap   &inpMap   ,  //input pixel map
                    int       N_pixels_h ,
                    int       N_pixels_v ){ //Number of pixels on a side
  std::cout << std::endl;
  for( int i=0; i<N_pixels_v; i++){
  for( int j=0; j<N_pixels_h; j++){
    printf("%12.3e", inpMap[j+i*N_pixels_h]);
  }
    std::cout<<std::endl;
  }
}




// Bon-Muller transformation to provide gaussian distribution
double gaussErr( userInfo   u ,
                 int     Ngal ){

  float x1, x2, w;

  do {

    x1 = 2.0 * randVal( 0.0, 1.0 ) - 1.0;
    x2 = 2.0 * randVal( 0.0, 1.0 ) - 1.0;

    w  = x1 * x1 + x2 * x2;

  } while ( w >= 1.0 );

  w = std::sqrt( ( -2.0 * std::log( w ) ) / w );

  return x1 * w * u.getShapeNoise() / std::sqrt( Ngal );
}





void jacknife( densProfile  *profile ,
               int         N_samples ){

  // Variances

  double M_var(0);
  double C_var(0);
  double A_var(0);


  double Mbar_o(0);
  double Cbar_o(0);
  double Abar_o(0);


  double Mbar_i[ N_samples ];
  double Cbar_i[ N_samples ];
  double Abar_i[ N_samples ];


  // Generate xbar_i
  for ( int i = 0; i < N_samples; ++i ){
      Abar_i[ i ]  = 0;
      Mbar_i[ i ]  = 0;
      Cbar_i[ i ]  = 0;
  for ( int j = 0; j < N_samples; ++j ){

    if ( j != i ){
      if ( profile[j+1].getType() == 2 )
      Abar_i[ i ] +=             profile[j+1].getAlpha()  ;
      Mbar_i[ i ] += std::log10( profile[j+1].getM_enc() );
      Cbar_i[ i ] +=             profile[j+1].getC    ()  ;
    }

  }
      Abar_i[ i ]  = Abar_i[ i ] / ( N_samples - 1 );
      Cbar_i[ i ]  = Cbar_i[ i ] / ( N_samples - 1 );
      Mbar_i[ i ]  = Mbar_i[ i ] / ( N_samples - 1 );

  }


  // Generate xbar_o
  for ( int i = 0; i < N_samples; ++i ){
    Abar_o += Abar_i[ i ] / N_samples;
    Cbar_o += Cbar_i[ i ] / N_samples;
    Mbar_o += Mbar_i[ i ] / N_samples;
  }

  double coeff = ( N_samples - 1.0 ) / N_samples;

  // Generate variance
  for ( int i = 0; i < N_samples; ++i ){
    A_var += ( Abar_i[i] - Abar_o ) * ( Abar_i[i] - Abar_o ) ;
    C_var += ( Cbar_i[i] - Cbar_o ) * ( Cbar_i[i] - Cbar_o ) ;
    M_var += ( Mbar_i[i] - Mbar_o ) * ( Mbar_i[i] - Mbar_o ) ;
  }

  A_var = A_var * ( N_samples - 1.0 ) / N_samples ;
  C_var = C_var * ( N_samples - 1.0 ) / N_samples ;
  M_var = M_var * ( N_samples - 1.0 ) / N_samples ;


  // Biases


  double Aavg  = 0 ;
  double Mavg  = 0 ;
  double Cavg  = 0 ;

  for ( int i = 0; i < N_samples; ++i ){
    if ( profile[i+1].getType() == 2 )
    Aavg +=             profile[i+1].getAlpha()   / N_samples ;
    Mavg += std::log10( profile[i+1].getM_enc() ) / N_samples ;
    Cavg +=             profile[i+1].getC    ()   / N_samples ;
  }


  double theta_A = 0 ;
  double theta_M = 0 ;
  double theta_C = 0 ;

  for ( int i = 0; i < N_samples; ++i ){
    if ( profile[0].getType() == 2 )
    theta_A += pow( (             profile[i+1].getAlpha()   - Aavg ), 2 ) / N_samples ;
    theta_M += pow( ( std::log10( profile[i+1].getM_enc() ) - Mavg ), 2 ) / N_samples ;
    theta_C += pow( (             profile[i+1].getC    ()   - Cavg ), 2 ) / N_samples ;
  }

  // With i removed
  double thetaHat_A[N_samples] ;
  double thetaHat_M[N_samples] ;
  double thetaHat_C[N_samples] ;

  double thetaBar_A(0);
  double thetaBar_M(0);
  double thetaBar_C(0);


  for ( int i = 0; i < N_samples; ++i ){
    thetaHat_A[i] = 0;
    thetaHat_M[i] = 0;
    thetaHat_C[i] = 0;
/*
    double tempAavg=0;
    double tempMavg=0;
    double tempCavg=0;

  for ( int j = 0; j < N_samples; ++j ){
    if ( i!=j ){
      if ( profile[0].getType() == 2 )
      tempAavg +=             profile[j+1].getAlpha()   / (N_samples-1) ;
      tempMavg += std::log10( profile[j+1].getM_enc() ) / (N_samples-1) ;
      tempCavg +=             profile[j+1].getC    ()   / (N_samples-1) ;
    }
  }
//*/
  for ( int j = 0; j < N_samples; ++j ){
    if ( i!=j ){
      if ( profile[0].getType() == 2 )
      thetaHat_A[i] += pow( (             profile[j+1].getAlpha()   - Aavg ), 2 ) / (N_samples-1) ;
      thetaHat_M[i] += pow( ( std::log10( profile[j+1].getM_enc() ) - Mavg ), 2 ) / (N_samples-1) ;
      thetaHat_C[i] += pow( (             profile[j+1].getC    ()   - Cavg ), 2 ) / (N_samples-1) ;

    }
  }

    thetaBar_A += thetaHat_A[i] / N_samples;
    thetaBar_M += thetaHat_M[i] / N_samples;
    thetaBar_C += thetaHat_C[i] / N_samples;

  }


  double Abias ;
  if ( profile[0].getType() == 2 )
         Abias = ( N_samples - 1) * ( thetaBar_A - theta_A );//            profile[0].getAlpha()   ) ;
  double Mbias = ( N_samples - 1) * ( thetaBar_M - theta_M );//std::log10( profile[0].getM_enc() ) ) ;
  double Cbias = ( N_samples - 1) * ( thetaBar_C - theta_C );//            profile[0].getC    ()   ) ;


printf("Averages : %10.6f %10.6f %14.4e\n", Aavg  , Cavg  ,          Mavg    );
printf("Variance : %10.6f %10.6f %14.4e\n", A_var , C_var ,          M_var   );
printf("X_o      : %10.6f %10.6f %14.4e\n", Abar_o, Cbar_o, pow( 10, Mbar_o) );
printf("Biases   : %10.6f %10.6f %14.4e\n", Abias , Cbias ,          Mbias   );
printf("Theta    : %10.6f %10.6f %14.4e\n",theta_A,theta_C,         theta_M  );
printf("Theta_Bar: %10.6f %10.6f %14.4e\n",thetaBar_A,thetaBar_C,         thetaBar_M  );
printf("Corrected: %10.6f %10.6f %14.4e\n", Abias, N_samples*theta_C-(N_samples-1)*thetaBar_C ,
                                          pow( 10, N_samples*theta_M-(N_samples-1)*thetaBar_M ) );


}

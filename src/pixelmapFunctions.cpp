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
                         double    center[2] ){

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


  // Calculate weighted average of each bin
  for(int i = 0; i < u.getNbins(); ++i){

    if (N_countsArr[i]==0){ // If nothing in bin, mark as -1
      avgArr[i] = -1.0;
    }

    else{

      errArr[i] = sqrt( 1   /   errArr[i]             );
      avgArr[i] = avgArr[i] * ( errArr[i] * errArr[i] );

    }
  }

}



//
//  Radially average pixelmap values for given source positions
//
void radialDistAverage( double       *avgArr ,  // Array to overwrite
                        double    *distances ,  // Distance array of sources
                        userInfo           u ,
                        double     center[2] ){

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

    // Distance of pixels from center, converted to "bin" units
    iBin = round( distances[i] / u.getPhysFOV() * 2 * u.getNbins() );

         avgArr[iBin] += distances[i];
    N_countsArr[iBin] += 1;
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

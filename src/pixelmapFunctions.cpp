#include <slsimlib.h>
#include <image_processing.h>
#include <gridmap.h>
#include <grid_maintenance.h>
#include <astro_constants.h>
#include <lensing_classes.h>
#include <pixelmap_functions.h>






//
//
//
void distMapCalc( PixelMap  &distMap ,  // pixelmap to output
                  int      N_pixelsH ,  // Number of pixels on a side
                  int      N_pixelsV ,  // Number of pixels on a side
                  double     inpSize ,  // Angular size of x side
                  double   center[2] ){ // Center values

  double posArr[2]={0,0};
  double step     = inpSize/N_pixelsH;
  #pragma omp parallel for private(posArr)
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


   kappaMap = inpGrid.writePixelMapUniform( center, N_pixels_h, N_pixels_v, KAPPA);
  gamma1Map = inpGrid.writePixelMapUniform( center, N_pixels_h, N_pixels_v, GAMMA1);
  gamma2Map = inpGrid.writePixelMapUniform( center, N_pixels_h, N_pixels_v, GAMMA2);
  invMagMap = inpGrid.writePixelMapUniform( center, N_pixels_h, N_pixels_v, INVMAG);


  distMapCalc( distMap, N_pixels_h, N_pixels_v, realSize, center);

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
    //glamer produces negative value for g1, swap signs in eqtns
    double a = 2*phi;
    g_tanMap[k] = ( gamma1Map[k]*cos(a)-gamma2Map[k]*sin(a)) / (1-kappaMap[k]);//gamma1 swapped
    g_secMap[k] = (-gamma1Map[k]*sin(a)-gamma2Map[k]*cos(a)) / (1-kappaMap[k]);
  }
  }

/*
  printPixelMap(   distMap, N_pixels_h, N_pixels_v );

  printPixelMap(  kappaMap, N_pixels_h, N_pixels_v );

  printPixelMap( gamma1Map, N_pixels_h, N_pixels_v );

  printPixelMap( gamma2Map, N_pixels_h, N_pixels_v );
*/

}
/*
//double temp = densProfile.getRho_o() * densProfile.getR_s  ()   * exp( 2./alpha ) * pow( alpha/2., 1./alpha - 1.0 ) * sqrt( M_PI ) / Sc;
//double kappa    = foxH2012( x, alpha );
//double kappaAVG = temp * foxH2123( x, alpha );

//sqrt( g_tanMap[k]* g_tanMap[k]+ g_aziMap[k]* g_aziMap[k])
//sqrt(gamma1Map[k]*gamma1Map[k]+gamma2Map[k]*gamma2Map[k])
//     gamma2Map[k]
//     SDNFW( distMap[k], densProfile)/Sc
// SDNFWFull( distMap[k], densProfile.getR_s(), densProfile.getRho_o())/Sc

  for (int i=0;i<N_pixels;++i){
  for (int j=0;j<N_pixels;++j){
    int k = j+i*N_pixels;
    printf("%12.3e ",
g_tanMap[k]
 );
  }
    printf("\n");
  }
    printf("\n");
    printf("\n");


  for (int i=0;i<N_pixels;++i){
  for (int j=0;j<N_pixels;++j){
    int k = j+i*N_pixels;
    printf("%12.3e ",
g_aziMap[k]
 );
  }
    printf("\n");
  }
    printf("\n");
    printf("\n");


  for (int i=0;i<N_pixels;++i){
  for (int j=0;j<N_pixels;++j){
    int k = j+i*N_pixels;
    printf("%12.3e ",
g_tanMap[k]/
sqrt( g_tanMap[k]* g_tanMap[k]+ g_aziMap[k]* g_aziMap[k])
 );
  }
    printf("\n");
  }
    printf("\n");
    printf("\n");

//*/



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






/*
//Get array of distances for sources
void distArrCalc( double *sourceDistArr ,
                  int          *indexes ,
                  userInfo            u ,
                  double          scale ,
                  double      center[2] ){
  double posArr[2];
  for(int i=0;i<u.N_sources;++i){

    //0 is y, or row
    //1 is x, or columns
    //Pixel numbers, from top corner
    posArr[0] =  indexes[i]              % u.N_pixels;
    posArr[1] = (indexes[i] - posArr[0]) / u.N_pixels;
    //Pixel coordinates relative to center
    posArr[0] = ( posArr[0] - 0.5 - u.N_pixels/2.0 ) - center[1]+1;
    posArr[1] = (-posArr[1] + 0.5 + u.N_pixels/2.0 ) - center[0]-1;
    //Get position from pixel coordinates
    sourceDistArr[i] = scale*sqrt(posArr[0]*posArr[0]+posArr[1]*posArr[1]);
  }

}
*/
/*
  Radially average pixelmap values for given source positions
*//*
void radialSourceAverage( double  *avgArr  ,
                          double  *errArr  ,
                          int    *indexes  ,
                          PixelMap  inpMap ,
                          double   *errors ,
                          userInfo       u ,
                          double center[2] ){

  int    iBin;
  double dist;
  double posArr[2]={0,0};

  int N_countsArr[u.N_bins];
  //Make sure average is 0 to start
  for(int i=0;i<u.N_bins;++i){
         avgArr[i]=0;
         errArr[i]=0;
    N_countsArr[i]=0;
  }

  for(int i=0;i<u.N_sources;++i){

    //0 is y, or row
    //1 is x, or columns
    //Pixel numbers, from top corner
    posArr[0] =  indexes[i]              % u.N_pixels;
    posArr[1] = (indexes[i] - posArr[0]) / u.N_pixels;
    //Pixel coordinates relative to center
    posArr[0] = ( posArr[0] - 0.5 - u.N_pixels/2.0 ) - center[1]+1;
    posArr[1] = (-posArr[1] + 0.5 + u.N_pixels/2.0 ) - center[0]-1;

    //Get position from pixel coordinates
    dist = sqrt(posArr[0]*posArr[0]+posArr[1]*posArr[1]);
    iBin = round(std::min(std::max(
            dist/u.N_pixels*2.0 * u.N_bins ,0.0),1.0*u.N_bins));

    //To the bin, add the map value of this index, sum errors in quadrature
         avgArr[iBin] += inpMap[indexes[i]];
         errArr[iBin] += errors[i]*errors[i];
    N_countsArr[iBin] += 1;
  }

  for(int i=0;i<u.N_bins;++i){
    if (N_countsArr[i]==0){
      avgArr[i] = -1.0;
    }
    else{
      avgArr[i] =      avgArr[i]/N_countsArr[i] ;
      errArr[i] = sqrt(errArr[i]/N_countsArr[i]);
    }
  }
}

*/
/*
  Get random indexes to use as source positions, sources not allowed within some pixels
*//*
void getRandomSourcesIndexes( int     *indexes ,
                              userInfo       u ){

  //Will use 2d points, convert to one index later
  int xPos[u.N_sources], yPos[u.N_sources], counter;
  //Initialize so we can compare, not use same index twice or too close
  for(int i=0;i<u.N_sources;++i){
       xPos[i] = -5;
       yPos[i] = -5;
    indexes[i] = -1;
  }
  bool accept;
  //Fill in each source
  //Assume point is good when pick point
  //Compare to all previous points
  //If any points too close
  //Dont accept, repeat
  for(int i=0;i<u.N_sources;++i){
    counter=0;
    do{
      accept=true;

      xPos[i] =floor(rand()/(float)RAND_MAX * (u.N_pixels-2*u.N_edgepixels))
                +u.N_edgepixels;
      yPos[i] =floor(rand()/(float)RAND_MAX * (u.N_pixels-2*u.N_edgepixels))
                +u.N_edgepixels;
      //Compare against previous points
      for(int j=i-1;j>=0;--j){
        if(  sqrt( (xPos[i]-xPos[j])*(xPos[i]-xPos[j])
                 + (yPos[i]-yPos[j])*(yPos[i]-yPos[j]) ) < u.nearestSourceNeighbor )
        accept=false;
      }
      //Abort if have too hard a time fitting sources into grid
      ++counter;
      if (counter>1e6){
        printf("Error: Unable to fit %4i sources in %4i square pixelmap\n",u.N_sources,
                                                                           u.N_pixels);
        printf(" Nearest neighbor tolerance: %5.2lf, Edgedist: %4i\n",
                                              u.nearestSourceNeighbor,
                                              u.N_edgepixels);
        exit(0);
      }
    }while(accept==false);
    indexes[i] = xPos[i] + yPos[i] * u.N_pixels;
  }
}

*/







/*Does radial averaging from center of double **
//
//input:
//         double ** to average, arrays of length Nbin for average and distances
//         N_pixels, length of size of square double **
//          angRange, angular FOV of double **
//          N_bins, number of radial bins and size of average/dist array
//
//output:
//           average, distance array overwritten
//
*/
/*
void findRadialAverage(  PixelMap      &inpMap ,
                         int          N_pixels ,
                         int            N_bins ,
                         double   *avgInpArray ,
                         double  *distInpArray ,
                         double       realSize ,
                         double      center[2] ){

  int iBin;
  double averageArray[N_bins], NinsideArray[N_bins], distArray[N_bins];
  double dist;
  PosType posArr[2]={0,0};

  //Initialize the average, count array
  for ( int i=0; i<N_bins; ++i){
       distArray[i]=0;
    NinsideArray[i]=0;
    averageArray[i]=0;
  }

  //Loops over pixels, finds bin based on distance from center
  // then sum for average distance and count
  for (int i=0; i<N_pixels; ++i){
    posArr[1] = (- i - 0.5 + N_pixels/2.0) * realSize/N_pixels;
  for (int j=0; j<N_pixels; ++j){
    posArr[0] = (  j + 0.5 - N_pixels/2.0) * realSize/N_pixels;

    dist = sqrt( pow( posArr[0]-center[0] ,2) + pow( posArr[1]-center[1] ,2) );
    if ( dist<(realSize/2.0) ){
      iBin = round( dist/(realSize/2.0) * N_bins )-1;
         distArray[iBin] += dist;
      NinsideArray[iBin] += 1;
      averageArray[iBin] += inpMap[j+i*N_pixels];
    }
  }
  }
  //average the values
  //distInpArray[0] = 0.0;
  for ( int i=0; i<N_bins; ++i){
    distInpArray[i] =    distArray[i] / NinsideArray[i];
     avgInpArray[i] = averageArray[i] / NinsideArray[i];
    if ( (i>0)  &&   (distInpArray[i]!= distInpArray[i]) ) //Handle NaNs
    distInpArray[i] =(distInpArray[i+1]-distInpArray[i-1])/2.0;
  }
}
*/

/*Prints double ** in 2d
//
//input:
//     square pixel map, number of pixels on a side
//
*/
/*
*/

/*Calculate surface density in one of the square pixel maps
//
//input:
//     inpMap for surface density, convergence double **, number of pixels per side, critSD
//
//output:
//     inpMap/SD
//
*/
/*
void sigmaMapCalc( PixelMap  &inpMap ,
                   PixelMap    &kMap ,
                   int      N_pixels ,
                   double Sigma_crit ){
  #pragma omp parallel for
  for (int i=0; i<N_pixels; ++i){
  for (int j=0; j<N_pixels; ++j){
    inpMap[j+i*N_pixels] = kMap[j+i*N_pixels]*Sigma_crit;
  }
  }
}
*/

/*Calculates grav potential based on mass distribution in image, phi=-GM/r
//
//input:
//     square potential map, surface dens map, number of pixels on side, ang&real FOV
//
//output:
//     void, potential Map
//
*/
/*
void phiMapCalc( PixelMap   &phiMap ,
                 PixelMap &sigmaMap ,
                 int       N_pixels ,
                 double    realSize ){

  //Physical size of pixel
  double pixelSize   = ( realSize * realSize ) / ( N_pixels * N_pixels );
  double pixelLength =   realSize / N_pixels;
  double pixel1Loc[2],pixel2Loc[2];
  //phi(x) = -GM/r
  #pragma omp parallel for private(pixel1Loc,pixel2Loc)
  for(int  i=0; i<N_pixels;++i ){
    pixel1Loc[1] = -i  - 0.5 + N_pixels/2.0;
  for(int  j=0; j<N_pixels;++j ){
    pixel1Loc[0] =  j  + 0.5 - N_pixels/2.0;
    phiMap[j+i*N_pixels]=0;

  for(int ii=0;ii<N_pixels;++ii){
    pixel2Loc[1] = -ii - 0.5 + N_pixels/2.0;
  for(int jj=0;jj<N_pixels;++jj){
    pixel2Loc[0] =  jj + 0.5 - N_pixels/2.0;

    if( !( (i==ii)&&(j==jj) ) ) {
      //Partial component of phi
      phiMap[j+i*N_pixels] += sigmaMap[jj+ii*N_pixels] / (sqrt(
        pow(pixel1Loc[0]-pixel2Loc[0],2)+pow(pixel1Loc[1]-pixel2Loc[1],2) ));
    }
    else
      phiMap[j+i*N_pixels] += 0.0;
  }
  }
    phiMap[j+i*N_pixels] = - phiMap[j+i*N_pixels] * pixelSize * G_COSMO;
    //rest of phi, remove rad^2 * G

  }
  }

}
*/

/*Generates mag map using lensing parameters
//
//input:
//          magMap to change, kappa & gamma double **
//
//output:
//          magMap
//
*/
/*
void magMapCalc( PixelMap  &magMap ,
                 PixelMap    &kMap ,
                 PixelMap    &gMap ,
                 int      N_pixels ){
  #pragma omp parallel for
  for (int i = 0; i<N_pixels; ++i){
  for (int j = 0; j<N_pixels; ++j){
    magMap[i*N_pixels+j] = 1.0/( ( 1-kMap[i*N_pixels+j] )*( 1-kMap[i*N_pixels+j] ) \
                                 - ( gMap[i*N_pixels+j]  *    gMap[i*N_pixels+j]));
    // 1/[ (1-k)^2 - g^2 ]
  }
  }
}
*/

/*Calculates the mass ENCLOSED TO CENTER BY A PIXEL
//
//input:
//     double **s to write mass, and containing potential, dist, grid size on a side
//
//output:
//     void, massMap
//
*/
/*
void massMapCalc( PixelMap   &massMap ,
                  PixelMap    &phiMap ,
                  PixelMap   &distMap ,
                  int        N_pixels ){
  #pragma omp for
  for(int i=0;i<N_pixels;++i){
  for(int j=0;j<N_pixels;++j){
    massMap[j+i*N_pixels]=-phiMap[j+i*N_pixels]*distMap[j+i*N_pixels]/G_COSMO;
  }
  }
}
*/


/*
  Calculates some values onto pixel maps, uses convergence at all pixels
  Due to nature, only can be done with simulated data/models, otherwise
    impossible to detect in observations
*/
/*
void calcMapsFromKappa(
  PixelMap  &kappaMap,  //Convergence map to pull from
  PixelMap    &phiMap,  //Potential in pixel
  PixelMap   &massMap,  //Mass estimate enclosed in pixel
  PixelMap   &distMap,  //Distance from center, in pixel value
  PixelMap  &sigmaMap,  //Surface density of pixel
  double    center[2],  //Center of map, in Mpc
  userInfo    inpInfo,  //User input information
  double    sigmaCrit,  //Sc of maps from source
  double     realSize){ //Physical size of map, in Mpc

  std::cout << "Constructing maps from kappa:" << std::endl;
   sigmaMapCalc( sigmaMap, kappaMap,          inpInfo.N_pixels, sigmaCrit        );
  std::cout << " -sigma map constructed" << std::endl;
    distMapCalc( distMap,                     inpInfo.N_pixels, realSize,  center);
  std::cout << " -dist  map constructed" << std::endl;
    phiMapCalc(  phiMap,   sigmaMap,          inpInfo.N_pixels, realSize );
  std::cout << " -phi   map constructed" << std::endl;
   massMapCalc( massMap,     phiMap, distMap, inpInfo.N_pixels );
  std::cout << " -mass  map constructed" << std::endl;

}
*/

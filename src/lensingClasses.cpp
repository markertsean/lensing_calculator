#include "lensing_classes.h"



haloInfo::haloInfo(){

       z=-1.0;        // Redshift
       m=-1.0;        // Mass
       c=-1.0;        // Concentration
    rmax=-1.0;        // Rvir, may be R200
      ba=-1.0;        // b/a axis ratio
      ca=-1.0;        // c/a axis ratio
     phi=-1.0;        // Orientation (out of the page?)
   theta=-1.0;        // Orientation (on the page?)
   alpha=-1.0;        // Orientation xy plane
   gamma=-1.0;        // Orientation along z
      id=-1  ;        // Halo ID number
}



userInfo::userInfo(){
        angFOV = -1.  ; //  Angular size of image
       physFOV = -1.  ; // Physical size of image
         R_max = -1.  ; // ?
   sourceR     =  1.0 ; // Radius of source in pixels
   sourceZ     = -1.  ; // Redshift of sources
   sourceDens  = 10.0 ; // Surface density of background sources, gal/arcmin^2
   shapeNoise  =  0.3 ; // Intrinsic shape noise in the sources
   N_pixels    = -1   ; // Number of pixels on the grid
   N_pixels_h  = -1   ; // Number of pixels on x-axis
   N_pixels_v  = -1   ; // Number of pixels on y-axis
   N_bins      = -1   ; // Number of bins for radial averaging
N_bins_R2D  = 10  ; // Number of bins for radial averaging
N_bins_A2D  = 10   ; // Number of bins for radial averaging
   N_sources   = -1   ; // Number of sources to generate
   N_particles = -1   ; // Number of particles in simulation
   num_threads =  1   ; // Number of threads for parallel processing
  N_edgepixels =  3   ; // Number of pixels to leave on an edge

  nearestSourceNeighbor = 1.5;   // Minimum distance between sources
  readFile     = " ";            // Fits file to read
  catType      = " ";            // Catalog file to read
  cosmo        = " ";            // Type of cosmology

  fox2012F = "src/foxH2012.dat"; // FoxH files to read
  fox2123F = "src/foxH2123.dat";

  outputPath = "data/";

      cMin =  2.5;     // Range of concentration values to fit
      cMax =  7.5;
      rMin =  0.3;     // Range of R_max values to fit, Mpc
      rMax =  2.0;
      mMin = 12.5;     // Range of mass values to fit
      mMax = 16.5;
  alphaMin =  9e-2;  // Range of alpha values to fit
  alphaMax =  0.68;

/*
  maxFitAttempts = 1e4   ;
   N_chromosomes = 1e4   ;
      N_chiTrack = 1e2   ;
      consistent = 2e2   ;
       tolerance = 1e-3  ;
       mutChance = 1e-2  ;
      avgTestVal = 0.5   ;
*/



  maxFitAttempts = 1e3   ; // Maximum number of steps to roll ball, or times to reproduce
   N_chromosomes = 1e3   ; // Number of chromosomes or balls
      consistent = 2e1   ; // Number of steps to converge before accepting
       tolerance = 1e-5  ; // If difference between steps less than this, converged


      N_chiTrack = 1e1   ; // Number of chromosomes to check for convergence
       mutChance = 1e-2  ; // Likelihood of converging
      avgTestVal = 0.7   ; // Criteria for reproducing

}



densProfile::densProfile(){
  concentration = -1.0;
  alpha         = -1.0;
  r_max         = -1.0;
  M_enc         = -1.0;
  type          =    1;
}


// If passed an argument, Einasto profile
// Alternatively can change type later
densProfile::densProfile( double inpA ){
  concentration = -1.0;
  alpha         = inpA;
  r_max         = -1.0;
  M_enc         = -1.0;
  type          =    2;
}



double densProfile::getRho_o() const {

  if ( M_enc > 0 && r_max > 0 && concentration > 0){ // All needed parameters def

    if ( type == 1 ){ // NFW
                                         return  M_enc / (4. * M_PI * r_max*r_max*r_max )  *
                                                            concentration   * concentration       * concentration /
                                                 ( log( 1 + concentration ) - concentration / ( 1 + concentration) );
    } else {          // Einasto

      double rs = r_max / concentration;

                                         return  M_enc * alpha / (   4. * M_PI    *    rs * rs * rs   *
                                                    exp(             2. / alpha ) *
                                                    pow( alpha / 2., 3. / alpha ) *
                                                    tgamma(          3. / alpha ) );
    }
  }

  return -1;
}








// Log file stuff



// Initializes log file, making directory and file name
void initLogFile(){

  // Sets up the times
  int execution_start = clock();
  time_t      nowTime =      time(        0 );
  tm       *startTime = localtime( &nowTime );

  // Logfile name & directory
  char logFileNameC[100];
  struct stat sb;
  char str[] = "mkdir logfiles";

  // Create logfiles directory, if one does not exist
  if ( stat( "logfiles/", &sb) != 0 ){
    system( str );
  }

  // Log file name is mostly the date
  sprintf( logFileNameC, "logfiles/lensCalc.%4i.%02i.%02i.%02i.%02i.%02i.log",
    (*startTime).tm_year+1900,
    (*startTime).tm_mon ,
    (*startTime).tm_mday,
    (*startTime).tm_hour,
    (*startTime).tm_min ,
    (*startTime).tm_sec );

  logFileName= std::string(logFileNameC) ;

  // First line of the log file is the time
  logMessage( (std::string(                "Code initialized at ")+
               std::to_string( (long long) (*startTime).tm_hour  )+
               std::string(                ":"                   )+
               std::to_string( (long long) (*startTime).tm_min   )+
               std::string(                ":"                   )+
               std::to_string( (long long) (*startTime).tm_sec   )));


}


// Generates log files, printing text to file
void logMessage( const std::string &text ){

    std::ofstream log_file(  logFileName, std::ios_base::out | std::ios_base::app );
    log_file << text << std::endl;

}


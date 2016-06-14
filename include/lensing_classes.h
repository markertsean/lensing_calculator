#ifndef MY_LENSING_CLASSES
#define MY_LENSING_CLASSES

#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cosmo.h>
#include <math.h>
#include <sys/stat.h>
#include <astro_constants.h>




// Log file name
extern std::string logFileName;


// Generates log files, printing text to file
inline void logMessage( const std::string &text ){

    std::ofstream log_file(  logFileName, std::ios_base::out | std::ios_base::app );
    log_file << text << std::endl;

}

// Initializes log file, making directory and file name
inline void initLogFile(){

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



//Holds info on lens, source, any halo for easy access
class haloInfo{
  public:

    inline haloInfo();

    void setZ     ( double inpZ ){ z     = inpZ; }
    void setM     ( double inpM ){ m     = inpM; }
    void setC     ( double inpC ){ c     = inpC; }
    void setRmax  ( double inpR ){ rmax  = inpR; }
    void setBA    ( double inpF ){ ba    = inpF; }
    void setCA    ( double inpF ){ ca    = inpF; }
    void setPhi   ( double inpF ){ phi   = inpF; }
    void setTheta ( double inpF ){ theta = inpF; }
    void setID    ( long   inpI ){ id    = inpI; }

    double      getZ     () const { return z     ; }
    double      getM     () const { return m     ; }
    double      getC     () const { return c     ; }
    double      getRmax  () const { return rmax  ; }
    double      getBA    () const { return ba    ; }
    double      getCA    () const { return ca    ; }
    double      getPhi   () const { return phi   ; }
    double      getTheta () const { return theta ; }
    long        getID    () const { return id    ; }


  private:
    double       z;             // redshift
    double       m;             // mass Msun
    double       c;             // concentration
    double    rmax;             // max radius Mpc
    double      ba;             // b/a axis ratio
    double      ca;             // c/a axis ratio
    double     phi;             // y/z orientation angle
    double   theta;             // x/z orientation angle
    long        id;             // id number of halo


    void   unDefVar( std::string inpS ){
      std::cerr<<"WARNING: variable "<<inpS<<" undefined in <haloInfo>"<<std::endl;
    }
};

haloInfo::haloInfo(){

       z=-1.0;
       m=-1.0;
       c=-1.0;
    rmax=-1.0;
      ba=-1.0;
      ca=-1.0;
     phi=-1.0;
   theta=-1.0;
      id=-1  ;
}


// Holds info on density profile, concentration, mass, shape param, etc
class densProfile{
  public:

    inline densProfile();
    inline densProfile( double inpA );

    // Modifiers, as parameters are modified need to adjust dependant values
    void setAlpha( double inpA ){         alpha = inpA;                         modRho_o(); }
    void setR_max( double inpR ){         r_max = inpR;   modR_s();             modRho_o(); }
    void setC    ( double inpC ){ concentration = inpC;   modR_s();             modRho_o(); }
    void setC    ( double inpC, int i ){ concentration = inpC;printf("Rs\n");   modR_s();printf("Rho\n");             modRho_o(); }
    void setR_s  ( double inpR ){           r_s = inpR;               modC();   modRho_o(); }
    void setRho_o( double inpR ){         rho_o = inpR;                         modM_enc(); }
    void setM_enc( double inpM ){         M_enc = inpM;                         modRho_o(); }
    void resetM_enc(           ){                                               modM_enc(); }

// Rs    -> C, Rm
// C     -> Rs, Rm
// Rho_o -> M, Rm, C, alpha
// M     -> Rho_o, Rm, C

    // If when getting an undefined variable, need to spit out a warning
    double getC    () const {  if ( concentration==-1.0 ) unDefVar("\"concentration\"");   return concentration;  }
    double getRho_o() const {  if (         rho_o==-1.0 ) unDefVar("\"rho_o\""        );   return         rho_o;  }
    double getR_max() const {  if (         r_max==-1.0 ) unDefVar("\"R_max\""        );   return         r_max;  }
    double getR_s  () const {  if (           r_s==-1.0 ) unDefVar("\"R_scale\""      );   return           r_s;  }
    double getM_enc() const {  if (         M_enc==-1.0 ) unDefVar("\"M_enc\""        );   return         M_enc;  }
    double getAlpha() const {  if (         alpha==-1.0 ) unDefVar("\"alpha\""        );   return         alpha;  }
    double getType () const {                                                              return          type;  }


  private:
    double concentration;
    double alpha        ;
    double rho_o        ;
    double r_s          ;
    double r_max        ;
    double M_enc        ;
    int    type         ; //1 NFW, 2 Einasto

    void   unDefVar( std::string inpS ) const {
      std::cerr<<"WARNING: variable "<<inpS<<" undefined in <lensProfile>"<<std::endl;
    }

    void modRho_o(){
      if ( M_enc > 0 && r_max > 0 && concentration > 0){ //All needed parameters def
        if ( type == 1 ){
          rho_o   =  M_enc / (4. * M_PI * r_max*r_max*r_max )  *
                     concentration*concentration*concentration /
                    ( log( 1 + concentration ) - concentration / ( 1 + concentration) );
        }
        else {
          rho_o   =  M_enc / ( 2. * M_PI  *            r_s * r_s * r_s * exp( 4./alpha ) *
                              pow( alpha / 2., 3. / alpha - 1 ) *    tgamma( 3./alpha ) );
        }
      }
    }

    void modM_enc(){

      if ( rho_o > 0 && r_max > 0 && concentration > 0){ // All needed parameters def
        if ( type == 1 ){
          M_enc   =  rho_o * (4. * M_PI * r_max*r_max*r_max )  /
                     concentration*concentration*concentration *
                    ( log( 1 + concentration ) - concentration / ( 1 + concentration) );
        }
        else {
          M_enc   = 2. * M_PI  *  rho_o *   r_s * r_s * r_s * exp( 4./alpha ) *
                    pow( alpha / 2., 3. / alpha - 1 ) *    tgamma( 3./alpha );
        }
      }

    }

    void modR_s(){  if (   concentration > 0 && r_max > 0 )           r_s = r_max / concentration;  }
    void modC  (){  if (             r_s > 0 && r_max > 0 ) concentration = r_max / r_s          ;  }

};

densProfile::densProfile(){
  concentration = -1.0;
  alpha         = -1.0;
  rho_o         = -1.0;
  r_s           = -1.0;
  r_max         = -1.0;
  M_enc         = -1.0;
  type          =    1;
}


densProfile::densProfile( double inpA ){
  concentration = -1.0;
  alpha         = inpA;
  rho_o         = -1.0;
  r_s           = -1.0;
  r_max         = -1.0;
  M_enc         = -1.0;
  type          =    2;
}

// Holds info on user input values
class userInfo{
  public:

    inline userInfo();

    void setFOV             ( double inpF ) {   angFOV       = inpF ; }
    void setNpix            ( int    inpI ) {   N_pixels     = inpI ; }
    void setNpart           ( int    inpI ) {   N_particles  = inpI ; }

    void setSourceZ         ( double inpF ) {        sourceZ = inpF ; }
    void setChiMin          ( double inpF ) {           cMin = inpF ; }
    void setChiMax          ( double inpF ) {           cMax = inpF ; }
    void setMaxFitNum       ( int    inpI ) { maxFitAttempts = inpI ; }
    void setNConsistent     ( int    inpI ) {  consistent    = inpI ; }
    void setTolerance       ( double inpF ) {      tolerance = inpF ; }
    void setMutChance       ( double inpF ) {      mutChance = inpF ; }
    void setTestVal         ( double inpF ) {     avgTestVal = inpF ; }
    void setNtrack          ( int    inpI ) {  N_chiTrack    = inpI ; }
    void setMassMin         ( double inpF ) {           mMin = inpF ; }
    void setMassMax         ( double inpF ) {           mMax = inpF ; }
    void setNchrome         ( int    inpI ) {  N_chromosomes = inpI ; }
    void setConMin          ( double inpF ) {           cMin = inpF ; }
    void setConMax          ( double inpF ) {           cMax = inpF ; }
    void setAlphaMin        ( double inpF ) {       alphaMin = inpF ; }
    void setAlphaMax        ( double inpF ) {       alphaMax = inpF ; }
    void setEdgePix         ( int    inpI ) {   N_edgepixels = inpI ; }
    void setMinNeighborDist ( int    inpI ) { nearestSourceNeighbor = inpI; }
    void setNbins           ( int    inpI ) {   N_bins       = inpI ; }
    void setNsrc            ( int    inpI ) {   N_sources    = inpI ; }
    void setNthreads        ( int    inpI ) {   num_threads  = inpI ; }
    void setNpixH           ( int    inpI ) {   N_pixels_h   = inpI ; }
    void setNpixV           ( int    inpI ) {   N_pixels_v   = inpI ; }
    void setNgridPoints     ( int    inpI ) {   N_gridPoints = inpI ; }
    void setIntegLength     ( double inpF ) {   integLength  = inpF ; }
    void setSourceRadius    ( double inpF ) {   sourceR      = inpF ; }
    void setPhysFOV         ( double inpF ) {   physFOV      = inpF ; }
    void setAngFOV          ( double inpF ) {   angFOV       = inpF ; }
    void setCatType         ( std::string inpS ) {   catType = inpS ; }
    void setCosmology       ( std::string inpS ) {     cosmo = inpS ; }

    double getIntegLength     () const { return  integLength   ; }
    int    getNpixH           () const { return  N_pixels_h    ; }
    int    getNpixV           () const { return  N_pixels_v    ; }
    double getPhysFOV         () const { return  physFOV       ; }
    double getAngFOV          () const { return  angFOV        ; }
    double getSourceRadius    () const { return  sourceR       ; }
    int    getNpix            () const { return   N_pixels     ; }
    int    getNsrc            () const { return   N_sources    ; }
    int    getNbins           () const { return   N_bins       ; }
    int    getNthreads        () const { return    num_threads ; }
    int    getNgridPoints     () const { return   N_gridPoints ; }
    double getMinNeighborDist () const { return nearestSourceNeighbor ; }
    int    getEdgePix         () const { return   N_edgepixels ; }
    std::string getCatType    () const { return  catType       ; }
    std::string getCosmology  () const { return  cosmo         ; }
    double getAlphaMin        () const { return       alphaMin ; }
    double getAlphaMax        () const { return       alphaMax ; }
    double getConMin          () const { return           cMin ; }
    double getConMax          () const { return           cMax ; }
    int    getNchrome         () const { return  N_chromosomes ; }
    double getMassMin         () const { return           mMin ; }
    double getMassMax         () const { return           mMax ; }
    int    getNtrack          () const { return  N_chiTrack    ; }
    double getTestVal         () const { return     avgTestVal ; }
    double getTolerance       () const { return      tolerance ; }
    double getMutChance       () const { return      mutChance ; }
    int    getMaxFitNum       () const { return maxFitAttempts ; }
    int    getNConsistent     () const { return  consistent    ; }
    double getChiMin          () const { return           cMin ; }
    double getChiMax          () const { return           cMax ; }
    double getSourceZ         () const { return       sourceZ  ; }




    double getFOV             () { return   angFOV       ; }
    int    getNpart           () { return   N_particles  ; }


  private:

    // Stuff to read in
    double      angFOV;  // Angular size of image
    double     physFOV;  // Physical size of image
    double integLength;  // Integration length of current image
    double       R_max;  // ?

    double    sourceR ;  // Radius of sources in pixels
    double    sourceZ ;

    int    N_pixels   ;  // Number of pixels on grid
    int    N_pixels_h ;  // Number of pixels on x-axis
    int    N_pixels_v ;  // Number of pixels on y-axis

    int    N_bins     ;  // Number of bins for radial averaging
    int    N_sources  ;  // Number of sources to generate
    int    N_particles;  // Number of particles in simulation
    int    N_gridPoints; // Number of grid points to interpolate mass to
    int    num_threads;  // Number of threads for parallel processing
    int   N_edgepixels;  // Number of pixels to leave on an edge

    double nearestSourceNeighbor; // Gap in distance between sources

  std::string readFile;
  std::string catType ;
  std::string cosmo   ;

  // Chi2 & genetic algorithm fitting values
  double     cMin;
  double     cMax;
  double     mMin;
  double     mMax;
  double alphaMin;
  double alphaMax;

  int maxFitAttempts; // Max attempts at fitting before abort
  int  N_chromosomes; // Number of chromosomes in population
  int     N_chiTrack; // Number of chi's to track for convergence
  int     consistent; // Number of times need avg below tolerance
  double   tolerance; // Average residual must be below tolerance
  double   mutChance; // Likelihood of mutation
  double  avgTestVal; // chiAvg*this is random range
};

userInfo::userInfo(){
        angFOV = -1.;
       physFOV = -1.;
         R_max = -1.;
   sourceR     = 1.0;
   sourceZ     = -1.;
   N_pixels    = -1;
   N_pixels_h  = -1;
   N_pixels_v  = -1;
   N_bins      = -1;
   N_sources   = -1;
   N_particles = -1;
   num_threads =  1;
  N_edgepixels =  3;

  nearestSourceNeighbor = 1.5;
  readFile = " ";
  catType  = " ";
  cosmo    = " ";

      cMin =  2.0;
      cMax =  9.0;
      mMin = 12.0;
      mMax = 16.0;
  alphaMin =  0.1;
  alphaMax =  1.0;

//5.06
//1.00e14
//10%
/*
  maxFitAttempts = 1e4   ;
   N_chromosomes = 1e4   ;
      N_chiTrack = 1e2   ;
      consistent = 2e2   ;
       tolerance = 1e-3  ;
       mutChance = 1e-2  ;
      avgTestVal = 0.5   ;
*/

/*
  maxFitAttempts = 1e3   ;
   N_chromosomes = 1e1   ;
      N_chiTrack = 4e0   ;
      consistent = 1e1   ;
       tolerance = 1e-1  ;
       mutChance = 1e-2  ;
      avgTestVal = 0.6   ;
*/

  maxFitAttempts = 1e4   ;
   N_chromosomes = 1e1   ;
      N_chiTrack = 2e1   ;
      consistent = 1e1   ;
       tolerance = 1e-3  ;
       mutChance = 1e-2  ;
      avgTestVal = 0.6   ;

}

#endif

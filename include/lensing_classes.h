#ifndef MY_LENSING_CLASSES
#define MY_LENSING_CLASSES

#include <iostream>
#include <cosmo.h>
#include <math.h>
#include <astro_constants.h>

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

    double      getZ     () { return z     ; }
    double      getM     () { return m     ; }
    double      getC     () { return c     ; }
    double      getRmax  () { return rmax  ; }
    double      getBA    () { return ba    ; }
    double      getCA    () { return ca    ; }
    double      getPhi   () { return phi   ; }
    double      getTheta () { return theta ; }
    long        getID    () { return id    ; }


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
    void setAlpha( double inpA ){         alpha = inpA;   modRho_o(); }
    void setR_max( double inpR ){         r_max = inpR;   modR_s  ();   modC    ();   modRho_o(); }
    void setC    ( double inpC ){ concentration = inpC;   modR_s  ();   modRho_o(); }
    void setR_s  ( double inpR ){           r_s = inpR;   modC    ();   modRho_o(); }
    void setRho_o( double inpR ){         rho_o = inpR;   modM_enc(); }
    void setM_enc( double inpM ){         M_enc = inpM;   modRho_o(); }
    void resetM_enc(           ){                         modM_enc(); }

    // If when getting an undefined variable, need to spit out a warning
    double getC    () {  if ( concentration==-1.0 ) unDefVar("\"concentration\"");   return concentration;  }
    double getRho_o() {  if (         rho_o==-1.0 ) unDefVar("\"rho_o\""        );   return         rho_o;  }
    double getR_max() {  if (         r_max==-1.0 ) unDefVar("\"R_max\""        );   return         r_max;  }
    double getR_s  () {  if (           r_s==-1.0 ) unDefVar("\"R_scale\""      );   return           r_s;  }
    double getM_enc() {  if (         M_enc==-1.0 ) unDefVar("\"M_enc\""        );   return         M_enc;  }
    double getAlpha() {  if (         alpha==-1.0 ) unDefVar("\"alpha\""        );   return         alpha;  }
    double getType () {                                                              return          type;  }


  private:
    double concentration;
    double alpha        ;
    double rho_o        ;
    double r_s          ;
    double r_max        ;
    double M_enc        ;
    int    type         ; //1 NFW, 2 Einasto

    void   unDefVar( std::string inpS ){
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


    void setMinNeighborDist ( int    inpI ) { nearestSourceNeighbor = inpI; }
    void setFOV             ( double inpF ) {   angFOV       = inpF ; }
    void setRmax            ( double inpF ) {   R_max        = inpF ; }
    void setNpix            ( int    inpI ) {   N_pixels     = inpI ; }
    void setNsrc            ( int    inpI ) {   N_sources    = inpI ; }
    void setNbins           ( int    inpI ) {   N_bins       = inpI ; }
    void setNpart           ( int    inpI ) {   N_particles  = inpI ; }
    void setEdgePix         ( int    inpI ) {   N_edgepixels = inpI ; }
    void setNthreads        ( int    inpI ) {   num_threads  = inpI ; }
    void setChiMin          ( double inpF ) {           cMin = inpF ; }
    void setChiMax          ( double inpF ) {           cMax = inpF ; }
    void setMassMin         ( double inpF ) {           mMin = inpF ; }
    void setMassMax         ( double inpF ) {           mMax = inpF ; }
    void setAlphaMin        ( double inpF ) {       alphaMin = inpF ; }
    void setAlphaMax        ( double inpF ) {       alphaMax = inpF ; }
    void setTolerance       ( double inpF ) {      tolerance = inpF ; }
    void setMutChance       ( double inpF ) {      mutChance = inpF ; }
    void setTestVal         ( double inpF ) {     avgTestVal = inpF ; }
    void setMaxFitNum       ( int    inpI ) { maxFitAttempts = inpI ; }
    void setNchrome         ( int    inpI ) {  N_chromosomes = inpI ; }
    void setNtrack          ( int    inpI ) {  N_chiTrack    = inpI ; }
    void setNConsistent     ( int    inpI ) {  consistent    = inpI ; }

    int    getMinNeighborDist () { return nearestSourceNeighbor ; }
    double getFOV             () { return   angFOV       ; }
    double getRmax            () { return   R_max        ; }
    int    getNpix            () { return   N_pixels     ; }
    int    getNsrc            () { return   N_sources    ; }
    int    getNbins           () { return   N_bins       ; }
    int    getNpart           () { return   N_particles  ; }
    int    getEdgePix         () { return   N_edgepixels ; }
    int    getNthreads        () { return    num_threads ; }
    double getChiMin          () { return           cMin ; }
    double getChiMax          () { return           cMax ; }
    double getMassMin         () { return           mMin ; }
    double getMassMax         () { return           mMax ; }
    double getAlphaMin        () { return       alphaMin ; }
    double getAlphaMax        () { return       alphaMax ; }
    double getTolerance       () { return      tolerance ; }
    double getMutChance       () { return      mutChance ; }
    double getTestVal         () { return     avgTestVal ; }
    int    getMaxFitNum       () { return maxFitAttempts ; }
    int    getNchrome         () { return  N_chromosomes ; }
    int    getNtrack          () { return  N_chiTrack    ; }
    int    getNConsistent     () { return  consistent    ; }


  private:

    // Stuff to read in
    double      angFOV;  // Angular size of image
    double       R_max;
    int    N_pixels   ;  // Number of pixels on grid
    int    N_bins     ;  // Number of bins for radial averaging
    int    N_sources  ;  // Number of sources to generate
    int    N_particles;  // Number of particles in simulation
    int    num_threads;  // Number of threads for parallel processing
    int   N_edgepixels;  // Gap in distance between sources
    double nearestSourceNeighbor; // Number of pixels to leave on an edge

  std::string readFile;

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
         R_max = -1.;
   N_pixels    = -1;
   N_bins      = -1;
   N_sources   = -1;
   N_particles = -1;
   num_threads =  1;
  N_edgepixels =  3;

  nearestSourceNeighbor = 1.5;
  readFile = " ";

      cMin =  2.0;
      cMax =  8.0;
      mMin = 12.0;
      mMax = 17.0;
  alphaMin =  0.1;
  alphaMax =  0.4;

  maxFitAttempts =  1e4    ;
   N_chromosomes = 1000    ;
      N_chiTrack =  100    ;
      consistent =   50    ;
       tolerance =    0.01 ;
       mutChance =    0.01 ;
      avgTestVal =    1.3  ;

}

#endif

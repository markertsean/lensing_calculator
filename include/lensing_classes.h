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
void logMessage( const std::string &text );

// Initializes log file, making directory and file name
void initLogFile();



// Holds info on lens, source, any halo for easy access
class haloInfo{
  public:

//    inline haloInfo();
    haloInfo();

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



// Holds info on user input values
class userInfo{
  public:

    userInfo();

    void setFOV             ( double inpF ) {   angFOV       = inpF ; }
    void setNpix            ( int    inpI ) {   N_pixels     = inpI ; }
    void setNpart           ( int    inpI ) {   N_particles  = inpI ; }

    void setSourceZ         ( double inpF ) {        sourceZ = inpF ; }
    void setSourceDens      ( double inpF ) {     sourceDens = inpF ; }
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
    void setRMinFit         ( double inpF ) {           rMin = inpF ; }
    void setRMaxFit         ( double inpF ) {           rMax = inpF ; }
    void setNchrome         ( int    inpI ) {  N_chromosomes = inpI ; }
    void setShapeNoise      ( double inpF ) {     shapeNoise = inpF ; }
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
    void setFoxH2012F       ( std::string inpS ) {  fox2012F = inpS ; }
    void setFoxH2123F       ( std::string inpS ) {  fox2123F = inpS ; }

    double getIntegLength     () const { return  integLength   ; }
    int    getNpixH           () const { return  N_pixels_h    ; }
    int    getNpixV           () const { return  N_pixels_v    ; }
    double getPhysFOV         () const { return  physFOV       ; }
    double getAngFOV          () const { return  angFOV        ; }
    double getSourceRadius    () const { return  sourceR       ; }
    double getSourceDensity   () const { return  sourceDens    ; }
    double getShapeNoise      () const { return  shapeNoise    ; }
    int    getNpix            () const { return   N_pixels     ; }
    int    getNbins           () const { return   N_bins       ; }
    int    getNthreads        () const { return    num_threads ; }
    int    getNgridPoints     () const { return   N_gridPoints ; }
    double getMinNeighborDist () const { return nearestSourceNeighbor ; }
    int    getEdgePix         () const { return   N_edgepixels ; }
    std::string getCatType    () const { return  catType       ; }
    std::string getCosmology  () const { return  cosmo         ; }
    std::string getFoxH2012F  () const { return  fox2012F      ; }
    std::string getFoxH2123F  () const { return  fox2123F      ; }
    double getAlphaMin        () const { return       alphaMin ; }
    double getAlphaMax        () const { return       alphaMax ; }
    double getConMin          () const { return           cMin ; }
    double getConMax          () const { return           cMax ; }
    double getRMinFit         () const { return           rMin ; }
    double getRMaxFit         () const { return           rMax ; }
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

    int    getNsrc            () const { if (    N_sources != -1 ) // If users provided a number, use that
                                          return N_sources     ;
                                          return (int) (           // Otherwise use source density and size
                                                 sourceDens *      //   of our FOV
                                                 physFOV    *
                                                 physFOV    )  ; }




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

    double    sourceDens ; // Source Number Density, in gal/arcmin^2
    double shapeNoise ;

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

  std::string fox2012F;
  std::string fox2123F;
  std::string readFile;
  std::string catType ;
  std::string cosmo   ;

  // Chi2 & genetic algorithm fitting values
  double     cMin;
  double     cMax;
  double     mMin;
  double     mMax;
  double     rMin;
  double     rMax;
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



// Stores values of FoxH functions used for generating
//  convergence of Einasto RTS
class einTable {

  public:

    einTable(){

    }


    void setBins  ( int a, int x ) { x_bins = x ;
                                     a_bins = a ;
                                     N_bins =a*x;
                                  initVals()    ; }

    void setX_min ( double     d ) { minX   = d ; }  // Boundaries of the table
    void setX_max ( double     d ) { maxX   = d ; }
    void setA_min ( double     d ) { minA   = d ; }
    void setA_max ( double     d ) { maxA   = d ; }


    double getX_min () { return   minX ; }
    double getX_max () { return   maxX ; }
    double getA_min () { return   minA ; }
    double getA_max () { return   maxA ; }

    int    getA_bins() { return a_bins ; }
    int    getX_bins() { return x_bins ; }
    int    getN_bins() { return N_bins ; }



    // Populate the array
    void setVal   ( int     aBin ,
                    int     xBin ,
                    double value ){

      if ( val != NULL ){

        val[ xBin + aBin * x_bins ] = value;

      }
    }


    double getVal ( int    aBin ,
                    int    xBin ){
      if ( val != NULL ){
        return val[ xBin + aBin * x_bins ];
      }
        return 0;
    }


  private:

    int    a_bins;
    int    x_bins;
    int    N_bins;
    double minA, minX;
    double maxA, maxX;

    double *val = NULL;

    void initVals () { val = new double[ N_bins ]; }

};


// Holds info on density profile, concentration, mass, shape param, etc
class densProfile{
  public:

    densProfile();
    densProfile( double inpA );

    // Modifiers, as parameters are modified need to adjust dependant values
    void setAlpha( double inpA ){         alpha = inpA; }
    void setR_max( double inpR ){         r_max = inpR; }
    void setC    ( double inpC ){ concentration = inpC; }
    void setM_enc( double inpM ){         M_enc = inpM; }
    void setType ( double inpT ){          type = inpT; }


    // If when getting an undefined variable, need to spit out a warning
    double getR_s  () const {  if ( concentration==-1.0 ) unDefVar("\"concentration\"");
                               if (         r_max==-1.0 ) unDefVar("\"R_max  \""      );   return r_max / concentration;  }
    double getC    () const {  if ( concentration==-1.0 ) unDefVar("\"concentration\"");   return         concentration;  }
    double getR_max() const {  if (         r_max==-1.0 ) unDefVar("\"R_max\""        );   return                 r_max;  }
    double getM_enc() const {  if (         M_enc==-1.0 ) unDefVar("\"M_enc\""        );   return                 M_enc;  }
    double getAlpha() const {  if (         alpha==-1.0 ) unDefVar("\"alpha\""        );   return                 alpha;  }
    double getType () const {                                                              return                  type;  }

    double getRho_o() const ;

  private:
    double concentration;
    double alpha        ;
    double r_max        ;
    double M_enc        ;
    int    type         ; //1 NFW, 2 Einasto

    void   unDefVar( std::string inpS ) const {
      std::cerr<<"WARNING: variable "<<inpS<<" undefined in <lensProfile>"<<std::endl;
      exit(0);
    }
};


#endif

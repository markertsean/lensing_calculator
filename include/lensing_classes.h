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

    void setZ      (      double inpZ ){  setAngDist(inpZ);
                                                 z = inpZ; }
    void setM      (      double inpM ){         m = inpM; }
    void setC      (      double inpC ){         c = inpC;
                                          if( rmax!=-1.0)
                                          setRscale(     );}
    void setRmax   (      double inpR ){      rmax = inpR;
                                          if(    c!=-1.0)
                                          setRscale(     );}

    void setLensZ  (      double inpZ ){     lensZ = inpZ; }
    void setCosmo  (   COSMOLOGY inpC ){     cosmo = inpC; }
    void setProfile( std::string inpP ){   profile = inpP; }

    double      getZ        () {
      if (      z==-1.0)
        unDefVar("\"z\"");
      return z;
    }
    double      getM        () {
      if (      m==-1.0)
        unDefVar("\"Mass\"");
      return       m;
    }
    double      getC        () {
      if (      c==-1.0)
        unDefVar("\"concentration\"");
      return       c;
    }
    double      getRmax     () {
      if (   rmax==-1.0)
        unDefVar("Rmax");
      return    rmax;
    }
    double      getRscale   () {
      if (     rs==-1.0)
        unDefVar("\"Rscale\"");
      return      rs;
    }
    double      getAngDist  () {
      if (angDist==-1.0 )
        unDefVar("\"angDist\"");
      return angDist;
    }
    std::string getProfile  () {
      if (profile=="NULL")
        unDefVar("\"profile\"");
      return profile;
    }

    double      getSigmaCrit() {
      if (      sigmaC!=-1.0 )               //If have sigmaC, return it
        return  sigmaC;
      else if ( (lensZ!=-1.0) && (z!=-1.0) ){//If first call, calc sigmaC
        sigmaC = cosmo.SigmaCrit( lensZ, z );
        return sigmaC;
      }
      else{                                   //Otherwise return -1
        std::cerr << "Error: sigmaC requires lensZ & z\n" <<std::endl;
        exit(0);
        return -1.0;
      }
    }

    double  getRealFOV ( double ang ) {
      if (angDist!=-1.0){
        return ang*angDist;
      }
      else{
        return -1.0;
      }
    }

private:
    double       z;             //redshift
    double       m;             //mass Msun
    double       c;             //concentration
    double      rs;             //scale radius of nfw
    double    rmax;             //max radius Mpc
    double   lensZ;             //redshift of lens, if this is a source
    double  sigmaC;             //Sc only for source
    double angDist;             //ang diam distance Mpc


    COSMOLOGY cosmo;
    std::string profile;      //If there is a profile IE NFW

    void setRscale (                  ){      rs = rmax/c             ; }
    void setAngDist(      double inpZ ){ angDist = cosmo.angDist(inpZ); }
    void   unDefVar( std::string inpS ){
      std::cerr<<"WARNING: variable "<<inpS<<" undefined in <haloInfo>"<<std::endl;
    }
};

haloInfo::haloInfo(){

  profile="NULL";
       z=-1.0;
       m=-1.0;
       c=-1.0;
      rs=-1.0;
    rmax=-1.0;
   lensZ=-1.0;
  sigmaC=-1.0;
 angDist=-1.0;

}


//Holds info on density profile, concentration, mass, shape param, etc
class lensProfile{
public:

  inline lensProfile();
  inline lensProfile( double inpA );

  void setAlpha( double inpA ){         alpha = inpA;
                                          modRho_o(); }


  void setR_max( double inpR ){         r_max = inpR;
                                          modR_s  ();
                                          modC    ();
                                          modRho_o(); }

  void setC    ( double inpC ){ concentration = inpC;
                                          modR_s  ();
                                          modRho_o(); }

  void setR_s  ( double inpR ){           r_s = inpR;
                                          modC    ();
                                          modRho_o(); }

  void setRho_o( double inpR ){         rho_o = inpR;
                                          modM_enc(); }
  void setM_enc( double inpM ){         M_enc = inpM;
                                          modRho_o(); }

  void resetM_enc(           ){           modM_enc(); }

  double getC() {
      if ( concentration==-1.0 )
        unDefVar("\"concentration\"");
      return concentration;
  }

  double getRho_o() {
      if ( rho_o==-1.0 )
        unDefVar("\"rho_o\"");
      return rho_o;
  }

  double getR_max() {
      if ( r_max==-1.0 )
        unDefVar("\"R_max\"");
      return r_max;
  }

  double getR_s() {
      if ( r_s==-1.0 )
        unDefVar("\"R_scale\"");
      return r_s;
  }

  double getM_enc() {
      if ( M_enc==-1.0 )
        unDefVar("\"M_enc\"");
      return M_enc;
  }

  double getAlpha() {
      if ( alpha==-1.0 )
        unDefVar("\"alpha\"");
      return alpha;
  }

  double getType() {
      return type;
  }


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

    if ( rho_o > 0 && r_max > 0 && concentration > 0){ //All needed parameters def
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

  void modR_s(){
    if (   concentration > 0 && r_max > 0 )
      r_s = r_max / concentration;
  }

  void modC(){
    if (  r_s > 0 && r_max > 0 )
      concentration = r_max / r_s ;
  }

};

lensProfile::lensProfile(){
  concentration = -1.0;
  alpha         = -1.0;
  rho_o         = -1.0;
  r_s           = -1.0;
  r_max         = -1.0;
  M_enc         = -1.0;
  type          =    1;
}


lensProfile::lensProfile( double inpA ){
  concentration = -1.0;
  alpha         = inpA;
  rho_o         = -1.0;
  r_s           = -1.0;
  r_max         = -1.0;
  M_enc         = -1.0;
  type          =    2;
}

//Holds info on what user wants, whether binning, running model or Nbody, etc
class userInfo{
public:

  inline userInfo();

  double      angFOV;//Angular size of image
  double       R_max;
  int    N_pixels   ;  //Number of pixels on grid
  int    N_bins     ;  //Number of bins for radial averaging
  int    N_sources  ;  //Number of sources to generate
  int    N_particles;  //Number of particles in simulation
  int    num_threads;  //Number of threads for parallel processing
  int          model;  //Flag to do model   analysis
  int          nbody;  //Flag to do nbody   analysis
  int        massMap;  //Flag to do massMap analysis
  int   N_partSmooth;  //Number of particles nearby for glamer smoothing
  int   N_edgepixels;  //Gap in distance between sources
  double nearestSourceNeighbor; //Number of pixels to leave on an edge

  std::string readFile;

  //Chi2 & genetic algorithm fitting values
  double     cMin;
  double     cMax;
  double     mMin;
  double     mMax;
  double alphaMin;
  double alphaMax;

  int maxFitAttempts; //Max attempts at fitting before abort
  int  N_chromosomes; //Number of chromosomes in population
  int     N_chiTrack; //Number of chi's to track for convergence
  int     consistent; //Number of times need avg below tolerance
  double   tolerance; //Average residual must be below tolerance
  double   mutChance; //Likelihood of mutation
  double  avgTestVal; //chiAvg*this is random range
};

userInfo::userInfo(){
      angFOV = -1.;
       R_max = -1.;
 N_pixels    = -1;
 N_bins      = -1;
 N_sources   = -1;
 N_particles = -1;
 num_threads =  1;
       model =  0;
       nbody =  0;
     massMap =  0;
N_partSmooth =  5;
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

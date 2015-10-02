#ifndef MY_LENSING_CLASSES
#define MY_LENSING_CLASSES

#include <iostream>
#include <cosmo.h>
#include <math.h>
#include <astro_constants.h>

//Holds info on lens, source, any halo for easy access
class haloInfo{
public:
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
    double       z=-1.0;             //redshift
    double       m=-1.0;             //mass Msun
    double       c=-1.0;             //concentration
    double      rs=-1.0;             //scale radius of nfw
    double    rmax=-1.0;             //max radius Mpc
    double   lensZ=-1.0;             //redshift of lens, if this is a source
    double  sigmaC=-1.0;             //Sc only for source
    double angDist=-1.0;             //ang diam distance Mpc


    COSMOLOGY cosmo;
    std::string profile="NULL";      //If there is a profile IE NFW

    void setRscale (                  ){      rs = rmax/c             ; }
    void setAngDist(      double inpZ ){ angDist = cosmo.angDist(inpZ); }
    void   unDefVar( std::string inpS ){
      std::cerr<<"WARNING: variable "<<inpS<<" undefined in <haloInfo>"<<std::endl;
    }
};

//Holds info on density profile, concentration, mass, shape param, etc
class lensProfile{
public:

  lensProfile  (             ){          type =    1; }
  lensProfile  ( int inpT    ){          type = inpT; }

  void setAlpha( double inpA ){         alpha = inpA; }


  void setR_max( double inpR ){         r_max = inpR;
                                          modR_s  ();
                                          modC    ();
                                          modM_enc();
                                          modRho_o(); }

  void setC    ( double inpC ){ concentration = inpC;
                                          modR_s  ();
                                          modM_enc();
                                          modRho_o(); }

  void setR_s  ( double inpR ){           r_s = inpR;
                                          modC    ();
                                          modM_enc();
                                          modRho_o(); }

  void setRho_o( double inpR ){         rho_o = inpR;
                                          modM_enc(); }
  void setM_enc( double inpM ){         M_enc = inpM;
                                          modRho_o(); }

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
  double concentration = -1.0;
  double alpha         = -1.0;
  double rho_o         = -1.0;
  double r_s           = -1.0;
  double r_max         = -1.0;
  double M_enc         = -1.0;
  int    type          =    1; //1 NFW, 2 Einasto

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
        rho_o   =   M_enc    *alpha / (4. * M_PI * r_s * r_s * r_s ) /
                  (   exp( 4./alpha ) * pow( alpha/2., 3./alpha)   ) /
                   tgamma( 3./alpha );
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
        M_enc   =  rho_o  * 4. * M_PI * r_s * r_s * r_s *    exp( 4./alpha ) *
                  pow( alpha/2., 3./alpha)   /     alpha * tgamma( 3./alpha );
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


//Holds info on what user wants, whether binning, running model or Nbody, etc
class userInfo{
public:
  double      angFOV = -1.0;//Angular size of image
  int    N_pixels    = -1;  //Number of pixels on grid
  int    N_bins      = -1;  //Number of bins for radial averaging
  int    N_sources   = -1;  //Number of sources to generate
  int    N_particles = -1;  //Number of particles in simulation
  int    num_threads =  1;  //Number of threads for parallel processing
  int          model =  0;  //Flag to do model analysis
  int          nbody =  0;  //Flag to do nbody analysis
  int   N_partSmooth =  5;  //Number of particles nearby for glamer smoothing
  int   N_edgepixels =  3;  //Gap in distance between sources
  double nearestSourceNeighbor = 1.5; //Number of pixels to leave on an edge

  std::string readFile = " ";

  //Chi2 & genetic algorithm fitting values
  double     cMin =  2.0;
  double     cMax =  7.0;
  double     mMin = 13.0;
  double     mMax = 16.0;
  double alphaMin =  0.1;
  double alphaMax =  0.4;

  int maxFitAttempts =  1e3    ; //Max attempts at fitting before abort
  int  N_chromosomes = 1000    ; //Number of chromosomes in population
  int     consistent =   10    ; //Number of times need avg below tolerance
  double   tolerance =    0.001; //Average residual must be below tolerance
  double   mutChance =    0.01 ; //Likelihood of mutation
  double  avgTestVal =    1.3  ; //chiAvg*this is random range
};


#endif

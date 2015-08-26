#include <math.h>
#include <omp.h>
#include <lensing_classes.h>
#include <lens_fitter.h>
#include <my_utilities.h>



//Surface density at radius for NFW profile
double    SDNFW( double r, double r_s, double rho_o ){

  double    x = r/r_s;
  double temp = 2*r_s*rho_o;
  if      ( x < (1.0-1e-6) ){
    return temp / ( x*x - 1 ) * ( 1 - 2./sqrt( 1 - x*x ) *
        atanh( sqrt(( 1 - x ) / ( 1 + x )) ) );
  }
  else if ( x > (1.0+1e-6) ){
    return temp / ( x*x - 1 ) * ( 1 - 2./sqrt( x*x - 1 ) *
        atan ( sqrt(( x - 1 ) / ( 1 + x )) ) );
  }
  else{
    return temp/3.;
  }
}


//Average surface density in radius for NFW
double SDAvgNFW( double r, double r_s, double rho_o ){

  double    x = r/r_s;
  double temp = 4*r_s*rho_o;

  if ( x < (1.0-1e-6) ){
    return temp / (x*x) * ( 2./sqrt(   1 - x*x ) *
                        atanh( sqrt((  1 - x )/( 1 + x )) ) + log(x/2.));
  }
  else if ( x > (1.0+1e-6) ){
    return temp / (x*x) * ( 2./sqrt( x*x - 1 ) *
                        atan ( sqrt((  x - 1 )/( 1 + x )) ) + log(x/2.));
  }
  else{
    return temp * ( 1 + log(0.5) );
  }
}


/*
Generates the radially averaged reduced tangential shear for
an NFW profile for given input
*/
void generateNFWRTS( double *gArr, lensProfile &lens, haloInfo &halo,                     userInfo   u, double *sourceSc, double *sourceDist ){

  double rho_o = lens.getM_enc() / (4. * M_PI * pow( halo.getRmax(), 3 ) ) *
                  pow(       lens.getC(), 3 ) /
                  ( log( 1 + lens.getC() )    - lens.getC() / ( 1 + lens.getC() ));

  int      N_bin[u.N_bins];
  double analRTS[u.N_bins];
  //Set avg array values to 0 initially, will be adding upon
  for (int i=0;i<u.N_bins;++i){
      N_bin[i] = 0;
    analRTS[i] = 0;
  }

  for (int i=0;i<u.N_sources;++i)
    gArr[i] = 0;

  //Distance converted to bins
  double stepSize = halo.getRealFOV( u.angFOV )/2.0/u.N_bins;

  //Loop over all sources, determining predicted rts for a given source
  for ( int i=0; i < u.N_sources; ++i ){
    double    SD =    SDNFW( sourceDist[i], lens.getR_s(), rho_o ); //At radius
    double avgSD = SDAvgNFW( sourceDist[i], lens.getR_s(), rho_o ); //Average
    double SigCr = sourceSc[i];

    int   binNum = std::min(std::max(
                  (int) round(sourceDist[i]/stepSize) ,0),u.N_bins-1);

    //RTS avg across bin
    analRTS[binNum] += ( avgSD - SD ) / ( SigCr - SD );
      N_bin[binNum] += 1;
  }
  //Returning gArr
  for ( int i=0; i < u.N_sources; ++i ){
    if ( N_bin[i] > 0) gArr[i] = analRTS[i] / N_bin[i];
  }
}

/*
  From Klypin 2014
  Redshift C_o gamma M_o/10^12
  0.00    7.40 0.120 5.5e5
  0.35    6.25 0.117 1.0e5
  0.50    5.65 0.115 2.0e4
  1.00    4.30 0.110 9.5e2

  C(M)= C_o (M/10^12)^-gamma [ 1 + (M/M_o)^0.4 ]

  z=0 M=13 C=5.68
  z=0 M=14 C=4.39
  z=0 M=15 C=3.49
  z=1 M=13 C=3.88
  z=1 M=14 C=3.64
  z=1 M=15 C=4.06

  NFW
  Mvir=M200 = 4 pi rho_o Rv^3 c^-3 [ ln(1+c) - c/1-c ]
  rho_o = Rvir^-3 c^3 /4pi[] * M

  C = 2-7

  M = 10^13 - 10^16

  alpha = 0.1-0.4

//*/
/*
  Takes input sources, and fits density profile to the rts values.
  Can do NFW, needs modifications for Einasto
  Potential to fit NFW to < 1%
*/
void fitDensProfile( lensProfile &densProfile, haloInfo &halo, userInfo u,
                double *gArr,  double *dArr, double *gErrArr, double *sourceSigC,
                double *sourceDist ){


  //Values used for comparing goodness of fits
  double        avgChi(0), oldAvg(0), testVal(0);

  int numGenes = 2;                  //2 variables for NFW
  if ( densProfile.getType() == 2.0)
      numGenes = 3;                  //3 variables for Einasto

  //The genetic algorithm part, and chi2 fitting stuffs
  lensProfile         parent[ u.N_chromosomes], child[ u.N_chromosomes ];
  double      gAnalyticArray[ u.N_bins       ],  chi2[ u.N_chromosomes ];

  double avgChiThreads[ u.num_threads ];
  for (int i = 0; i  <  u.num_threads  ; ++i)
    avgChiThreads[i]=0;




  //Generate first parent values
  #pragma omp parallel for
  for ( int i=0; i < u.N_chromosomes; ++i ){

    //Generate initial parent values
    parent[i].setC    (         randVal( u.cMin, u.cMax )  );
    parent[i].setM_enc(pow( 10, randVal( u.mMin, u.mMax ) ));
    parent[i].setR_s  ( halo.getRmax() / parent[i].getC()  );

    //Gives analytic value for NFW halo, for source locations, radially averaged
    generateNFWRTS( gAnalyticArray, parent[i], halo, u, sourceSigC, sourceDist);

    chi2[i] = chiSquared( gAnalyticArray, gArr, gErrArr, u.N_bins );

    avgChiThreads[omp_get_thread_num()]+=chi2[i];
  }
  //Avg used for reproduction, randomly sample halos based on goodness / 1.5avg
  for (int i  = 0; i < u.num_threads; ++i){
    avgChi           += avgChiThreads[i];
    avgChiThreads[i]  = 0;
  }
  avgChi  = avgChi / u.N_chromosomes;

  int loopCounter(0); //Counter of number of iterations, over 1e6 then stop
  int     counter(0); //Counter used to count number of times within tolerance
  double   totAvg(0); //Keeps running average, stop when converges

  do {
    //Reset chi^2 values
    //Tot average is a running average over all runs, eventually should converge
    // to a fairly constant value
    //Old average saves the previous tot, to measure how much it is changing
    oldAvg  =   totAvg;
    totAvg  = ( totAvg * loopCounter + avgChi) / ( loopCounter + 1 ) ;
    testVal =   avgChi * u.avgTestVal;
    avgChi  =        0;


    //Generate new children
    #pragma omp parallel for
    for ( int i=0; i < u.N_chromosomes; ++i ){

      int parentIndex[2];

      //Select mates from parent pool
      for ( int j=0; j < 2; ++j ){
        do{
          //Randomly select an index
          parentIndex[j] = round( randVal( 0, u.N_chromosomes-1) );
          //If parent fit good enough (small) select w/ uniform distribution
        } while ( chi2[parentIndex[j]] > randVal( 0, testVal) );
      }


      int aIndex, cIndex, mIndex;

      //Get indexes to use from parents, since parents are random doesn't
      //  matter which we select from, but if Einasto need to decide
      //  whether alpha coming from parent 1 or 2
      if ( densProfile.getType() == 2.0 )
      aIndex = parentIndex[ (int) round( randVal( 0, 1) )  ];
      cIndex = parentIndex[0];
      mIndex = parentIndex[1];


      //Set child values from parents
      child[i].setC    ( parent[cIndex].getC()     );
      child[i].setM_enc( parent[mIndex].getM_enc() );


      //Mutate values, with random chance. Needed to keep genes fresh
      if ( u.mutChance > randVal(0,1) )
        child[i].setC    ( child[i].getC()     * randVal( 0.9, 1.1 )  );

      if ( u.mutChance > randVal(0,1) )
        child[i].setM_enc( child[i].getM_enc() * randVal( 0.9, 1.1 ) );

      child[i].setR_s  ( halo.getRmax()/child[i].getC()   );


      //Generate analytic NFW for child
      generateNFWRTS( gAnalyticArray, child[i], halo, u, sourceSigC, sourceDist);

      //chi2 for the children
      chi2[i] = chiSquared( gAnalyticArray, gArr, gErrArr, u.N_bins );

      avgChiThreads[ omp_get_thread_num() ] += chi2[i];
    }

    //Avg used for reproduction, randomly sample halos based on goodness / 1.5avg
    for (int i  = 0; i < u.num_threads; ++i){
      avgChi           += avgChiThreads[i];
      avgChiThreads[i]  = 0;
    }
    avgChi  =  avgChi / u.N_chromosomes;

    //Children replace parents
    for ( int i=0; i<u.N_chromosomes; ++i ){
      parent[i] = child[i];
    }

    //Need counter to be greater than u.consistent,
    // a count of how many times average has been
    // consistently below tolerance. Tests for convergence
    if (fabs(totAvg-oldAvg)/oldAvg < u.tolerance){
      counter += 1;
    }
    else{
      counter  = 0;
    }
    ++loopCounter;

  } while ( u.consistent > counter && loopCounter < u.maxFitAttempts ) ;

  int    minIndex =              0; //index of lowest chi2
  double minChi   = chi2[minIndex];

  //Find lowest chi2 index
  for ( int i = 1; i < u.N_chromosomes; ++i ){
    if ( minChi > chi2[i] ){
      minIndex = i;
      minChi   = chi2[i];
    }
  }

  densProfile = child[minIndex];
}

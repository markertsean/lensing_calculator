#include <math.h>
#include <cmath>
#include <omp.h>
#include <astro_constants.h>
#include <lensing_classes.h>
#include <lens_fitter.h>
#include <my_utilities.h>


//   2 0[              (0,1)     ]    1     ^ G(2/as)G(s-1/2)           -2s         (-1)^k G(1/a-2/a k)  2k        a(-1)^k G(-1/2-a/2 k)  ak+1
// H    [                        ]= -----  |  ------------------------ X   ds=  Sum ------ ------------ X    + Sum ------- ------------- X
//   1 2[(0,2/a)(-1/2,1)         ]  2pi i U            G( s )                   k=0 (k  )! G(1/2-    k)        k=1 2(k)!   G(    -a/2 k)
double foxH2012(
                double         z ,  // Z from fox H function
                double     alpha ,  // Shape parameter Ein profile
                double tolerance ){ // Tolerance level for convergence

  double sum1(0), sum2(0), oldSum(0), totSum(0), s1(0), s2(0);

  int converge(0), k(-1);

  //k = 0 term
//     sum1 += tgamma( 1. / alpha ) / tgamma( 0.5 );
//  oldSum  += sum1;
//printf("%12.3e,%12.3e %12.3e,%12.3e %12.3e,%12.3e    %12.3e\n",s1,sum1,s2,sum2,totSum,oldSum,fabs((totSum-oldSum)/oldSum));
  do {
      ++k;

      //Terms in summation
      s1 =    pow(       - 1     , k ) / factorial(k) *    pow( z, 2 * k ) *
           tgamma(  ( 1. - 2.    * k ) / alpha      ) / tgamma(  0.5 - k );
      s2 =    pow(       - 1     , k ) / factorial(k) *    pow( z, 1 + k * alpha ) *
           tgamma( -( 1. + alpha * k ) / 2.         ) / tgamma(      - k * alpha / 2.0);

      if ( ! (std::isinf( s1 ) || s1!=s1) )  sum1 += s1;
      if ( ! (std::isinf( s2 ) || s2!=s2) )  sum2 += s2;

      //totSum is term for current k
      totSum = sum1 + alpha/2. * sum2;
      //oldSum an average over all old sums
      oldSum = ( totSum + oldSum * k ) / ( k + 1 );
//printf("%12.3e,%12.3e %12.3e,%12.3e %12.3e,%12.3e    %12.3e\n",s1,sum1,s2,sum2,totSum,oldSum,fabs((totSum-oldSum)/oldSum));
      if ( fabs((totSum-oldSum)/oldSum) < tolerance ){
        converge +=1;
      }
      else{
        converge = 0;
      }
  } while ( converge < 10 && k < 1e2 );

  sum2 *= alpha/2.;

  return sum1+sum2;
}



//  2 1[     (-1/2,1) (0,1)     ]    1     ^ G(2/as)G(s-1/2) G(3/2-s)  -2s        (-1)^k G(1/a-2/a k)  2k+2      a(-1)^k G(-3/2-a/2 k)  ak+3
// H   [                        ]= -----  |  ------------------------ X   ds= Sum ------ ------------ X    - Sum ------- ------------- X
//  2 3[(0,2/a)(-1/2,1)(-3/2,1) ]  2pi i U   G(5/2-s) G( s )                  k=0 (k+1)! G(1/2-1/a  )        k=1 2(k)!   G(    -a/2 k)
double foxH2123(
                double         z ,  // Z from fox H function
                double     alpha ,  // Shape parameter Ein profile
                double tolerance ){ // Tolerance level for convergence

  double sum1(0), sum2(0), oldSum(0), totSum(0), s1(0), s2(0);

  int converge(0), k(-1);

  //k = 0 term
//     sum1 += tgamma( 1. / alpha ) / tgamma( 0.5 - 1. / alpha) * ( z * z );
//  oldSum  += sum1;
//printf("%12.3e,%12.3e %12.3e,%12.3e %12.3e,%12.3e    %12.3e\n",s1,sum1,s2,sum2,totSum,oldSum,fabs((totSum-oldSum)/oldSum));
  do {
      ++k;

      //Terms in summation
      s1 =    pow(       - 1     , k ) / factorial(k+1) *    pow( z, 2 + k * 2     ) *      // (-1)^k/(k+1)! x^(2k + 2)
           tgamma(  ( 1. - 2.    * k ) / alpha        ) / tgamma(  0.5 - 1 / alpha );       // G(  (1-2k)/a ) / G( 1/2 - 1/a )

      s2 =    pow(       - 1     , k ) / factorial(k  ) *    pow( z, 3 + k * alpha ) *      // (-1)^k/(k  )! x^(ak + 3)
           tgamma( -( 3. + alpha * k ) / 2.           ) / tgamma(      - k * alpha / 2.0);  // G( -(3+ak)/2 ) / G( -ak/2 )

      if ( ! (std::isinf( s1 ) || s1!=s1) )  sum1 += s1;
      if ( ! (std::isinf( s2 ) || s2!=s2) )  sum2 += s2;

      //totSum is term for current k
      totSum = sum1 - alpha/2. * sum2;
      //oldSum an average over all old sums
      oldSum = ( totSum + oldSum * k ) / ( k + 1 );
//printf("%12.3e,%12.3e %12.3e,%12.3e %12.3e,%12.3e    %12.3e\n",s1,sum1,s2,sum2,totSum,oldSum,fabs((totSum-oldSum)/oldSum));
      if ( fabs((totSum-oldSum)/oldSum) < tolerance ){
        converge +=1;
      }
      else{
        converge = 0;
      }
  } while ( converge < 10 && k < 1e2 );

  sum2 *= alpha/2.;

  return sum1-sum2;
}


void generateEinRTS(
                    double          *gArr ,  // RTS array to output
                    lensProfile     &lens ,  // Input density profile to generate profile for
                    haloInfo        &halo ,  // Actual information from the halo
                    userInfo            u ,  // Information from the user
                    double    *sourceSc   ,  // Critical surface density for the sources
                    double    *sourceDist ){ // Projected distances between source and lens

  //Arrays for binning
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


  //Constant part of kappa_c, divided by gamma(1/alpha) * sqrt pi, a constant for easier kappas
  double temp = lens.getRho_o() * lens.getR_s  ()   * exp( 2./lens.getAlpha() )       *
                             pow( lens.getAlpha()/2.,      1./lens.getAlpha() - 1.0 ) *
                             sqrt( M_PI );

  for ( int i = 0; i<u.N_bins; ++i )
  {
    double x        = sourceDist[i] / lens.getR_s();

    //Modified kappa_c, kappa_c * sqrt(pi) / Gamma(1/alpha), just need to multiply by H function
    double kappa_cM = temp / sourceSc[i];

    double kappa    = kappa_cM * foxH2012( x*x, lens.getAlpha() );
    double kappaAVG = kappa_cM * foxH2123( x*x, lens.getAlpha() );

    int    binNum   = std::min(std::max(
                      (int) round(sourceDist[i]/stepSize) ,0),u.N_bins-1);

    //RTS avg across bin
    analRTS[binNum] += ( kappaAVG - kappa ) / ( 1 - kappa );
      N_bin[binNum] += 1;

//printf("%12.3e %12.3e %12.3e\n", kappa, kappaAVG, (kappaAVG-kappa)/(1-kappa));
  }

  //Returning gArr
  for ( int i=0; i < u.N_sources; ++i ){
    if ( N_bin[i] > 0) gArr[i] = analRTS[i] / N_bin[i];
  }

  //Kappas are fox H functions

  // |   | m, n [ (a_1,A_1) ... (a_n,A_n) ... (a_p,A_p)  |    ]
  // |---|      [                                        |  Z ]
  // |   | p, q [ (b_1,B_1) ... (b_m,B_m) ... (b_q,B_q)  |    ]

  //            ^ Pj=  1, m Gamma(b_j + B_j s) Pj=  1, n Gamma(1 - a_j + A_j s)
  //  1        |  -------------------------------------------------------------  Z^-s ds
  // ----      |  Pj=m+1, p Gamma(a_j + A_j s) Pj=m+1, q Gamma(1 - b_j + A_j s)
  // 2ipi     U L

  //              2 0[              (0,1)     ]  const   ^ G(2/as)G(s-1/2)           -2s         [     (-1)^k G(1/a-2/a k)  2k        a(-1)^k G(-1/2-a/2 k)  ak+1 ]
  // k = const* H    [                        ]= -----  |  ------------------------ X   ds= const[ Sum ------ ------------ X    + Sum ------- ------------- X     ]
  //              1 2[(0,2/a)(-1/2,1)         ]  2pi i U            G( s )                       [ k=0 (k  )! G(1/2-    k)        k=1 2(k)!   G(    -a/2 k)       ]

  // _            2 1[     (-1/2,1) (0,1)     ]  const   ^ G(2/as)G(s-1/2) G(3/2-s)  -2s         [     (-1)^k G(1/a-2/a k)  2k+2      a(-1)^k G(-3/2-a/2 k)  ak+3 ]
  // k = const* H    [                        ]= -----  |  ------------------------ X   ds= const[ Sum ------ ------------ X    - Sum ------- ------------- X     ]
  //              2 3[(0,2/a)(-1/2,1)(-3/2,1) ]  2pi i U   G(5/2-s) G( s )                       [ k=0 (k+1)! G(1/2-1/a  )        k=1 2(k)!   G(    -a/2 k)       ]
}


//Surface density at radius for NFW profile integrated to infinity along LOS
double    SDNFWFull(
                double   r   ,  // Input radius
                double   r_s ,  // Scale radius, r_-2
                double rho_o ){ // Initial density

  double    x = r / r_s         ;
  double temp = 2 * r_s * rho_o ;

  if      ( x < (1.0-1e-6) ){
    return temp / ( x*x - 1 ) * ( 1 - 2./sqrt(  1  - x*x ) *
        atanh( sqrt(( 1 - x ) / ( 1 + x )) ) );

  }
  else if ( x > (1.0+1e-6) ){
    return temp / ( x*x - 1 ) * ( 1 - 2./sqrt( x*x -  1  ) *
        atan ( sqrt(( x - 1 ) / ( 1 + x )) ) );

  }
  else{
    return temp/3.;
  }
}

//Surface density at input radius for NFW profile, integrated to R_max
double    SDNFW(
                double               r ,  //Input radius to calc SD at
                lensProfile inpProfile ){ //Input NFW profile

  double      x = fabs( r / inpProfile.getR_s() );
  double      c =           inpProfile.getC  ()  ;
  double factor =       2 * inpProfile.getR_s() * inpProfile.getRho_o();

  //Outside of our integration radius
  if ( x > inpProfile.getC() ){
    return 0;
  }
  //If within r_s, arctanh solution has imaginary component we ignore
  else if ( x < 1 ){
    return factor/( x*x - 1 ) * ( sqrt( c*c - x*x ) / ( c + 1 ) - 0.5 / sqrt( 1 - x*x ) * ( log( ( 1. / c * sqrt( ( c*c - x*x )/( 1. - x*x )   ) + 1 ) /
                                                                                                 ( 1. / c * sqrt( ( c*c - x*x )/( 1. - x*x )   ) - 1 ) ) -
                                                                                            log( (          sqrt( ( c*c - x*x )/( 1. - x*x )   ) + 1 ) /
                                                                                                 (          sqrt( ( c*c - x*x )/( 1. - x*x )   ) - 1 ) )   ) ); //+i pi to the logs
  }
  //Outside of r_s
  else if ( x > 1 ) {
    return factor/( x*x - 1 ) * ( sqrt( c*c - x*x ) / ( c + 1)  + 1.0 / sqrt( x*x - 1 ) * ( atan(  1. / c * sqrt( ( c*c - x*x )/( x*x - 1. )         ) ) -
                                                                                            atan(           sqrt( ( c*c - x*x )/( x*x - 1. )         ) )   ) );
  }
  //If close enough to r_s
  else {
    return factor/3.*pow( c*c - 1., -1.5 ) * ( c * ( c*c - 1) - 2 * c + 2 );
  }

}


//Surface density at input radius for NFW profile, integrated to R_max
double    SDAvgNFW(
                double               r ,  //Input radius to calc SD at
                lensProfile inpProfile ){ //Input NFW profile

  double      x = fabs( r / inpProfile.getR_s() );
  double      c =           inpProfile.getC  ()  ;
  double factor =       4 * inpProfile.getR_s() * inpProfile.getRho_o() / ( x*x );

  //Outside of our integration radius, integral from 0 to 1 + integral from 1 to max radius, and take out their factors to use the current one
  if ( x > inpProfile.getC() ){
    return factor * ( SDAvgNFW( inpProfile.getR_s  (), inpProfile ) / ( 4 * inpProfile.getR_s  () * inpProfile.getRho_o() )
                    + SDAvgNFW( inpProfile.getR_max(), inpProfile ) / ( 4 * inpProfile.getR_max() * inpProfile.getRho_o() / ( inpProfile.getC() * inpProfile.getC() ) ) );

  }

  //If within r_s, this is analytic solution integrating to a max value of the concentration
  else if ( x < 1 ){

    return factor * ( ( sqrt( c*c - x*x )  - c )/( c + 1 )

                    +         log( c + 1 )

                    - atanh( sqrt( 1 - x*x/( c*c )) )

                    + 0.5 /  sqrt( 1 - x*x ) * ( log( ( 1. / c * sqrt( ( c*c - x*x )/( 1. - x*x )   ) + 1 )
                                             /        ( 1. / c * sqrt( ( c*c - x*x )/( 1. - x*x )   ) - 1 ) )
                                             -   log( (          sqrt( ( c*c - x*x )/( 1. - x*x )   ) + 1 )
                                             /        (          sqrt( ( c*c - x*x )/( 1. - x*x )   ) - 1 ) )   ) );
  }
  //Outside of r_s, integral from 0 to 1 + new component to a max of rmax/rs
  else if ( x > 1 ){
    return factor * (  SDAvgNFW( inpProfile.getR_s(), inpProfile  ) / ( 4 * inpProfile.getR_s() * inpProfile.getRho_o() )
                    + ( ( sqrt( c*c - x*x ) - 2 * sqrt( c*c - 1 ) ) / ( c + 1 )

                    + 0.5 * log( ( 1. / c * sqrt( c*c -  1  ) + 1 ) / ( 1 - 1. / c * sqrt( c*c -  1  ) ) )
                    - 0.5 * log( ( 1. / c * sqrt( c*c - x*x ) + 1 ) / ( 1 - 1. / c * sqrt( c*c - x*x ) ) )

                    - 1.0 / sqrt( x*x - 1 ) * ( atan(  1. / c * sqrt( ( c*c - x*x )/( x*x - 1. )         ) )
                                            -   atan(           sqrt( ( c*c - x*x )/( x*x - 1. )         ) )   ) ));
  }
  //If close enough to r_s will integrate to 1, hyperbolic previous arctan component turns into sqrt component
  else {
    return factor * ( ( 2 * sqrt( c*c - 1 ) - c ) / ( c + 1 )

                  +         log( c + 1 )

                  - atanh( sqrt( 1 - x*x / ( c*c ) ) )  );
  }

}


//Average surface density within radius for NFW
double SDAvgNFWFull(
                double     r ,  // Input radius
                double   r_s ,  // Scale radius, r_-2
                double rho_o ){ // Initial density

  double    x = r / r_s         ;
  double temp = 4 * r_s * rho_o ;

  if ( x < (1.0-1e-6) ){
    return temp / (x*x) * ( 2./sqrt(  1 - x*x )*
                        atanh( sqrt(( 1 - x   )/( 1 + x )) ) + log(x/2.));
  }

  else if ( x > (1.0+1e-6) ){
    return temp / (x*x) * ( 2./sqrt(  x*x - 1 )*
                        atan ( sqrt((   x - 1 )/( 1 + x )) ) + log(x/2.));
  }

  else{
    return temp * ( 1.0 + log(0.5) );
  }
}


/*
Generates the radially averaged reduced tangential shear for
an NFW profile for given input
*/
void generateNFWRTS(
                    double          *gArr ,  // RTS array to output
                    lensProfile     &lens ,  // Input density profile to generate profile for
                    haloInfo        &halo ,  // Actual information from the halo
                    userInfo            u ,  // Information from the user
                    double    *sourceSc   ,  // Critical surface density for the sources
                    double    *sourceDist ){ // Projected distances between source and lens

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
    double    SD =    SDNFW( sourceDist[i], lens ); //At radius
    double avgSD = SDAvgNFW( sourceDist[i], lens ); //Average
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
void fitDensProfile(
                    lensProfile &densProfile ,  // Density profile we are outputting
                    haloInfo    &       halo ,  // Info about parent halo
                    userInfo               u ,  // Info from the user
                    double       *      gArr ,  // RTS binned array we "observed"
                    double       *      dArr ,  // Distance binned array
                    double       *   gErrArr ,  // Error array in RTS
                    double       *sourceSigC ,  // Crit surface densities of sources
                    double       *sourceDist ){ // 2D distance from source to lens


  //Values used for comparing goodness of fits
  double        avgChi(0), oldAvg(0), testVal(0);

  int numGenes = 2;                  //2 variables for NFW
  if ( densProfile.getType() == 2.0)
      numGenes = 3;                  //3 variables for Einasto

  //The genetic algorithm part, and chi2 fitting stuffs
  lensProfile         parent[ u.N_chromosomes], child[ u.N_chromosomes ];
  double      gAnalyticArray[ u.N_bins       ],  chi2[ u.N_chromosomes ];



  //Initialize average arrays for diff threads
  double avgChiThreads[ u.num_threads ];
  for (int i = 0; i  <  u.num_threads   ; ++i)
    avgChiThreads[i]=0;

  //Set all parent R_maxes to the lens
  for (int i = 0; i  <  u.N_chromosomes ; ++i){
    parent[i].setR_max( halo.getRmax() );
     child[i].setR_max( halo.getRmax() );
  }

  //Generate first parent values
  #pragma omp parallel for
  for ( int i=0; i < u.N_chromosomes; ++i ){

    //Generate initial parent values, if Ein profile need alpha parameter
    if ( densProfile.getType() == 2 )
    parent[i].setAlpha(         randVal( u.alphaMin, u.alphaMax )  );
    parent[i].setC    (         randVal(     u.cMin,     u.cMax )  );
    parent[i].setM_enc(pow( 10, randVal(     u.mMin,     u.mMax ) ));



    //Gives analytic value for RTS, for source locations, radially averaged
    if ( densProfile.getType() == 1 ){
      generateNFWRTS( gAnalyticArray, parent[i], halo, u, sourceSigC, sourceDist);
    } else{
      generateEinRTS( gAnalyticArray, parent[i], halo, u, sourceSigC, sourceDist);
    }

    //Calc chi2 between theoretical predictions and real values
    chi2[i] = chiSquared( gAnalyticArray, gArr, gErrArr, u.N_bins );
    avgChiThreads[omp_get_thread_num()]+=chi2[i];
  }

  //Avg used for reproduction, randomly sample halos based on goodness / 1.5avg
  for (int i  = 0; i < u.num_threads; ++i){
    avgChi           += avgChiThreads[i];
    avgChiThreads[i]  = 0;
  }
  avgChi  = avgChi / u.N_chromosomes;

  //Reproduction time
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
      cIndex = parentIndex[0];
      mIndex = parentIndex[1];

      //If Einasto, another parameter needed
      if ( densProfile.getType() == 2 ){
      aIndex = parentIndex[ (int) round( randVal( 0, 1) )  ];


      child[i].setAlpha( parent[aIndex].getAlpha() );
      }

      //Set child values from parents
      child[i].setC    ( parent[cIndex].getC()     );
      child[i].setM_enc( parent[mIndex].getM_enc() );




      //Mutate values, with random chance. Needed to keep genes fresh
      if ( u.mutChance > randVal(0,1) )
        child[i].setC    ( child[i].getC()     * randVal( 0.9, 1.1 )  );

      if ( u.mutChance > randVal(0,1) )
        child[i].setM_enc( child[i].getM_enc() * randVal( 0.9, 1.1 ) );




      //Generate analytic RTS for child
      if ( densProfile.getType() == 1 ){
        generateNFWRTS( gAnalyticArray, child[i], halo, u, sourceSigC, sourceDist);
      } else{
        generateEinRTS( gAnalyticArray, child[i], halo, u, sourceSigC, sourceDist);
      }


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
//printf("%5i %12.3e         %12.3e %12.3e       %12.3e %5i\n",loopCounter,avgChi,totAvg,oldAvg,fabs(totAvg-oldAvg)/oldAvg,counter);
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

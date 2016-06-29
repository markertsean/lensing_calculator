#include <math.h>
#include <cmath>
#include <omp.h>
#include <astro_constants.h>
#include <lensing_classes.h>
#include <lens_fitter.h>
#include <my_utilities.h>



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
//
// Does the rolling part of the lens fitting
//
void rollBall(        densProfile   &ball ,  // Ball to roll
                      double        &chi2 ,  // Chi2 value
               const  double        *gArr ,  // g values we are fitting
               const  double        *dArr ,  // distance array corresponding to g values
               const  double     *gErrArr ,  // Errors associated with g
               const  double       sigmaC ,  // Critical surface density to use, of sources
               const  userInfo          u ){


  double smallStep = 200;  // Number of possible steps in a direction
  double   bigStep = 100;

  double decrement = 1/4.; // Amount to decrease step size by


  // Step sizes partially random
  double cStep = ( u.getConMax  () - u.getConMin  () ) / randVal( bigStep, smallStep ) ;
  double mStep = ( u.getMassMax () - u.getMassMin () ) / randVal( bigStep, smallStep ) ;
  double aStep = ( u.getAlphaMax() - u.getAlphaMin() ) / randVal( bigStep, smallStep ) ;

  // For kicking the ball back
  double cStep_0 = cStep;
  double mStep_0 = mStep;
  double aStep_0 = aStep;


  // Previous values, to test for convergence
  double prevM=1;
  double prevC=1;
  double prevA=1;


  int loopCounter = 0;                      // Number of times we have rolled
  int convCounter = 0;                      // Number of times we have been in tolerance

  // Roll the ball
  do{

      // Used to check which direction is best
      densProfile testProfile;
      testProfile.setR_max( ball.getR_max() );
      double mVals[3], cVals[3], aVals[3];


      // Sets the parameters for moving along the hill
      for ( int j = -1; j < 2; ++j ){
        aVals[j+1] =          std::min(std::max(        ball.getAlpha()   + aStep * j, u.getAlphaMin() ), u.getAlphaMax() );
        mVals[j+1] = pow( 10, std::min(std::max( log10( ball.getM_enc() ) + mStep * j, u.getMassMin () ), u.getMassMax () ) );
        cVals[j+1] =          std::min(std::max(        ball.getC    ()   + cStep * j, u.getConMin  () ), u.getConMax  () );
      }

      double chiChoose(999), compChi(998);

      double gAnalyticArray[  u.getNbins()   ];

      for ( int ii = 0; ii < 3; ++ii ){

        testProfile.setC(     cVals[ii] );

      for ( int jj = 0; jj < 3; ++jj ){

        testProfile.setM_enc( mVals[jj] );

//      generateEinRTS( gAnalyticArray, parent[i], halo, u, sourceSigC, sourceDist);
// Add Einasto
        generateNFWRTS( gAnalyticArray, testProfile, u.getNbins(), dArr, sigmaC );
        compChi = chiSquared( gAnalyticArray, gArr, gErrArr, u.getNbins() );

        // If a better fit, keep it
        if (compChi < chiChoose){
          chiChoose = compChi;
          ball.setC    ( cVals[jj] );
          ball.setM_enc( mVals[jj] );

          chi2 = compChi;
        }


      }
      }

      // Kick ball if stuck on edge

      if ( fabs( ball.getC    () - u.getConMax () ) / u.getConMax () < 1e-4 ||
           fabs( ball.getC    () - u.getConMin () ) / u.getConMin () < 1e-4 ){
           ball.setC(              randVal( u.getConMin  (), u.getConMax  () ) );
           cStep = cStep_0;
      }

      if ( fabs( ball.getM_enc() - u.getMassMax() ) / u.getMassMax() < 1e-4 ||
           fabs( ball.getM_enc() - u.getMassMin() ) / u.getMassMin() < 1e-4 ){
           ball.setM_enc( pow( 10, randVal( u.getMassMin (), u.getMassMax () ) ) );
           mStep = mStep_0;
      }
// Einasto


      // If didn't move in one direction, decrease step size
      if ( ball.getC    () == cVals[1] &&
           ball.getM_enc() == mVals[1] ){
        cStep *= decrement;
        mStep *= decrement;
      }
// Add einasto




      // Need counter to be greater than u.consistent,
      //  a count of how many times average has been
      //  consistently below tolerance. Tests for convergence
      if (  fabs( ball.getM_enc() - prevM ) / prevM < u.getTolerance() &&
            fabs( ball.getC    () - prevC ) / prevC < u.getTolerance() ){
//Einasto
        convCounter += 1;
      } else {
        convCounter  = 0;
      }



// Einasto
      prevM = ball.getM_enc();
      prevC = ball.getC    ();

      loopCounter +=1;
    } while ( loopCounter < u.getMaxFitNum()   &&
              convCounter < u.getNConsistent() );

}


//
//  Takes input sources, and fits density profile to the rts values.
//  Can do NFW, needs modifications for Einasto
//  Potential to fit NFW to < 1%
//
void rollingFitDensProfile(
                          densProfile   &  profile ,  // Density profile we are outputting
                    const haloInfo      &     halo ,  // Info about parent halo
                    userInfo               u ,  // Info from the user
//                    const userInfo               u ,  // Info from the user
                    const double       *      gArr ,  // RTS binned array we "observed"
                    const double       *      dArr ,  // Distance binned array
                    const double       *   gErrArr ,
                    const COSMOLOGY          cosmo ){ // Error array in RTS


  densProfile ball[ u.getNchrome() ];
  double      chi2[ u.getNchrome() ];

  for (int i = 0; i  <  u.getNchrome() ; ++i)
    ball[i].setR_max( halo.getRmax() );


  // Do for each rolling ball
  #pragma omp parallel for
  for ( int i = 0; i < u.getNchrome(); ++i ){


    // Set starting parameters, location on hill
    if ( profile.getType() == 2 );
    ball[i].setAlpha(          randVal( u.getAlphaMin(), u.getAlphaMax() )   );
    ball[i].setM_enc( pow( 10, randVal( u.getMassMin (), u.getMassMax () ) ) );
    ball[i].setC    (          randVal( u.getConMin  (), u.getConMax  () )   );

    rollBall( ball[i], chi2[i], gArr, dArr, gErrArr, cosmo.SigmaCrit( halo.getZ(), u.getSourceZ() ), u );

  }

  int    minIndex =              0; //index of lowest chi2
  double minChi   = chi2[minIndex];
  double cAvg(0), mAvg(0), weightedChi(0);

  //Find lowest chi2 index
  for ( int i = 0; i < u.getNchrome(); ++i ){
    if ( minChi > chi2[i] ){
      minIndex = i;
      minChi   = chi2[i];
    }
//         avgArr[iBin] += avgVal / n_srcs / ( errors[i]*errors[i] ); // Weighted average sum

    weightedChi += 1/chi2[i] ;

    cAvg +=             ball[i].getC()       /  chi2[i]  ;
    mAvg += std::log10( ball[i].getM_enc() ) /  chi2[i]  ;
//Einasto
  }

  weightedChi = 1 / weightedChi;
  cAvg = cAvg * weightedChi ;
  mAvg = mAvg * weightedChi ;

  profile.setC    (  ball[minIndex].getC    ()  );
  profile.setM_enc(  ball[minIndex].getM_enc()  );
//Einasto

printf("\n");
//printf("%14.4e %7.5f %14.4e\n", weightedChi, cAvg, pow(10, mAvg ) );
printf("%14.4e ", minChi);


}



//
//  Takes input sources, and fits density profile to the rts values.
//  Can do NFW, needs modifications for Einasto
//  Potential to fit NFW to < 1%
//
void fitDensProfile(
                          densProfile   &  profile ,  // Density profile we are outputting
                    const haloInfo      &     halo ,  // Info about parent halo
                    userInfo               u ,  // Info from the user
//                    const userInfo               u ,  // Info from the user
                    const double       *      gArr ,  // RTS binned array we "observed"
                    const double       *      dArr ,  // Distance binned array
                    const double       *   gErrArr ,
                    const COSMOLOGY          cosmo ){ // Error array in RTS

  // Values used for comparing goodness of fits
  double  avgChi(0), oldAvg(0), testVal(0);

  int numGenes = 2;                  //2 variables for NFW
  if ( profile.getType() == 2.0)
      numGenes = 3;                  //3 variables for Einasto


  // The genetic algorithm part, and chi2 fitting stuffs
  densProfile parent[ u.getNchrome() ], child[ u.getNchrome() ];
  double        chi2[ u.getNchrome() ];



  // Initialize average arrays for diff threads
  double avgChiThreads[ u.getNthreads() ];
  for (int i = 0; i  <  u.getNthreads()  ; ++i)
    avgChiThreads[i]=0;


  // Set all parent R_maxes to the lens
  for (int i = 0; i  <  u.getNchrome() ; ++i){
    parent[i].setR_max( halo.getRmax() );
     child[i].setR_max( halo.getRmax() );
  }

  // Generate first parent values
  #pragma omp parallel for
  for ( int i=0; i<u.getNchrome(); ++i ){

    double gAnalyticArray[ u.getNbins()   ];

    // Generate initial parent values, if Ein profile need alpha parameter
    if ( profile.getType() == 2 )
    parent[i].setAlpha(          randVal( u.getAlphaMin(), u.getAlphaMax() )   );
    parent[i].setM_enc( pow( 10, randVal( u.getMassMin (), u.getMassMax () ) ) );
    parent[i].setC    (          randVal( u.getConMin  (), u.getConMax  () )   );



    // Gives analytic value for RTS, for source locations, radially averaged
    if ( profile.getType() == 1 ){
      generateNFWRTS( gAnalyticArray, parent[i], u.getNbins(), dArr, cosmo.SigmaCrit( halo.getZ(), u.getSourceZ() ) ); // Note this distance is physical seperation
    } else{
//      generateEinRTS( gAnalyticArray, parent[i], halo, u, sourceSigC, sourceDist);
    }

    // Calc chi2 between theoretical predictions and real values
    chi2[i] = chiSquared( gAnalyticArray, gArr, gErrArr, u.getNbins() );

    avgChiThreads[omp_get_thread_num()] += chi2[i];
  }

  // Avg used for reproduction, randomly sample halos based on goodness / 1.5avg
  for (int i  = 0; i < u.getNthreads(); ++i){
    avgChi           += avgChiThreads[i];
    avgChiThreads[i]  = 0;
  }
  avgChi  = avgChi / u.getNchrome();



  // Reproduction time
  int loopCounter(0); // Counter of number of iterations, over 1e6 then stop
  int     counter(0); // Counter used to count number of times within tolerance
  double   totAvg(0); // Keeps running average, stop when converges

  double convChiArray[ u.getNtrack() ];


  // Reset chi^2 values
  // Tot average is a running average over recent runs, eventually should converge
  //  to a fairly constant value
  // Old average saves the previous tot, to measure how much it is changing

  do {
    oldAvg  =   totAvg;
    totAvg  =  0;

    // If tracking array not full, need to only loop over number we have
    //  Otherwise, we only keep track of N_chiTrack
    int N_cT= std::min( u.getNtrack(), loopCounter );

    for (int i=0; i < N_cT; ++i)         // Determines current totAvg over recent generations
    totAvg += convChiArray[i];

    totAvg  = totAvg / N_cT;


    testVal =   avgChi * u.getTestVal(); // Some value of goodness over chi, like 1.5*chi
    avgChi  =        0;                  // Chi over a generation

    // Stores the child's chi2, before transfers to their children
    double newChi2[ u.getNchrome() ];

    // Generate new children
    #pragma omp parallel for
    for ( int i=0; i < u.getNchrome(); ++i ){

      double gAnalyticArray[ u.getNbins()   ];

      int parentIndex[2];

      // Select mates from parent pool
      for ( int j=0; j < 2; ++j ){
        do{
          parentIndex[j] = round( randVal( 0, u.getNchrome()-1) ); // Randomly select an index

        } while ( chi2[parentIndex[j]] > randVal( 0, testVal) );   // If parent fit good enough (small) select w/ uniform distribution
      }



      int aIndex, cIndex, mIndex;

      // Get indexes to use from parents, since parents are random doesn't
      //   matter which we select from, but if Einasto need to decide
      //   whether alpha coming from parent 1 or 2
      cIndex = parentIndex[0];
      mIndex = parentIndex[1];

      // If Einasto, another parameter needed
      if ( profile.getType() == 2 ){
      aIndex = parentIndex[ (int) round( randVal( 0, 1) )  ];


      child[i].setAlpha( parent[aIndex].getAlpha() );
      }

      // Set child values from parents
      child[i].setC    ( parent[cIndex].getC()     );
      child[i].setM_enc( parent[mIndex].getM_enc() );

      // Mutate values, with random chance. Don't exceed possible max or min
      // Needed to keep genes fresh

      double upperRand = 1.1;
      double lowerRand = 0.9;

      if ( u.getMutChance() > randVal(0,1) )
        child[i].setC    ( std::max( std::min( child[i].getC    () * randVal( lowerRand, upperRand ) ,        u.getConMax()   ) ,        u.getConMin()   ) );

      if ( u.getMutChance() > randVal(0,1) )
        child[i].setM_enc( std::max( std::min( child[i].getM_enc() * randVal( lowerRand, upperRand ) , pow(10,u.getMassMax()) ) , pow(10,u.getMassMin() )) );

      if ( u.getMutChance() > randVal(0,1) && profile.getType() == 2 )
        child[i].setAlpha( std::max( std::min( child[i].getAlpha() * randVal( lowerRand, upperRand ) ,        u.getAlphaMax() ) ,        u.getAlphaMin() ) );


      // Generate analytic RTS for child
      if ( profile.getType() == 1 ){
        generateNFWRTS( gAnalyticArray, child[i], u.getNbins(), dArr, cosmo.SigmaCrit( halo.getZ(), u.getSourceZ() ) ); // Note this distance is physical seperation

      } else{
//        generateEinRTS( gAnalyticArray, child[i], halo, u, sourceSigC, sourceDist);
      }


      // chi2 for the children
      newChi2[i] = chiSquared( gAnalyticArray, gArr, gErrArr, u.getNbins() );


      avgChiThreads[ omp_get_thread_num() ] += chi2[i];
    }

    // Avg used for reproduction, randomly sample halos based on goodness / avg
    for (int i  = 0; i < u.getNthreads(); ++i){

      avgChi           += avgChiThreads[i];
      avgChiThreads[i]  = 0;

    }
    avgChi  =  avgChi / u.getNchrome();


    // Tracks the most recent N_ChiTrack generation's chi, to test
    //  whether converged to constant value
    convChiArray[ loopCounter % u.getNtrack() ] = avgChi;

    // Children replace parents
    for ( int i=0; i<u.getNchrome(); ++i ){
      parent[i] =   child[i];
        chi2[i] = newChi2[i];
    }

    // Need counter to be greater than u.consistent,
    //  a count of how many times average has been
    //  consistently below tolerance. Tests for convergence
    if ( fabs(totAvg-oldAvg)/oldAvg < u.getTolerance() ){
      counter += 1;
    }
    else{
      counter  = 0;
    }
if ( loopCounter % 10 == 0 )
printf("%7i    %14.4e %14.4e %14.4e   %3i\n",loopCounter, totAvg, oldAvg, fabs(totAvg-oldAvg)/oldAvg, counter);
    ++loopCounter;
  } while ( u.getNConsistent() > counter && loopCounter < u.getMaxFitNum() ) ;


  // Once converged, just find best child

  int    minIndex =              0; //index of lowest chi2
  double minChi   = chi2[minIndex];
  double cAvg(0), mAvg(0), weightedChi(0);

  //Find lowest chi2 index
  for ( int i = 0; i < u.getNchrome(); ++i ){
    if ( minChi > chi2[i] ){
      minIndex = i;
      minChi   = chi2[i];
    }
//         avgArr[iBin] += avgVal / n_srcs / ( errors[i]*errors[i] ); // Weighted average sum

    weightedChi += 1/chi2[i] ;

    cAvg +=             child[i].getC()       /  chi2[i]  ;
    mAvg += std::log10( child[i].getM_enc() ) /  chi2[i]  ;
  }

  weightedChi = 1 / weightedChi;
  cAvg = cAvg * weightedChi ;
  mAvg = mAvg * weightedChi ;

printf("\n");
printf("%14.4e %7.5f %14.4e\n", weightedChi, cAvg, pow(10, mAvg ) );
printf("%14.4e ", minChi);




  profile = child[minIndex];

}





//
// Generates the radially averaged reduced tangential shear for
//  a NFW profile for given input
//
void generateNFWRTS(
                          double          *gArr ,  // RTS array to output
                    const densProfile     &lens ,  // Input density profile to generate profile for
                    const double         N_bins ,  // Actual information from the halo
                    const double          *dist ,  // Projected distances between source and lens
                    const double           SigC ){ // Critical surface density of sources


  // Loop over all distances, determining predicted rts for a given dist
  for ( int i = 0; i < N_bins ; ++i ){

    double    SD =    SDNFW( dist[i], lens ); // At radius
    double avgSD = SDAvgNFW( dist[i], lens ); // Average

    gArr[i] = ( avgSD - SD ) / ( SigC - SD );
  }
}



// Surface density at input radius for NFW profile, integrated to R_max
double    SDNFW( const double               r ,  //Input radius to calc SD at
                 const densProfile inpProfile ){ //Input NFW profile

  double      x = fabs( r / inpProfile.getR_s() );
  double      c =           inpProfile.getC  ()  ;
  double factor =       2 * inpProfile.getR_s() * inpProfile.getRho_o();

  // Outside of our integration radius
  if ( x > inpProfile.getC() ){
    return 0;
  }
  // If within r_s, arctanh solution has imaginary component we ignore
  else if ( x < 1 ){
    return factor/( x*x - 1 ) * ( sqrt( c*c - x*x ) / ( c + 1 ) - 0.5 / sqrt( 1 - x*x ) * ( log( ( 1. / c * sqrt( ( c*c - x*x )/( 1. - x*x )   ) + 1 ) /
                                                                                                 ( 1. / c * sqrt( ( c*c - x*x )/( 1. - x*x )   ) - 1 ) ) -
                                                                                            log( (          sqrt( ( c*c - x*x )/( 1. - x*x )   ) + 1 ) /
                                                                                                 (          sqrt( ( c*c - x*x )/( 1. - x*x )   ) - 1 ) )   ) ); //+i pi to the logs
  }
  // Outside of r_s
  else if ( x > 1 ) {
    return factor/( x*x - 1 ) * ( sqrt( c*c - x*x ) / ( c + 1)  + 1.0 / sqrt( x*x - 1 ) * ( atan(  1. / c * sqrt( ( c*c - x*x )/( x*x - 1. )         ) ) -
                                                                                            atan(           sqrt( ( c*c - x*x )/( x*x - 1. )         ) )   ) );
  }
  // If close enough to r_s
  else {
    return factor/3.*pow( c*c - 1., -1.5 ) * ( c * ( c*c - 1) - 2 * c + 2 );
  }

}


//Surface density at input radius for NFW profile, integrated to R_max
double    SDAvgNFW( const double               r ,  //Input radius to calc SD at
                    const densProfile inpProfile ){ //Input NFW profile

  double      x = fabs( r / inpProfile.getR_s() );
  double      c =           inpProfile.getC  ()  ;
  double factor =       4 * inpProfile.getR_s() * inpProfile.getRho_o() / ( x*x );



  // If very close to our max integration radius C, errors can occur.
  // Therefore, need explicit solution for C
  if ( fabs(  x - inpProfile.getC    () ) < 1e-4){
    return factor * ( SDAvgNFW( inpProfile.getR_s  (), inpProfile ) / ( 4 * inpProfile.getR_s  () * inpProfile.getRho_o() )
                  - 2.0 * sqrt( c*c - 1 ) / ( c + 1 )
                  + 0.5 * log( ( 1  + 1   /   c * sqrt( c*c - 1 ) )
                        /      ( 1  - 1   /   c * sqrt( c*c - 1 ) ) ) );
  }

  // Outside of C, integral from 0 to 1 + integral from 1 to C, and take out their factors to use the current one
  else if ( x > inpProfile.getC() ){
    return factor * ( SDAvgNFW( inpProfile.getR_s  (), inpProfile ) / ( 4 * inpProfile.getR_s  () * inpProfile.getRho_o() )
                    + SDAvgNFW( inpProfile.getR_max(), inpProfile ) / ( 4 * inpProfile.getR_max() * inpProfile.getRho_o() / ( inpProfile.getC() * inpProfile.getC() ) ) );
  }

  // If within r_s, this is analytic solution integrating to a max value of r_s
  else if ( x < 1 ){

    return factor * ( ( sqrt( c*c - x*x )  - c )/( c + 1 )

                    +         log( c + 1 )

                    - atanh( sqrt( 1 - x*x/( c*c )) )

                    + 0.5 /  sqrt( 1 - x*x ) * ( log( ( 1. / c * sqrt( ( c*c - x*x )/( 1. - x*x )   ) + 1 )
                                             /        ( 1. / c * sqrt( ( c*c - x*x )/( 1. - x*x )   ) - 1 ) )
                                             -   log( (          sqrt( ( c*c - x*x )/( 1. - x*x )   ) + 1 )
                                             /        (          sqrt( ( c*c - x*x )/( 1. - x*x )   ) - 1 ) )   ) );
  }

  // Outside of r_s but less than c, integral from 0 to 1 + new component to a max of rmax/rs
  else if ( x > 1 ){
    return factor * (  SDAvgNFW( inpProfile.getR_s(), inpProfile  ) / ( 4 * inpProfile.getR_s() * inpProfile.getRho_o() )
                    + ( ( sqrt( c*c - x*x ) - 2 * sqrt( c*c - 1 ) ) / ( c + 1 )

                    + 0.5 * log( ( 1. / c * sqrt( c*c -  1  ) + 1 ) / ( 1 - 1. / c * sqrt( c*c -  1  ) ) )
                    - 0.5 * log( ( 1. / c * sqrt( c*c - x*x ) + 1 ) / ( 1 - 1. / c * sqrt( c*c - x*x ) ) )

                    - 1.0 / sqrt( x*x - 1 ) * ( atan(  1. / c * sqrt( ( c*c - x*x )/( x*x - 1. )         ) )
                                            -   atan(           sqrt( ( c*c - x*x )/( x*x - 1. )         ) )   ) ));
  }

  // If close enough to r_s will integrate to 1, hyperbolic previous arctan component turns into sqrt component
  else {
    return factor * ( ( 2 * sqrt( c*c - 1 ) - c ) / ( c + 1 )

                  +         log( c + 1 )

                  - atanh( sqrt( 1 - x*x / ( c*c ) ) )  );
  }

}


// Average surface density within radius for NFW
double SDAvgNFWFull(
                const double     r ,  // Input radius
                const double   r_s ,  // Scale radius, r_-2
                const double rho_o ){ // Initial density

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

// Surface density at radius for NFW profile integrated to infinity along LOS
double    SDNFWFull(
                const double   r   ,  // Input radius
                const double   r_s ,  // Scale radius, r_-2
                const double rho_o ){ // Initial density

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







//   2 0[              (0,1)     ]    1     ^ G(2/as)G(s-1/2)           -2s         (-1)^k G(1/a-2/a k)  2k        a(-1)^k G(-1/2-a/2 k)  ak+1
// H    [                        ]= -----  |  ------------------------ X   ds=  Sum ------ ------------ X    + Sum ------- ------------- X
//   1 2[(0,2/a)(-1/2,1)         ]  2pi i U            G( s )                   k=0 (k  )! G(1/2-    k)        k=1 2(k)!   G(    -a/2 k)
double foxH2012(
                double         x ,  // Z from fox H function
                double     alpha ,  // Shape parameter Ein profile
                double tolerance ){ // Tolerance level for convergence

  // Sum3 only for specific conditions, when second order poles
  double sum1(0), sum2(0), oldSum1(0), oldSum2(0), s1(0), s2(0);
  double sum3(0);

  int converge(0), k(0);

  bool secondOrder = false;

  //k = 0 term
  sum1 =      tgamma(   1. / alpha  ) / tgamma( 0.5 );
//  sum3 = 1. / tgamma( - 1. / alpha+1) * ( x*x / 4 ) ;//* ( - log( x / 2.) + diGamma(1.) + 1./alpha * diGamma(-1./alpha)-diGamma(-1.) );

printf("%4s %4s %14s %14s %14s %14s\n", "k", "conv", "sum1", "sum2", "sum3", "Gammas");
  do {
      ++k;

      int            sign = pow(   -1,   k );
      double       x2kPow = pow(    x, 2*k );
      double       x1kPow = pow(    x, 1+k*alpha );
      double invFactorial = 1./tgamma( 1+k );

      // Terms in summation
      s1   = sign * x2kPow * invFactorial *
           tgamma(  ( 1. - 2.    * k ) / alpha      ) / tgamma(  0.5 - k );


      s2   = sign * x1kPow * invFactorial *
           tgamma( -( 1. + alpha * k ) / 2.         ) / tgamma(      - k * alpha / 2.0);

      sum3+= sign * x2kPow * invFactorial * invFactorial / pow( 2, 2*k ) *
           tgamma( 2 * k + 1 ) / tgamma( ( 2 * k - 1 )   / alpha + 1 ) *
           ( -log( x / 2.) - 1 / (2.0*k) + diGamma( k+1. ) + diGamma( (2.*k-1.)/alpha ) / alpha - diGamma( 2*k-1. ) );

printf("%4i %4i %14.7e %14.7e %14.7e   %14.7e %14.7e \n", k, converge, (double) sum1, (double) sum2, (double) sum3, tgamma((1.-2.*k)/alpha)/tgamma(0.5-k), tgamma(-(1.+alpha*k)/2.0)/tgamma(-k*alpha/2.0));

      //Don't include inf or NaN, which only occurs for second order poles
      if ( ! (std::isinf( s1 ) || s1!=s1) ) {
        oldSum1  = sum1;
           sum1 += s1;
      } else {
        secondOrder = true;
      }

      if ( ! (std::isinf( s2 ) || s2!=s2) ) {
        oldSum2  = sum2;
           sum2 += s2;
      } else {
        secondOrder = true;
      }

      //totSum is term for current k
      if ( fabs((sum1-oldSum1)/oldSum1) < tolerance &&
           fabs((sum2-oldSum2)/oldSum2) < tolerance ){
        converge +=1;
      }
      else{
        converge = 0;
      }
  } while ( converge < 10 && k < 1e1 );

  if ( !secondOrder ) sum3 = 0; // If we are not including secondOrder, set sum to 0

  sum2 *= alpha/2.;
  sum3 *= alpha/tgamma( 0.5 );
printf("%14.7e   %14.7e + %14.7e + %14.7e =  %14.7e\n\n      ", x, (double) sum1, (double) sum2, (double) sum3, sum1+sum2+sum3);
//std::cout<<sum3<<std::endl;;
  return sum1+sum2+sum3;
}



//  2 1[     (-1/2,1) (0,1)     ]    1     ^ G(2/as)G(s-1/2) G(3/2-s)  -2s        (-1)^k G(1/a-2/a k)  2k+2      a(-1)^k G(-3/2-a/2 k)  ak+3
// H   [                        ]= -----  |  ------------------------ X   ds= Sum ------ ------------ X    - Sum ------- ------------- X
//  2 3[(0,2/a)(-1/2,1)(-3/2,1) ]  2pi i U   G(5/2-s) G( s )                  k=0 (k+1)! G(1/2-1/a  )        k=1 2(k)!   G(    -a/2 k)
double foxH2123(
                double         x ,  // Z from fox H function
                double     alpha ,  // Shape parameter Ein profile
                double tolerance ){ // Tolerance level for convergence

  double sum1(0), sum2(0), oldSum1(0), oldSum2(0), s1(0), s2(0);

  int converge(0), k(0);

  //k = 0 term
     sum1 += tgamma( 1. / alpha ) / tgamma( 0.5 ) * ( x * x );
//  oldSum  += sum1;
//printf("%4s %4s %14s %14s\n", "k", "conv", "sum1", "sum2");
  do {
      ++k;

      int            sign = pow(   -1,     k         );
      double       x2kPow = pow(    x, 2 * k + 2     );
      double       x3kPow = pow(    x, 3 + k * alpha );
      double invFactorial = 1./tgamma( 1 + k         );
      double invFactorial1= 1./tgamma( 2 + k         );


      //Terms in summation
      s1 = sign * x2kPow * invFactorial1* //   pow(       - 1     , k ) / factorial(k+1) *    pow( z, 2 + k * 2     ) *      // (-1)^k/(k+1)! x^(2k + 2)
           tgamma(  ( 1. - 2.    * k ) / alpha        ) / tgamma(  0.5 - k              );       // G(  (1-2k)/a ) / G( 1/2 - 1/a )

      s2 = sign * x3kPow * invFactorial * //   pow(       - 1     , k ) / factorial(k  ) *    pow( z, 3 + k * alpha ) *      // (-1)^k/(k  )! x^(ak + 3)
           tgamma( -( 3. + alpha * k ) / 2.           ) / tgamma(      - k * alpha / 2.0);  // G( -(3+ak)/2 ) / G( -ak/2 )

//printf("%4i %4i %14.7e %14.7e    %14.7e %14.7e \n", k, converge, (double) sum1, (double) sum2, tgamma((1.-2.*k)/alpha)/tgamma(0.5-k), tgamma(-(3.+alpha*k)/2.)/tgamma(-k*alpha / 2.0));


      //Don't include inf or NaN
      if ( ! (std::isinf( s1 ) || s1!=s1) ) {
        oldSum1  = sum1;
           sum1 += s1;
      }

      if ( ! (std::isinf( s2 ) || s2!=s2) ) {
        oldSum2  = sum2;
           sum2 += s2;
      }

      //totSum is term for current k
      if ( fabs((sum1-oldSum1)/oldSum1) < tolerance &&
           fabs((sum2-oldSum2)/oldSum2) < tolerance ){
        converge +=1;
      }
      else{
        converge = 0;
      }

  } while ( converge < 10 && k < 1e2 );

  sum2 *= alpha/2.;
//printf("%14.7e   %14.7e + %14.7e =  %14.7e\n      ", x, (double) sum1, (double) sum2, sum1-sum2);
//exit(0);
  return sum1-sum2;
}


void generateEinRTS(
                    double          *gArr ,  // RTS array to output
                    densProfile     &lens ,  // Input density profile to generate profile for
                    userInfo            u ,  // Information from the user
                    double     sourceSc   ,  // Critical surface density for the sources
                    double    *sourceDist ){ // Projected distances between source and lens


  for (int i=0;i<u.getNbins();++i)
    gArr[i] = 0;


  //Constant part of kappa_c, divided by gamma(1/alpha) * sqrt pi, a constant for easier kappas
  //Modified kappa_c, kappa_c * sqrt(pi) / Gamma(1/alpha), just need to multiply by H function

  double kappa_c = lens.getRho_o() * lens.getR_s  ()   * exp( 2./lens.getAlpha() )       *
                         pow(        lens.getAlpha()/2.,      1./lens.getAlpha() - 1.0 ) *
                         tgamma( 1./ lens.getAlpha() ) / sourceSc;
//printf("%7s %14s     %14s %14s      %14s %14s\n", "dist", "gArr", "k_EIN", "k_NFW", "k_AVG_EIN", "k_AVG_NFW");
lens.setAlpha( 1.0 );
  for ( int i = 0; i < u.getNbins(); ++i )
  {
    double x        = sourceDist[i] / lens.getR_s();

//kappa(x)=kc*x*k1(x)=kc*x*G2002=kc*k*2.55912

    double kappa    ;//= kappa_c * std::sqrt( M_PI ) / tgamma( 1. / lens.getAlpha() ) * foxH2012( x, lens.getAlpha() );
    double kappaAVG ;//= kappa_c * std::sqrt( M_PI ) / tgamma( 1. / lens.getAlpha() ) * foxH2123( x, lens.getAlpha() );
kappa    = foxH2012( x, lens.getAlpha() )*std::sqrt(M_PI)/x;
kappaAVG = foxH2123( x, lens.getAlpha() )*x*kappa_c/2;
    // RTS avg across bin
    gArr[i] = ( kappaAVG - kappa ) / ( 1 - kappa );
printf("%7s %14s     %14s %14s      %14s %14s\n", "dist", "R_s", "k_EIN", "k_NFW", "k_AVG_EIN", "k_AVG_NFW");
printf("%7.2f %14.7e     %14.7e %14.7e      %14.7e %14.7e\n",sourceDist[i],x,
                                                             kappa   , SDNFWFull   ( sourceDist[i], lens.getR_s(), lens.getRho_o())/sourceSc,
                                                             kappaAVG, SDAvgNFWFull( sourceDist[i], lens.getR_s(), lens.getRho_o())/sourceSc);
exit(0);
//printf("%12.3e %12.3e %12.3e\n", kappa/kappa_cM, kappaAVG/kappa_cM, gArr[i]);
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





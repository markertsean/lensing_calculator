#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <omp.h>
#include <vector>
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


  double smallStep = 500;  // Number of possible steps in a direction
  double   bigStep = 250;

  double decrement = 1/10.; // Amount to decrease step size by

  double boundDist = 1e-1; // Percent distance can be close to edge

  // Step sizes partially random
//  double rStep = ( u.getRMaxFit () - u.getRMinFit () ) / randVal( bigStep, smallStep ) ;
  double cStep = ( u.getConMax  () - u.getConMin  () ) / randVal( bigStep, smallStep ) ;
  double mStep = ( u.getMassMax () - u.getMassMin () ) / randVal( bigStep, smallStep ) ;
  double aStep = ( u.getAlphaMax() - u.getAlphaMin() ) / randVal( bigStep, smallStep ) ;

  // For kicking the ball back
//  double rStep_0 = rStep;
  double cStep_0 = cStep;
  double mStep_0 = mStep;
  double aStep_0 = aStep;


  // Previous values, to test for convergence
  double prevM = 1;
  double prevC = 1;
  double prevA = 1;
//  double prevR = 1;


  int loopCounter = 0;  // Number of times we have rolled
  int convCounter = 0;  // Number of times we have been in tolerance


  // Roll the ball
  do{

      // Used to check which direction is best
      densProfile testProfile;

      if ( ball.getType() == 2 ) testProfile.setType( 2 );


      double mVals[3], cVals[3], aVals[3];

      // Sets the parameters for moving along the hill
      for ( int j = -1; j < 2; ++j ){
        if ( ball.getType() == 2 )
        aVals[j+1] =          std::min(std::max(        ball.getAlpha()   + aStep * j, u.getAlphaMin() ), u.getAlphaMax()   );
        cVals[j+1] =          std::min(std::max(        ball.getC    ()   + cStep * j, u.getConMin  () ), u.getConMax  ()   );
        mVals[j+1] = pow( 10, std::min(std::max( log10( ball.getM_enc() ) + mStep * j, u.getMassMin () ), u.getMassMax () ) );
      }


      double chiChoose(999), compChi(998);

      double gAnalyticArray[  u.getNbins()   ];

      int maxK = 1;
      if ( ball.getType() == 2 )
        maxK = 3;

      testProfile.setR_max( ball.getR_max() );

      for ( int ii = 0; ii <    3; ++ii ){
                                          testProfile.setC    ( cVals[ii] );
      for ( int jj = 0; jj <    3; ++jj ){
                                          testProfile.setM_enc( mVals[jj] );
      for ( int kk = 0; kk < maxK; ++kk ){


        // Generate RTS
        if ( ball.getType() == 2 ){
          testProfile.setAlpha( aVals[kk] );

          generateEinRTS     ( gAnalyticArray, testProfile, u,            dArr, sigmaC );
        } else
        if ( ball.getType() == 0 ){
          generateNFWTruncRTS( gAnalyticArray, testProfile, u.getNbins(), dArr, sigmaC );
        } else{
          generateNFWRTS     ( gAnalyticArray, testProfile, u.getNbins(), dArr, sigmaC );
        }


        // Determine goodness of fit
        compChi = chiSquared( gAnalyticArray, gArr, gErrArr, u.getNbins() );


        // Checks each possible direction for ball to roll
        // Will roll in direction of best fit
        if (compChi < chiChoose){
          chiChoose = compChi;
          if ( ball.getType() == 2 )
          ball.setAlpha( aVals[kk] );
          ball.setC    ( cVals[ii] );
          ball.setM_enc( mVals[jj] );


          chi2 = compChi;
        }

      }
      }
      }

      if ( fabs( ball.getC    () - u.getConMax () ) / u.getConMax () < boundDist ||
           fabs( ball.getC    () - u.getConMin () ) / u.getConMin () < boundDist ){
           ball.setC(     randVal( u.getConMin ()   , u.getConMax () ) );
           cStep = cStep_0;
      }

      if ( fabs( log10(ball.getM_enc()) - u.getMassMax() ) / u.getMassMax() < boundDist*boundDist ||
           fabs( log10(ball.getM_enc()) - u.getMassMin() ) / u.getMassMin() < boundDist*boundDist ){
           ball.setM_enc( pow( 10, randVal( u.getMassMin (), u.getMassMax () ) ) );
           mStep = mStep_0;
      }

      if ( ball.getType() == 2 ) {
        if(
           fabs( ball.getAlpha() - u.getAlphaMax () ) / u.getAlphaMax () < boundDist ||
           fabs( ball.getAlpha() - u.getAlphaMin () ) / u.getAlphaMin () < boundDist ){
           ball.setAlpha( randVal( u.getAlphaMin ()   , u.getAlphaMax () ) );
           aStep = aStep_0;
        }
      }

      // If didn't move in one direction, decrease step size
      if ( ball.getType () !=       2  &&
           ball.getC    () == cVals[1] &&
           ball.getM_enc() == mVals[1] ){
        cStep *= decrement;
        mStep *= decrement;
      }else
      if (   ball.getType () ==       2  &&
             ball.getC    () == cVals[1] &&
             ball.getM_enc() == mVals[1] ){
        if ( ball.getAlpha() == aVals[1] ){
          cStep *= decrement;
          mStep *= decrement;
          aStep *= decrement;
        }
      }



      // Need counter to be greater than u.consistent,
      //  a count of how many times average has been
      //  consistently below tolerance. Tests for convergence
      if (  fabs( ball.getM_enc() - prevM ) / prevM  < u.getTolerance() &&
            fabs( ball.getC    () - prevC ) / prevC  < u.getTolerance() ){

        if ( ball.getType() != 2 ){
          convCounter += 1;
        } else
        if ( fabs( ball.getAlpha() - prevA ) / prevA  < u.getTolerance() ){  // Either alpha converged, or not einasto
          convCounter += 1;
        }
      } else {
        convCounter  = 0;
      }


      if ( ball.getType() == 2 )
      prevA = ball.getAlpha();
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
                    const userInfo               u ,  // Info from the user
                    const double       *      gArr ,  // RTS binned array we "observed"
                    const double       *      dArr ,  // Distance binned array
                    const double       *   gErrArr ,
                    const COSMOLOGY          cosmo ){ // Error array in RTS


  densProfile ball[ u.getNchrome() ];
  double      chi2[ u.getNchrome() ];

  for (int i = 0; i  <  u.getNchrome() ; ++i){
    if ( profile.getType() == 2 )
    ball[i].setType (                  2 );
    ball[i].setR_max( profile.getR_max() );
    chi2[i] = 1e4;
  }

  // Do for each rolling ball
  #pragma omp parallel for
  for ( int i = 0; i < u.getNchrome(); ++i ){


    // Set starting parameters, location on hill
    if ( profile.getType() == 2 )
    ball[i].setAlpha(          randVal( u.getAlphaMin() , u.getAlphaMax() )   );
    ball[i].setC    (          randVal( u.getConMin  () , u.getConMax  () )   );
    ball[i].setM_enc( pow( 10, randVal( u.getMassMin () , u.getMassMax () ) ) );


    rollBall( ball[i], chi2[i], gArr, dArr, gErrArr, cosmo.SigmaCrit( halo.getZ(), u.getSourceZ() ), u );
  }

  // Now take weighted average of the balls, weighted by their chi2


  int    minIndex =              0; //index of lowest chi2
  double minChi   = chi2[minIndex];

  //Find lowest chi2 index
  for ( int i = 0; i < u.getNchrome(); ++i ){
    if ( minChi > chi2[i] ){
      minIndex = i;
      minChi   = chi2[i];
    }

  }

  if ( profile.getType() == 2 )
  profile.setAlpha(  ball[minIndex].getAlpha()  );
  profile.setC    (  ball[minIndex].getC    ()  );
  profile.setM_enc(  ball[minIndex].getM_enc()  );
  profile.setR_max(  ball[minIndex].getR_max()  );


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
                    const double       *      gArr ,  // RTS binned array we "observed"
                    const double       *      dArr ,  // Distance binned array
                    const double       *   gErrArr ,
                    const COSMOLOGY          cosmo ){ // Error array in RTS

  // Percent closeness to edge we can be
  double edgeTolerance = 1e-2;


  // Values used for comparing goodness of fits
  double  avgChi(0), oldAvg(0), testVal(0);


  // The genetic algorithm part, and chi2 fitting stuffs
  densProfile parent[ u.getNchrome() ], child[ u.getNchrome() ];
  double        chi2[ u.getNchrome() ];



  // Initialize average arrays for diff threads
  double avgChiThreads[ u.getNthreads() ];
  for (int i = 0; i  <  u.getNthreads()  ; ++i)
    avgChiThreads[i]=0;


  // Set all parent R_maxes to the lens
  for (int i = 0; i  <  u.getNchrome() ; ++i){

    if ( profile.getType() == 2 ){

      parent[i].setType ( 2 );
       child[i].setType ( 2 );


    }

//    parent[i].setR_max( halo.getRmax() );
//     child[i].setR_max( halo.getRmax() );
  }

  double sigmaC = cosmo.SigmaCrit( halo.getZ(), u.getSourceZ() );

  // Generate first parent values
//  #pragma omp parallel for
  for ( int i=0; i<u.getNchrome(); ++i ){

    double gAnalyticArray[ u.getNbins()   ];

    // Generate initial parent values, if Ein profile need alpha parameter
    if ( profile.getType() == 2 )
    parent[i].setAlpha(          randVal( u.getAlphaMin(), u.getAlphaMax() )   );
    parent[i].setC    (          randVal( u.getConMin  (), u.getConMax  () )   );
    parent[i].setR_max(          randVal( u.getRMinFit (), u.getRMaxFit () )   );
    parent[i].setM_enc( pow( 10, randVal( u.getMassMin (), u.getMassMax () ) ) );



    // Gives analytic value for RTS, for source locations, radially averaged
    if ( profile.getType() == 1 ){
      generateNFWRTS( gAnalyticArray, parent[i], u.getNbins(), dArr, sigmaC ); // Note this distance is physical seperation
    } else{
      generateEinRTS( gAnalyticArray, parent[i], u           , dArr, sigmaC );
    }

    // Calc chi2 between theoretical predictions and real values
    chi2[i] = chiSquared( gAnalyticArray, gArr, gErrArr, u.getNbins() );

    if ( chi2[i] != chi2[i]  ||
         isinf(     chi2[i]) ){
      chi2[i] = 1e5;
    }
//printf("%14.4e    %10.6f %14.4e %10.6f \n",chi2[i],parent[i].getC(),parent[i].getM_enc(), parent[i].getR_max());

    avgChiThreads[omp_get_thread_num()] += std::log10(chi2[i]);
  }



  // Avg used for reproduction, randomly sample halos based on goodness / 1.5avg
  for (int i  = 0; i < u.getNthreads(); ++i){
    avgChi           += avgChiThreads[i];
    avgChiThreads[i]  = 0;
  }
  avgChi  = pow( 10, avgChi / u.getNchrome());

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
    int N_invalid(0);

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

      int selectCounter(0);

      // Select mates from parent pool
      for ( int j=0; j < 2; ++j ){
        do{
                       parentIndex[j] = round( randVal( 0, u.getNchrome()-1) ); // Randomly select an index

          if ( ++selectCounter > u.getNchrome() ) {                     // If having trouble finding good values, increase range
            avgChi = avgChi * 1.1;
            break;
          }

        } while ( chi2[parentIndex[j]] >       randVal( 0, testVal ) ); // If parent fit good enough (small) select w/ uniform distribution
      }



      int aIndex, cIndex, mIndex, rIndex;

      // Get indexes to use from parents, since parents are random doesn't
      //   matter which we select from, but if Einasto need to decide
      //   whether alpha coming from parent 1 or 2
      cIndex = parentIndex[0];
      mIndex = parentIndex[1];
      rIndex = parentIndex[ (int) round( randVal( 0, 1) )  ];

      // If Einasto, another parameter needed
      if ( profile.getType() == 2 ){
      aIndex = parentIndex[ (int) round( randVal( 0, 1) )  ];

      child[i].setAlpha( parent[aIndex].getAlpha() );
      }

      // Set child values from parents
      child[i].setC    ( parent[cIndex].getC()     );
      child[i].setR_max( parent[rIndex].getR_max() );
      child[i].setM_enc( parent[mIndex].getM_enc() );




      // Mutate values, with random chance. Don't exceed possible max or min
      // Needed to keep genes fresh

      double upperRand = 1.05;
      double lowerRand = 0.95;

      if ( u.getMutChance() > randVal(0.0,1.0) )
        child[i].setC    ( std::max( std::min( child[i].getC    () * randVal( lowerRand, upperRand ) , u.getConMax()   * (1-edgeTolerance) ) ,
                                                                                                       u.getConMin()   * (1+edgeTolerance) ) );


      if ( u.getMutChance() > randVal(0.0,1.0) )
        child[i].setR_max( std::max( std::min( child[i].getR_max() * randVal( lowerRand, upperRand ) , u.getRMaxFit()  * (1-edgeTolerance) ) ,
                                                                                                       u.getRMinFit()  * (1+edgeTolerance) ) );


      if ( u.getMutChance() > randVal(0.0,1.0) && profile.getType() == 2 )
        child[i].setAlpha( std::max( std::min( child[i].getAlpha() * randVal( lowerRand, upperRand ) , u.getAlphaMax() * (1-edgeTolerance) ) ,
                                                                                                       u.getAlphaMin() * (1+edgeTolerance) ) );


      if ( u.getMutChance() > randVal(0.0,1.0) )
        child[i].setM_enc( std::max( std::min( child[i].getM_enc() * randVal( lowerRand, upperRand ) , pow(10,u.getMassMax() * (1-edgeTolerance) )) ,
                                                                                                       pow(10,u.getMassMin() * (1-edgeTolerance) )) );


      // Generate analytic RTS for child
      if ( profile.getType() != 2 ){
        generateNFWRTS( gAnalyticArray, child[i], u.getNbins(), dArr, sigmaC ); // Note this distance is physical seperation
      } else{
        generateEinRTS( gAnalyticArray, child[i], u           , dArr, sigmaC ); // Note this distance is physical seperation
      }

      newChi2[i] = chiSquared( gAnalyticArray, gArr, gErrArr, u.getNbins() );

      // Avoid counting NaNs
      if ( newChi2[i] == newChi2[i] ){
        avgChiThreads[ omp_get_thread_num() ] += newChi2[i];
      } else {
        ++N_invalid;
      }

    }

// Make one for outside function
    // Find the best fits
    double minChis[ u.getNtrack() ];

    for ( int i = 0; i < u.getNtrack(); ++i )
      minChis[i] = 1e9 * (i+1);



    // Find the top fits
    for ( int i = 0; i < u.getNchrome(); ++i ){         // Check each child
      if ( chi2[i] < minChis[ u.getNtrack()-1 ] ){      // If child better than at least the worst we track

        int startIndex = u.getNtrack()-1;               // Know it's at least as good as last index
        for ( int j = u.getNtrack()-2; j > -1 ; --j ){  // Once we find better chi, know it's place
          if ( chi2[i] > minChis[j] ){                  // If we check against worse, use previous index
            break;
          }
            startIndex = j;
        }

        // Replace the values
        double  replaceVal = chi2[i];
        double replacedVal;

        for ( int j = startIndex; j < u.getNtrack(); ++j ){

          replacedVal = minChis[j];
          minChis[j]  = replaceVal;
          replaceVal  = replacedVal;

        }
      }
    }


    // Avg used for reproduction, randomly sample halos based on goodness / avg
    for (int i  = 0; i < u.getNthreads(); ++i){

      avgChi           += avgChiThreads[i];
      avgChiThreads[i]  = 0;

    }
    avgChi  =  avgChi / u.getNthreads() ;


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

/*

    oldAvg  =   totAvg;
    totAvg  =  0;

    // If tracking array not full, need to only loop over number we have
    //  Otherwise, we only keep track of N_chiTrack
    int N_cT= std::min( u.getNtrack(), loopCounter );
    int N_invalid(0);

    for (int i=0; i < N_cT; ++i)         // Determines current totAvg over recent generations
    totAvg += convChiArray[i];

    totAvg  = totAvg / N_cT;


    testVal =   avgChi * u.getTestVal(); // Some value of goodness over chi, like 1.5*chi
    avgChi  =        0;                  // Chi over a generation

//*/



    ++loopCounter;

  } while ( u.getNConsistent() > counter && loopCounter < u.getMaxFitNum() ) ;


  // Once converged, just find best child

  int    minIndex =              0; //index of lowest chi2
  double minChi   = chi2[minIndex];
  double cAvg(0), mAvg(0), aAvg(0), rAvg(0), weightedChi(0);

  //Find lowest chi2 index
  for ( int i = 0; i < u.getNchrome(); ++i ){
    if ( chi2[i] == chi2[i] &&
         chi2[i]  >      0  ){
    if ( minChi > chi2[i] ){
      minIndex = i;
      minChi   = chi2[i];
    }
//         avgArr[iBin] += avgVal / n_srcs / ( errors[i]*errors[i] ); // Weighted average sum

    weightedChi += 1/chi2[i] ;

    if ( profile.getType() == 2 )
    aAvg +=             child[i].getAlpha()   /  chi2[i]  ;

    cAvg +=             child[i].getC()       /  chi2[i]  ;
    rAvg +=             child[i].getR_max()   /  chi2[i]  ;
    mAvg += std::log10( child[i].getM_enc() ) /  chi2[i]  ;
    }
  }

  weightedChi =    1 / weightedChi ;
  if ( profile.getType() == 2 )
  aAvg        = aAvg * weightedChi ;
  rAvg        = rAvg * weightedChi ;
  cAvg        = cAvg * weightedChi ;
  mAvg        = mAvg * weightedChi ;

printf("\n");
printf("%14.4e %7.5f %14.4e ", weightedChi, cAvg, pow(10, mAvg ) );
if ( profile.getType() == 2 ){
printf("%7.5f ", aAvg);}
else{
  printf("        ");}
printf("%7.5f\n", rAvg );
printf("%14.4e ", minChi);




  profile = child[minIndex];

}



//
// Generates the radially averaged reduced tangential shear for
//  a NFW profile for given input
//
void generateNFWTruncRTS(
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

    double    SD =    SDNFWFull( dist[i], lens.getR_s(), lens.getRho_o() ); // At radius
    double avgSD = SDAvgNFWFull( dist[i], lens.getR_s(), lens.getRho_o() ); // Average

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



void generateEinRTS(
                          double          *gArr ,  // RTS array to output
                    const densProfile     &lens ,  // Input density profile to generate profile for
                    const userInfo            u ,  // Information from the user
                    const double    *sourceDist ,  // Projected distances between source and lens
                    const double     sourceSc   ){ // Critical surface density for the sources


  for (int i=0;i<u.getNbins();++i)
    gArr[i] = 0;


  //Constant part of kappa_c, divided by gamma(1/alpha) * sqrt pi, a constant for easier kappas
  //Modified kappa_c, kappa_c * sqrt(pi) / Gamma(1/alpha), just need to multiply by H function

  double kappa_c = lens.getRho_o() * lens.getR_s  ()   * exp( 2./lens.getAlpha() )       *
                         pow(        lens.getAlpha()/2.,      1./lens.getAlpha() - 1.0 ) *
                         tgamma( 1./ lens.getAlpha() ) / sourceSc;


  double modKappa_c = kappa_c / std::sqrt( M_PI ) / tgamma( 1. / lens.getAlpha() );

  for ( int i = 0; i < u.getNbins(); ++i ){

    double        x = sourceDist[i]  /  lens.getR_s();

    double kappa    = modKappa_c     * pow( 10, interpolateEinRTS( x, lens.getAlpha(), einKappa    ) ); // Interpolate table of Kappa    values
    double kappaAvg = modKappa_c * x * pow( 10, interpolateEinRTS( x, lens.getAlpha(), einKappaAvg ) ); // Interpolate table of KappaAvg values

    gArr[i] = ( kappaAvg - kappa ) / ( 1 - kappa );

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


double interpolateEinRTS(  double        x ,  // r/r_s
                           double        a ,  // alpha
                           einTable  table ){ // Table to interpolate on

  double x_step = table.getX_bins() / ( table.getX_max() - table.getX_min() );  // Bin / value conversion
  double a_step = table.getA_bins() / ( table.getA_max() - table.getA_min() );

  double  x_bin = ( x - table.getX_min() ) * x_step ; // Decimal bin positions
  double  a_bin = ( a - table.getA_min() ) * a_step ;

  int    x1_bin = (int) floor( x_bin );               // Bin locations of the points
  int    a1_bin = (int) floor( a_bin );


  int    x2_bin = (int) ceil ( x_bin );
  int    a2_bin = (int) ceil ( a_bin );

  double x1     = x1_bin / x_step + table.getX_min(); // Actual position in x and a coordinates
  double x2     = x2_bin / x_step + table.getX_min();

  double a1     = a1_bin / a_step + table.getA_min();
  double a2     = a2_bin / a_step + table.getA_min();

  return 1.0 / ( ( x2 - x1 ) * ( a2 - a1 ) ) * (
               ( ( x2 - x  ) * ( a2 - a  )   * table.getVal( a1_bin, x1_bin ) ) +
               ( ( x  - x1 ) * ( a2 - a  )   * table.getVal( a1_bin, x2_bin ) ) +
               ( ( x2 - x  ) * ( a  - a1 )   * table.getVal( a2_bin, x1_bin ) ) +
               ( ( x  - x1 ) * ( a  - a1 )   * table.getVal( a2_bin, x2_bin ) ) );
}


/*



bool num_den(  int       &p ,  // Finds ratio of p/q, if rational
               int       &q ,
               double alpha ,
               int     maxK ){


    int k = 0;
    int q1(0), q2(0);

    bool  secondOrder = false; // Whether rational
    bool first_s1_nan = true;  // The NaNs are used to find ratio p/q = n = 1/alpha
    bool   sec_s1_nan = true;
    bool first_s2_nan = true;

    do{
      ++k;

      double testVal1 = tgamma( (1.0 - 2. * k) / alpha ) ;
      double testVal2 ;

      bool nan_s1 = !( testVal1 == testVal1 ); // If s1 is NaN
      if ( first_s1_nan && nan_s1 ){
        q1 = k;
        first_s1_nan = false;
      } else
      if (   sec_s1_nan && nan_s1 && !first_s1_nan ){
        q2 = k;
        sec_s1_nan = false;
      }

      testVal1 = tgamma( -( 1. + alpha * k ) / 2. );
      testVal2 = tgamma(       - alpha * k   / 2. );

      bool nan_s2 = !( testVal1 == testVal1 && // If s1 is NaN
                       testVal2 == testVal2 );
      if ( first_s2_nan && nan_s2 ){
        p = k;
        first_s2_nan = false;
      }

      if ( p!=0  && q1!=0 &&  q2!=0 ) {
        secondOrder = true;
        break;
      }
    } while ( k < maxK ); // Abort if won't be second order in our integration range

    q = q2-q1;

  if( q == 0 || p == 0 ){

    q=0;
    p=0;
    secondOrder = false;
  }

  return secondOrder;
}



//   2 0[              (0,1)     ]    1     ^ G(2/as)G(s-1/2)           -2s         (-1)^k G(1/a-2/a k)  2k        a(-1)^k G(-1/2-a/2 k)  ak+1
// H    [                        ]= -----  |  ------------------------ X   ds=  Sum ------ ------------ X    + Sum ------- ------------- X
//   1 2[(0,2/a)(-1/2,1)         ]  2pi i U            G( s )                   k=0 (k  )! G(1/2-    k)        k=1 2(k)!   G(    -a/2 k)
void   foxH2012(
                double  f2012Arr[]     ,  // Values to return
                double  f2123Arr[]     ,  // Values to return
                double         x[]     ,  // Z from fox H function
                double    N_bins       ,  // Number of bins of z
                double     alpha       ,  // Shape parameter Ein profile
                double tolerance       ){ // Tolence level for convergence



  // Sum3 only for specific conditions, when second order poles

  int maxK = 5e2;           // Maximum number of terms in sum
  int maxC = 30 ;           // Number of times to converge

  double maxG = 1e-300;     // Max ratio of Gamma/Gamma to accept, otherwise drop


  bool secondOrder = false; // For all second order terms
  int p(0), q(0);           // Nonzero for second order


  // Sets ratio p/q = n = 1/alpha, both 0 if nothing found
  if ( num_den( p, q, alpha, 10*maxK ) ) {
    secondOrder = true;
  } else{
    p = q = 0;
  }

  int k0 = ( q - 1 ) / 2;          // For checking sums

  short sConv[6];                  // Whether to keep summing each individual sum

  short converge [ (int)N_bins ] ; // Count number of times converged
  short converge2[ (int)N_bins ] ; // Count number of times converged
  short firstIndex = 0           ; // So don't loop over converged indexes of sum1-3
  short seconIndex = 0           ; // So don't loop over converged indexes of sum4-6


  mpf_class   alph( alpha );       // Keep high precision
  mpf_class     km( 0     );
  mpf_class   sign( 1     );       // -1 or 1
  mpf_class s3Sign( 1     );

  if ( p % 2 != 0 ) s3Sign = -1;   // Negative if p is odd

  mpf_class     s1( 0     );       // Individual terms in sums
  mpf_class     s2( 0     );
  mpf_class     s3( 0     );


  mpf_class  invJunk ( 0 );         // Sum3 inverse x term
  mpf_class  invJunk2( 0 );         // Sum6 inverse x term
  mpf_class  diG1    ( 0 );         // Sum3 gamma terms
  mpf_class  diG2    ( 0 );
  mpf_class  diG3    ( 0 );
  mpf_class  diG4    ( 0 );
  mpf_class  s       ( 0 );

  mpf_class invFactorial( 1 );     // Inverse of the factorial


  typedef std::vector<mpf_class> mpf_v;

  // Sum terms
  mpf_v sum1(N_bins), oldSum1( N_bins ); // Oldsum for convergence
  mpf_v sum2(N_bins), oldSum2( N_bins );
  mpf_v sum3(N_bins), oldSum3( N_bins );
  mpf_v sum4(N_bins), oldSum4( N_bins );
  mpf_v sum5(N_bins), oldSum5( N_bins );
  mpf_v sum6(N_bins), oldSum6( N_bins );


  mpf_v      z ( N_bins );               // High precision x
  mpf_v      z2( N_bins );
  mpf_v x2kPow ( N_bins );               // Powers, x^something, with the inverse factorial term
  mpf_v x1kPow ( N_bins );
  mpf_v x4kPow ( N_bins );
  mpf_v xakPow ( N_bins );
  mpf_v logJunk( N_bins );

  // k = 0 term
  sum1[0] = gamma( 1 / alph ) / gamma( mpf_class( 0.5 ) );

  // Initialize array elements
  for ( int i  = 0; i < N_bins; ++i ){
    sum1[i] = sum1[0];
    sum2[i] =      0 ;
    sum3[i] =      0 ;
    sum4[i] = sum1[0];
    sum5[i] =      0 ;
    sum6[i] =      0 ;

    oldSum1[i] = sum1[0];
    oldSum2[i] =      0 ;
    oldSum3[i] =      0 + tolerance/2;
    oldSum4[i] = sum1[0];
    oldSum5[i] =      0 ;
    oldSum6[i] =      0 + tolerance/2;

    z [i] = x[i] * pow( 2./alpha, 1./alpha );
    z2[i] = z[i] * z[i];

    x1kPow[i] = z[i];
    x2kPow[i] =   1;
    x4kPow[i] =   1 ;
    xakPow[i] =  pow( z[i], alpha );

    logJunk[i] = -ln( z[i] / 2 );

    converge[i] = 0;

  }

  sConv[0] = 0;
  sConv[1] = 0;
  sConv[2] = 0;
  sConv[3] = 0;
  sConv[4] = 0;
  sConv[5] = 0;

  mpf_class rootPi( gamma( mpf_class(0.5) ) ); // Used a lot

  // Sum over km
  do {

      km    = km + 1;

      int k = (int) km.get_d();

      sign  = sign * -1         ; // -1^k

      // Raise x's to the k's
      for ( int i = 0; i < N_bins; ++i ){
        x4kPow[i] = x4kPow[i] / km / km * z2[i] / 4  ;  // (x/2)^2k    /k!/k!
        x2kPow[i] = x2kPow[i] / km      * z2[i]      ;  //  x   ^2k    /k!
        x1kPow[i] = x1kPow[i] / km      * xakPow[i]  ;  //  x   ^1+ak  /k!
      }

      double testVal1 = tgamma( (1.0 - 2. * k) / alpha ) ; // Tests whether NaNs
      double testVal2 = tgamma(  0.5 -      k          ) ;

      // First sum
      if ( (     !secondOrder     ||   // If second order, need special condition
           ( (k+k0) % q != 0 ) )  &&   // Still need to check for testVal values

           ( testVal1 == testVal1 &&   // Make sure no NaNs
             testVal2 == testVal2)&&

           ((sConv[0] == 0       )||   // Sum1 or 4 has not converged
            (sConv[3] == 0        ))
                                  ){

            mpf_class g1( spouges(  ( 1 - 2 * km ) / alph ) ); // Gamma function to 50 digits of precision
            mpf_class g2( spouges(  0.5 -     km ) );

            mpf_class g( g1/g2 );



            if ( abs( g ) < maxG  ) {                     // If out of our range, stop summing
              sConv[0] = 1;
              sConv[3] = 1;
            }

            s1   = sign * g;                              // Term in the sum

            for ( int i = firstIndex; i < N_bins; ++i ){  // Only take sums that haven't converged
              sum1[i] += s1 * x2kPow[i];
            }

            for ( int i = seconIndex; i < N_bins; ++i ){  // Only take sums that haven't converged
              sum4[i] += s1 * x2kPow[i] / ( km + 1 );
            }

      }


      testVal1 = tgamma( -( 1. + alpha * k ) / 2. );
      testVal2 = tgamma(       - alpha * k   / 2. );

      // Sum 2
      if ( (  !secondOrder        ||
             ( k % p != 0 )     ) &&

           ( testVal1 == testVal1 &&
             testVal2 == testVal2)&&

           ((sConv[1] == 0       )||   // Sum2 or 5 has not converged
            (sConv[4] == 0        ))
                                  ){

            mpf_class g1( spouges( -( 1 + alph * km ) / 2 ) ); // Gamma functions
            mpf_class g2( spouges( - km * alph / 2  ) );

            mpf_class g( g1/g2 );

            if ( abs( g ) < maxG  ) // If out of our range, stop summing
              sConv[1] = 1;


            s2  = sign * g;      // Term in sum

        for ( int i = firstIndex; i < N_bins; ++i ){
          sum2[i] += s2 * x1kPow[i];
        }
        for ( int i = seconIndex; i < N_bins; ++i ){
          sum5[i] += s2 * x1kPow[i] / ( - 1.5 - alph * km / 2 );
        }
      }


      testVal1 = tgamma(   2 * km.get_d() + 1 );
      testVal2 = tgamma( ( 2 * km.get_d() - 1 ) / alph.get_d() + 1 );

      if (        secondOrder &&  // If second order, need extra sum
           ( (k+k0) % q == 0) &&

           (( sConv[2]  == 0) ||
            ( sConv[5]  == 0))
                              ){

        invJunk  =         - 1 / ( 2 * km     ) ;
        invJunk2 = invJunk - 1 / ( 2 * km + 2 ) ;

        s =       km + 1.;
        diG1 = diGamma( s );

        s = (2. * km - 1.)/alph;
        diG2 = diGamma( s );

        s =  2. * km - 1.;
        diG3 = diGamma( s );


        diG4 = diG1 + 1 / ( km + 1 ) ; // phi( z+1 ) = phi( z ) + 1/z


        if ( abs( mpf_class(1) - alph ) <= 1e-8 ){
          diG2 = 0;
          diG3 = 0;
        }

        mpf_class g1( spouges(   2 * km + 1 ) );
        mpf_class g2( spouges( ( 2 * km - 1 ) / alph + 1 ) );

        mpf_class g( g1/g2 );

        if ( abs( g ) < maxG  ) // If out of our range, stop summing
           sConv[2] = 1;


        s3   = s3Sign * g;

        for ( int i   = firstIndex; i < N_bins; ++i ){
             sum3[i] += s3 * x4kPow[i] *
             ( logJunk[i] +
               invJunk    +
                  diG1    +
                  diG2    / alph -
                  diG3    );
        }

        for ( int i   = seconIndex; i < N_bins; ++i ){
             sum6[i] += s3 * x4kPow[i] / ( km + 1 ) *
             ( logJunk[i] +
               invJunk2   +
                  diG4    +
                  diG2    / alph -
                  diG3    );
        }


     }

      // Check each for convergence
      // Find what the total sum will be, compare it to previous
      // If converged maxC times, can stop
      for ( int i = firstIndex; i < N_bins; ++i ){

        mpf_class    sum(    sum1[i] + alph/2 *    sum2[i] ); // Add the two sums
        mpf_class oldSum( oldSum1[i] + alph/2 * oldSum2[i] );


        if ( secondOrder ){
             sum =    sum + alph/rootPi *    sum3[i];         // If necessary add third
          oldSum = oldSum + alph/rootPi * oldSum3[i];
        }


        mpf_class diff( abs( ( sum - oldSum ) / sum ) ) ;     // How much it has changes

        if ( diff < tolerance ){
          converge[i] += 1;
        } else {
          converge[i]  = 0;
        }

        if ( converge [i] >= maxC       && // If first sum  we checked converged,
                       i  == firstIndex ){ //    don't need to check anymore
            firstIndex += 1;
        }

        oldSum1[i] = sum1[i]; // For next iteration, this counts NaNs as converging
        oldSum2[i] = sum2[i];
        oldSum3[i] = sum3[i];
      }


      for ( int i = seconIndex; i < N_bins; ++i ){

        mpf_class    sum(    sum4[i] - alph/2 *    sum5[i] ); // Add the two sums
        mpf_class oldSum( oldSum4[i] - alph/2 * oldSum5[i] );


        if ( secondOrder ){
             sum =    sum + alph / rootPi *    sum6[i];       // If necessary add third
          oldSum = oldSum + alph / rootPi * oldSum6[i];
        }

           sum =    sum / z[i];
        oldSum = oldSum / z[i];


        mpf_class diff( abs( ( sum - oldSum ) / sum ) ) ;     // How much it has changes

        if ( diff < tolerance ){
          converge2[i] += 1;
        } else {
          converge2[i]  = 0;
        }

        if ( converge2[i] >= maxC       && // If first sum  we checked converged,
                       i  == seconIndex ){ //    don't need to check anymore
            seconIndex += 1;
        }

        oldSum4[i] = sum4[i]; // For next iteration, this counts NaNs as converging
        oldSum5[i] = sum5[i];
        oldSum6[i] = sum6[i];
      }


//printf("\n");
  } while ( firstIndex != N_bins &&  // firstIndex will increase when first one converges
                    km  < maxK   );


  mpf_class returnValues( 0 );


  // The final sum, sum3 is 0 if not second order
  for ( int i = 0; i < N_bins; ++i ){

    returnValues = sum1[i] + alph/mpf_class(2) * sum2[i] + alph / rootPi * sum3[i];

    f2012Arr[i]  = returnValues.get_d();

    returnValues = sum4[i] - alph/mpf_class(2) * sum5[i] + alph / rootPi * sum6[i];

    f2123Arr[i]  = returnValues.get_d() / x[i];

  }

}


// Done in other function
//  2 1[     (-1/2,1) (0,1)     ]    1     ^ G(2/as)G(s-1/2) G(3/2-s)  -2s        (-1)^k G(1/a-2/a k)  2k+2      a(-1)^k G(-3/2-a/2 k)  ak+3
// H   [                        ]= -----  |  ------------------------ X   ds= Sum ------ ------------ X    - Sum ------- ------------- X
//  2 3[(0,2/a)(-1/2,1)(-3/2,1) ]  2pi i U   G(5/2-s) G( s )                  k=0 (k+1)! G(1/2-1/a  )        k=1 2(k)!   G(    -a/2 k)
double foxH2123(
                double    retArr[]     ,  // Values to return
                double         z[]     ,  // Z from fox H function
                double    N_bins       ,  // Number of bins of z
                double     alpha       ,  // Shape parameter Ein profile
                double tolerance       ){ // Tolerance level for convergence
}


//*/

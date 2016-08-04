#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <iostream>


//Returns value of chi^2 fit, ignores points w/0 error
double chiSquared(
                   const double    *array1 ,   // First  array to compare
                   const double    *array2 ,   // Second array to compare
                   const double  *errorArr ,   // Error values in measurement
                   const int    N_elements ){  // Number of elements in arrays, if different size might break

  //Chi^2 = sum ( predicted - data / error )^2
  double chi2(0), temp(0);


  for ( int i=0; i<N_elements; ++i ){

    if ( errorArr[i] != 0.0 ){             //If no error, ignore

      temp  = ( array1[i] - array2[i] ) / errorArr[i] ;
      chi2 += temp * temp;
    }

  }


  return chi2;
}

//Generates uniformly distribution between low and high
double randVal(
                double  low ,  //Lower limit of random number
                double high ){ //Upper limit of random number

  return ((double) rand() / RAND_MAX) * (high - low) + low;

}


//Like randVal, but with integers
double randVal(
                int     low ,  //Lower limit of random number
                double high ){ //Upper limit of random number

  return randVal( (double) low, high);
}


//Like randVal, but with integers
double randVal(
                double  low ,  //Lower limit of random number
                int    high ){ //Upper limit of random number

  return randVal( low, (double) high);
}



//Like randVal, but with integers. This version returns an integer.
int      randVal(
                int     low ,  //Lower limit of random number
                int    high ){ //Upper limit of random number

  return (int) (randVal( (double) low, (double) high));
}

//Simple factorial function
int     factorial(
                  int  facVal){ //Value to calculate factorial of

  int product(1);

  for (int i=facVal; i>0; --i)
    product *= i;

  return product;
}

/*
// Gamma'(z) / Gamma(z)
long double diGamma(  long double  z ){


  long double  gamma = 0.57721566490153286060651209008240243104215933593992359880576L; // Euler-Mascheroni constant

  // Integral solution, but also follows relation phi(x+1) = phi(x) + 1/x

  if ( z >= 1.L &&
       z <= 3.L ){   // We can integrate

    long double stepSize = 1.0e-7L;
    long double      s   = z-1.L;
    long double      sum = 0.L;

    // Midpoint Integration, offset x by half a step
    for ( long double x = stepSize/2.0; x < 1.0  ; x+= stepSize ){

      sum += stepSize *              // dx
             ( 1. - pow( x, s ) ) /  // ( 1-x^s )
             ( 1. -      x      ) ;  // ( 1-x)

    }

    return sum - gamma;

  } else
  if ( z > 3.L ){
    return diGamma( z - 1 ) + 1.L/(z-1.L);
  } else {
    return diGamma( z + 1 ) - 1.L/z;
  }
}

//*/

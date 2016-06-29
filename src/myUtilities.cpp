#include <cstdlib>
#include <cmath>
#include <stdio.h>


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


// Gamma'(z) / Gamma(z)
double diGamma(  double  z        ,
                 double tolerance ,
                 int    N_consis  ){

  double    sum = 0;
  double oldSum = 1000;

  int    consis = 0;
  int         n = 0;

  double  gamma = 0.577215664901532860606512090082402; // Euler-Mascheroni constant

  do {

    ++n;

    sum += 1./( n + z ) - 1./n;

    if ( fabs((sum-oldSum)/oldSum) < tolerance ) {
      ++consis;
    } else {
      consis = 0;
    }
    oldSum = sum;

  } while ( consis < N_consis && n < 1e6);

  return - ( 1./z + gamma + sum );

}

#include <cstdlib>


//Returns value of chi^2 fit, ignores points w/0 error
double chiSquared(
                    double    *array1 ,   // First  array to compare
                    double    *array2 ,   // Second array to compare
                    double  *errorArr ,   // Error values in measurement
                    int    N_elements ){  // Number of elements in arrays, if different size might break

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

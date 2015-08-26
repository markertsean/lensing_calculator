#include <cstdlib>

//Returns value of chi^2 fit, ignores points w/0 error
double chiSquared( double *array1, double *array2, double *errorArr, int N_elements){
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
double randVal( double low, double high ){
  return ((double) rand() / RAND_MAX) * (high - low) + low;
}

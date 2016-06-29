#ifndef MY_UTILITIES
#define MY_UTILITIES



double chiSquared(
		  const double    *array1,  //One set of values to compare
		  const double    *array2,  //Set to compare against
		  const double  *errorArr,  //Errors to divide by
		  const int    N_elements); //Size of the arrays


double    randVal(
		  double  low,  //Lower limit on randVal
		  double high); //Upper limit on randVal

double    randVal(
		  int     low,  //Lower limit on randVal
		  double high); //Upper limit on randVal

double    randVal(
		  double  low,  //Lower limit on randVal
		  int    high); //Upper limit on randVal

int       randVal(
		  int     low,  //Lower limit on randVal
		  int    high); //Upper limit on randVal

int     factorial(
                  int  facVal); //Value to calculate factorial of

double diGamma(  double  z               ,
                 double tolerance = 1e-9,
                 int    N_consis  =  20  );


#endif

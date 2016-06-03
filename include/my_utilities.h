#ifndef MY_UTILITIES
#define MY_UTILITIES



double chiSquared(
		  double    *array1,  //One set of values to compare
		  double    *array2,  //Set to compare against
		  double  *errorArr,  //Errors to divide by
		  int    N_elements); //Size of the arrays


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


#endif

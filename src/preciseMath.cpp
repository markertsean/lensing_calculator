#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>

#include "precise_math.h"
#include "mpcomplex.h"
#include "matrix.h"


extern "C" {
  #include "bernoulli.h"
}



using namespace std;


mpdouble Comb(int n, int k)
{
	mpdouble r;

	if (k < 0) return r;
	r = 1;
	for (int l = 1; l <= n-k; l++)
	{
		r *= n + 1 - l;
		r /= l;
	}
	return r;
}

void makeB(matrix<mpdouble> &B, int n)
{
	int r, c;
	for (c = 0; c < n; c++) B[0][c] = 1;
	for (r = 1; r < n; r++)
	{
		for (c = 0; c < n; c++)
		{
			B[r][c] = Comb(r+c-1, c-r);
			if ((r+c)&1) B[r][c] = -B[r][c];
		}
	}
}

void makeC(matrix<mpdouble> &C, int n)
{
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k <= i; k++)
			{
				C[i][j] += Comb(2*i, 2*k) * Comb(k, k+j-i);
			}
			if ((i-j)&1) C[i][j] = -C[i][j];
		}
	}
	C[0][0] = 0.5;
}

void makeD(matrix<mpdouble> &D, int n)
{
	D[0][0] = 1;
	D[1][1] = -1;
	for (int i = 2; i < n; i++)
	{
		D[i][i] = D[i-1][i-1] * (2 * (2*i-1));
		D[i][i] /= (i - 1);
	}
}

void makeF(matrix<mpdouble> &F, mpdouble g, int n)
{
	for (int a = 0; a < n; a++)
	{
		int i;
		F[a][0] = 2;
		for (i = a + 1; i <= 2 * a; i++)
		{
			F[a][0] *= i;
			F[a][0] /= 4;
		}
		F[a][0] *= exp(mpdouble(a) + g + .5);
		F[a][0] /= pow(mpdouble(a) + g + .5, a);
		F[a][0] /= sqrt(mpdouble(a) + g + .5);
	}
}




// Uses Lanczos approx
mpf_class gamma( mpf_class x ){

  double  g( 20.3209821879863739013671875 );
  int     n( 24 );

	mpdouble zRe, zIm;

  zIm = 0;
  zRe = mpdouble (x.get_d());

	matrix<mpdouble>  B(n, n);
	matrix<mpdouble>  C(n, n);
	matrix<mpdouble>  D(n, n);
	matrix<mpdouble>  F(n, 1);
	matrix<mpdouble>  P(n, 1);
	matrix<mpcomplex> Z(1, n);
	matrix<mpcomplex> R(1, 1);

	makeB(B, n);
	makeC(C, n);
	makeD(D, n);
	makeF(F, g, n);
	P = D*B*C*F;

	bool s = false;
  if (zRe < 0)
	{
		zRe = 1 - zRe;
		zIm = -zIm;
		s = true;
	}
	mpcomplex z(zRe - 1, zIm), r;
	Z[0][0] = mpdouble(1);
	for (int i = 1; i < n; i++)
	{
		Z[0][i] = 1/(z+i);
	}
	R = Z * P;
	r = z + g + .5;
	r = log(R[0][0]) + (z+.5) * log(r) - r;
	if (s)
	{
		r = log(mpdouble::pi()) - r - log(-sin(mpdouble::pi() * z));
	}
	//cout << "lnGma(z)=" << r << endl << "Gamma(z)=" << exp(r) << endl;

  std::string myString =  static_cast<ostringstream*>( &(ostringstream() <<  std::real( exp(r) )  ) )->str();

  mpf_class returnVal( myString );

  return returnVal;
}



// Natural logarithm
mpf_class ln( mpf_class inpVal ){

  if ( inpVal <= 0.0 ){
      printf("Error in ln of value: ");
      std::cout << inpVal << std::endl;
      exit(1);
  }

  // use property of ln( a x 10^n ) = ln( a ) + n ln( 10 )
  mpz_class n(      0 );
  mpf_class a( inpVal );

  mpf_class ln10Val(ln_ln10s);
  mpf_class    eVal(ln_es   );

  // Sets a as coefficient, n as exponent
  while ( a >= mpf_class( 10 ) ){
      a = a /  mpf_class( 10 ) ;
      n = n +  mpz_class(  1 ) ;
  }
  while ( a <  mpf_class(  1 ) ){
      a = a *  mpf_class( 10 ) ;
      n = n -  mpz_class(  1 ) ;
  }

  // Now have ln( a ) + n ln( 10 ), 1 < a < 10, remove any extra e's as well
  // If a > e^2, b = a / e^2, ln( a ) = ln( b ) + 2
  // If a > e^1, b = a / e^1, ln( a ) = ln( b ) + 1
  // If a > e^m, b = a / e^m, ln( a ) = ln( b ) + m

  mpz_class m( 0 );

  // Reusing a is simpler...
  if        ( a > eVal * eVal ) {
    m = m + mpz_class( 2 );
    a = a /     ( eVal * eVal );
  } else if ( a > eVal ){
    m = m + mpz_class( 1 );
    a = a /       eVal;
  }

  // Drop another power of 10 for converging the series
  n = n + mpz_class(  1 );
  a = a / mpf_class( 10 );

  // a should be < 1
  // for x < 1:
  // y = (x-1)/(x+1)
  // ln(x) = 2 * y SUM k = 1 to inf of 1/(2k+1) * y^(2k)
  mpf_class   y( ( a - mpf_class( 1 ) )  /
                 ( a + mpf_class( 1 ) )  );
  mpf_class  yp(                  1      );
  mpf_class sum(                  1      );
  mpz_class   k(                  0      );

  do {
    k   = k   + mpz_class( 1 )              ; // iterate k
    yp  = yp  * y  * y                      ; // = y^(2k)
    sum = sum + yp / mpf_class( 2 * k + 1 ) ;
  } while ( k < 2e2 );

  return mpf_class(2) * y * sum + m + n * ln10Val;
}



// Exponential function
mpf_class exp( mpf_class inpVal ){

// e^8.14, e^178.78, e^-43.5

  // e^( a + b ) = e^a * e^b
  // a whole number
  // b decimal

  mpf_class  eVal(         ln_es     );  // Value of Eulers Constant
  mpz_class  sign(   sgn( inpVal )   );  // Sign of the exponent
  mpf_class     a( trunc( inpVal )   );  // Whole number portion
  mpf_class     b(        inpVal - a );  // Decimal portion
  mpf_class    ea(               1.0 );  // e^a


  // Calc e^a
  for ( int i = 0 ; i < abs( a ); i = i + 1 )
    ea = ea * eVal;

  if ( sign == -1 )
    ea = 1/ea;

  mpz_class   k( 0 ); // Index
  mpf_class sum( 1 ); // Running sum, = e^b
  mpz_class fac( 1 ); // Factorial
  mpf_class  xk( 1 ); // x^k, where x = b

  do {

    k   = k   + 1; // Increment index
    xk  = xk  * b; // Calculates x^k
    fac = fac * k; // Calculates the factoral

    sum = sum + xk / fac;
  } while ( k < 3e1 );


  return ea * sum;
}



// Power
mpf_class pow( mpf_class  a ,
               mpf_class  b ){

  /*
      c =     a^b
   ln c =  ln a^b
        =       b ln a
      c = exp ( b ln a)
  */

  return exp( b * ln( a ) );

}



//need special power function
mpf_class diGamma( mpf_class z ){


  int  lowVal = 17;
  int highVal = lowVal+1;

  // Integral solution, but also follows relation phi(x+1) = phi(x) + 1/x

  if ( z == lowVal ) {

    return mpf_class("2.803513328327460368386716903146");

  } else
  if ( z >   lowVal &&
       z <  highVal ){   // We can integrate

    //mpf_class gammaC( ln_ems ); // Euler-Mascheroni constant

    mpf_class nPhi ( ln( z ) - mpf_class(1) / ( mpz_class(2) * z )  ); // To add to sum

    mpf_class sum ( 0 );
    mpf_class zp  ( 1 );


    // Bernoulli numbers, external function
    mpq_t rop;
    mpz_t n, d;
    mpq_init(rop);
    mpz_inits(n, d, NULL);



    unsigned int kui (0);
    mpf_class    km  (0);

    do{

      kui = kui + 1 ;
      km  = km  + 1 ;
      zp  = zp  * z * z;

      bernoulli( rop, 2 * kui );
      mpq_get_num( n, rop );
      mpq_get_den( d, rop );

      mpq_class  rat( rop );

      mpz_class num( rat.get_num_mpz_t() );
      mpz_class den( rat.get_den_mpz_t() );

      sum = sum + num / ( den * 2 * km * zp );

    } while ( kui < 15 );

    return nPhi-sum;

  } else
  if ( z >= highVal ){
    return diGamma( z - 1 ) + 1.0/( z - 1.0);
  } else {
    return diGamma( z + 1 ) - 1.0/  z;
  }
}


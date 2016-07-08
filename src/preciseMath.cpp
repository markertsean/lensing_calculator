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


  mpf_class  gamma( ln_ems ); // Euler-Mascheroni constant

  // Integral solution, but also follows relation phi(x+1) = phi(x) + 1/x
/*
-0.577215664901532860606512090082402431042159335939923598805767234884867726
1.644934066848226436472415166646025189218949901206798437735558229370007470 pi2/6
-2.404113806319188570799476323022899981529972584680997763584543110683676411
6.493939402266829149096022179247007416648505711512361446097857292664723697 pi4/15
-24.88626612344087823195277167496882003336994206804590748738062426969128615
122.0811674381338967657421515749104633482180988039424274210890396805198619 8pi6/63
-726.0114797149844353246542358918536669119017636069718686204438583315530643
5060.549875237639470468573602083608424905162384875244102458715421369018080 8pi8/15
-40400.97839874763488532782365545085427877962792824191466403941555348253335
363240.91142238262680714352556574776489114436761411760146250647037361744687 128pi10/33
9^
-3.630593311606628712990618842832054105457279034079387429801940443319876172e6
3.992662298773108670232707324047201489779273691976043029806448907811259421e7
-4.790603798898314524268767644990636347189709185641186486362327984306841255e8
6.227402193410971764192853408947415910060510469664605698043745885736784393e9
-8.718095783017206784519122031036435756610645482962514563576904352293937892e10
1.307694352218913820890099907485110170285328633962085017522629158340478485e12
-2.092294967948151090663165568811151437911025141531737318144060238699595562e13
3.556887858592237159756123967161824480929801507384519141393675325474109458e14
-6.402385922818921400735649453323975529947187425748758832484198207540196311e15
1.216452164536393966698766962740413842383112994013078739263026453952777028e17
-2.432903168507861321737256818243197497086463760919206750959866308632857752e18
20^

*/
/*
  mpf_class phi011("1.644934066848226436472415166646025189218949901206798437735558229370007470");
  mpf_class phi021("-2.404113806319188570799476323022899981529972584680997763584543110683676411");
  mpf_class phi031("6.493939402266829149096022179247007416648505711512361446097857292664723697");
  mpf_class phi041("-24.88626612344087823195277167496882003336994206804590748738062426969128615");
  mpf_class phi051("122.0811674381338967657421515749104633482180988039424274210890396805198619");
  mpf_class phi061("-726.0114797149844353246542358918536669119017636069718686204438583315530643");
  mpf_class phi071("5060.549875237639470468573602083608424905162384875244102458715421369018080");
  mpf_class phi081("-40400.97839874763488532782365545085427877962792824191466403941555348253335");
  mpf_class phi091("363240.91142238262680714352556574776489114436761411760146250647037361744687");

  mpf_class phi101("-3.63059331160662871299061884283205410545727903407938742980194044331987617e6");
  mpf_class phi111("3.992662298773108670232707324047201489779273691976043029806448907811259421e7");
  mpf_class phi121("-4.79060379889831452426876764499063634718970918564118648636232798430684125e8");
  mpf_class phi131("6.227402193410971764192853408947415910060510469664605698043745885736784393e9");
  mpf_class phi141("-8.71809578301720678451912203103643575661064548296251456357690435229393789e10");
  mpf_class phi151("1.307694352218913820890099907485110170285328633962085017522629158340478485e12");
  mpf_class phi161("-2.09229496794815109066316556881115143791102514153173731814406023869959556e13");
  mpf_class phi171("3.556887858592237159756123967161824480929801507384519141393675325474109458e14");
  mpf_class phi181("-6.40238592281892140073564945332397552994718742574875883248419820754019631e15");
  mpf_class phi191("1.216452164536393966698766962740413842383112994013078739263026453952777028e17");
  mpf_class phi201("-2.43290316850786132173725681824319749708646376091920675095986630863285775e18");
*/

  mpf_class pPhi( -gamma );

  if ( z >= 1 &&
       z <  2 ){   // We can integrate


    mpf_class stepSize(          1.0e-6  );
    mpf_class      s  (          z-1.0   );
    mpf_class      sum(            0.0   );
    mpf_class      x  ( stepSize / 2.0   );

    // Midpoint Integration, offset x by half a step
    for ( ; x < 1.0  ; x+= stepSize ){

      sum += stepSize *              // dx
             ( 1. - pow( x, s ) ) /  // ( 1-x^s )
             ( 1. -      x      ) ;  // ( 1-x)

    }

/*
  mpf_class eps( z - 1 );


  pPhi = pPhi + phi011 * pow( eps, mpf_class(  1 ) ) +
                phi021 * pow( eps, mpf_class(  2 ) ) / mpz_class( 2 ) +
                phi031 * pow( eps, mpf_class(  3 ) ) / mpz_class( 6 ) +
                phi041 * pow( eps, mpf_class(  4 ) ) / mpz_class( 24 ) +
                phi051 * pow( eps, mpf_class(  5 ) ) / mpz_class( 120 ) +
                phi061 * pow( eps, mpf_class(  6 ) ) / mpz_class( 720 ) +
                phi071 * pow( eps, mpf_class(  7 ) ) / mpz_class( 5040 ) +
                phi081 * pow( eps, mpf_class(  8 ) ) / mpz_class( 40320 ) +
                phi091 * pow( eps, mpf_class(  9 ) ) / mpz_class( 362880 ) +
                phi101 * pow( eps, mpf_class( 10 ) ) / mpz_class( 3628800 ) +
                phi111 * pow( eps, mpf_class( 11 ) ) / mpz_class( 39916800 ) +
                phi121 * pow( eps, mpf_class( 12 ) ) / mpz_class( 479001600 ) +
                phi131 * pow( eps, mpf_class( 13 ) ) / mpz_class( 6227020800 ) +
                phi141 * pow( eps, mpf_class( 14 ) ) / mpz_class( 87178291200 ) +
                phi151 * pow( eps, mpf_class( 15 ) ) / mpz_class( 1307674368000 ) +
                phi161 * pow( eps, mpf_class( 16 ) ) / mpz_class( 20922789888000 ) +
                phi171 * pow( eps, mpf_class( 17 ) ) / mpz_class( 355687428096000 ) +
                phi181 * pow( eps, mpf_class( 18 ) ) / mpz_class( 6402373705728000 ) +
                phi191 * pow( eps, mpf_class( 19 ) ) / mpz_class( 121645100408832000 ) +
                phi201 * pow( eps, mpf_class( 20 ) ) / mpz_class( 2432902008176640000 );
*/
/*
  mpf_class      sum(            0.0   );
  mpf_class nPhi( ln( z ) - mpz_class(1) / ( mpz_class(2) * z )  );

  // Bernoulli numbers
  mpq_t rop;
  mpz_t n, d;
  mpq_init(rop);
  mpz_inits(n, d, NULL);


  sum = 0;

  unsigned int kui (0);
  mpf_class    km  (0);
std::cout.precision(30);
std::scientific;

  do{

    kui = kui + 2;
    km  = km  + 2;
    bernoulli( rop, kui );
    mpq_get_num( n, rop );
    mpq_get_den( d, rop );

    mpq_class rat( rop );

    mpz_class num( rat.get_num_mpz_t() );
    mpz_class den( rat.get_den_mpz_t() );

    sum = sum + num / ( den * km * pow( z, km ) );
*/
/*
std::cout << std::setw(4) <<
kui/2 << " " << std::setw(20) <<
num << " " << std::setw(10) <<
den*km << " " << std::setw(30) <<
pow( z, km ) << " " << std::setw(30) <<
sum << std::endl;

  } while ( kui < 200 );

  nPhi = nPhi + sum;

std::cout.precision(30);
std::scientific;
std::cout<< std::setw(30) << z << std::endl;
//std::cout<< std::setw(30) << eps << std::endl;
std::cout<< std::setw(30) << nPhi << std::endl;
//std::cout<< std::setw(30) << sum-gamma << std::endl;
*/
    return sum - gamma;

  } else
  if ( z > 2 ){
    return diGamma( z - 1 ) + 1.0/( z - 1.0);
  } else {
    return diGamma( z + 1 ) - 1.0/  z;
  }
}

#ifndef PRECISE_MATH

#define PRECISE_MATH

#include <gmpxx.h>

#define ln_ln10   2.302585092994045684017991454684364207601101488628
#define ln_ln10s "2.302585092994045684017991454684364207601101488628"
#define ln_e      2.718281828459045235360287471352662497757247093699
#define ln_es    "2.718281828459045235360287471352662497757247093699"
#define ln_emc    0.57721566490153286060651209008240243104215933593992359880576
#define ln_ems   "0.57721566490153286060651209008240243104215933593992359880576"


mpf_class diGamma( mpf_class z );

//
mpf_class gamma( mpf_class z );


// Roughly 1e2 - 1e4 times slower than exp
mpf_class exp( mpf_class inpVal );

// Hundreds of times slower than pow
mpf_class pow( mpf_class  a , // a^b
               mpf_class  b );


// Roughly 1e5 times slower than log
mpf_class ln( mpf_class inpVal );

mpf_class spouges( mpf_class z );


#endif // PRECISE_MATH

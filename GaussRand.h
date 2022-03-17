#ifndef ISC2019_ARB_GAUSSRAND_H
#define ISC2019_ARB_GAUSSRAND_H

#endif //ISC2019_ARB_GAUSSRAND_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <flint/fmpz_poly.h>
#define PI 3.141592654
void gaussrand(fmpz_poly_t e)
{
    double res_standard;
    int deviation=8;
    int res;
    int len;
    len=fmpz_poly_length(e);
    for(int i=0;i<len;i++)
    {
        res_standard = sqrt(-2.0 * log(rand() / (RAND_MAX + 1.0)))* sin(2.0 * PI * rand() / (RAND_MAX + 1.0));
        res=res_standard*deviation;
        fmpz_poly_set_coeff_si(e,i,res);
    }
}

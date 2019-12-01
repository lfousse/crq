#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

typedef mpz_t *listz_t;

int main (int argc, char *argv[])
{
    unsigned long n;
    int i, st, logn, j, meth;
    mpz_t *coeffs;
    mpz_t d;
    meth = atoi (argv[1]);
    n = atoi (argv[2]);
    coeffs = malloc (((n + 1) / 2) * sizeof (*coeffs));
    for (i = 0; i < (n + 1) / 2; i++)
	mpz_init (coeffs[i]);
    mpz_init (d);

    st = cputime ();
    if (meth == 1)
	compute_bnk (coeffs, d, n - 1);
    else
	compute_bnk2 (coeffs, d, n);
    st = cputime () - st;
    printf ("%d [n]\n%d [time]\n", n, st);
    printf ("[ ");
    for (i = 0; i < (n+1) / 2; i++) {
	mpz_out_str (stdout, 10, coeffs[i]);
	printf (" / ");
	mpz_out_str (stdout, 10, d);
	if (i < (n+1)/2 - 1)
	    printf (", ");
    }
    printf ("];\n");
    return 0;
}

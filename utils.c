/* utils.c -- misc. utility functions.

Copyright 2005 Laurent Fousse <laurent@komite.net>

This file is part of the CRQ Library.

The CRQ Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 2.1
as published by the Free Software Foundation.

The CRQ Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the CRQ Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <sys/time.h>
#include <sys/resource.h>

#include "crq-impl.h"

int cputime (void)
{
    struct rusage rus;

    getrusage (0, &rus);
    return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

void sign_set (char *signmask, int pos, int value) {
    char tmp;
    tmp = ~(1 << pos);
    tmp = *signmask & tmp;
    if (value < 0) {
	tmp |= 1 << pos;
	*signmask = tmp;
    }
    else
	*signmask = tmp;
}

int sign_get (const char signmask, int pos) {
    if (signmask & (1 << pos))
	return -1;
    else
	return 1;
}

#if 1
int sign_eval (mpz_t c, mpz_t *P, unsigned long k, unsigned long n)
{
    mpz_t v;
    mpz_t t;
    unsigned int i;
    int ret;
    mpz_init (v);
    mpz_init (t);

    mpz_set (v, P[n]);
    i = n;
    do {
	i--;
	mpz_mul (v, v, c);
	mpz_mul_2exp (t, P[i], k * (n - i));
	mpz_add (v, v, t);
    } while (i > 0);
    ret = mpz_sgn (v);
    mpz_clear (t);
    mpz_clear (v);
    return ret;
}

#else

void hornerdac (mpz_t r, mpz_t *p, unsigned long n, unsigned long k, mpz_t c)
{
    if (n == 0)
	mpz_set (r, p[0]);
    else {
	/* write P = Pl(x) + x^l*Ph(x):
	 2^(k*(n-1))*P(c/2^k) =
	 2^(k*(h+l-1))*Pl(c/2^k) + 2^(k*(h+l-1))*(c/2^k)^l*Ph(c/2^k)
	 = 2^(k*h) * eval(p, l, k, c) 
	 + c^l     * eval(p + l, h, k, c) */
	unsigned long l = (n+1) / 2; /* Pl has l coeffs */
	unsigned long h = n + 1 - l; /* Ph has h coeffs */
	mpz_t u;
	mpz_init (u);
	hornerdac (r, p + l, h - 1, k, c);
	mpz_pow_ui (u, c, l);
	mpz_mul (r, r, u);
	hornerdac (u, p, l - 1, k, c);
	mpz_mul_2exp (u, u, k * h);
	mpz_add (r, r, u);
	mpz_clear (u);
    }
}

int sign_eval (mpz_t c, mpz_t *p, unsigned long k, unsigned long n)
{
    mpz_t r;
    int ret;
    mpz_init (r);
    hornerdac (r, p, n, k, c);
    ret = mpz_sgn (r);
    mpz_clear (r);
    return ret;
}
#endif

/* Compute sign of P''(c/2^k) */
int sign_eval_p2 (mpz_t c, mpz_t *P, unsigned long k, unsigned long n)
{
    mpz_t v;
    mpz_t t;
    unsigned int i;
    int ret;

    if (n < 2)
	return 0;

    mpz_init (v);
    mpz_init (t);

    mpz_set (v, P[n]);
    mpz_mul_ui (v, v, n);
    mpz_mul_ui (v, v, n - 1);
    i = n;
    do {
	i--;
	mpz_mul (v, v, c);
	mpz_mul_2exp (t, P[i], k * (n - i));
	mpz_mul_ui (t, t, i);
	mpz_mul_ui (t, t, i - 1);
	mpz_add (v, v, t);
    } while (i > 1);
    ret = mpz_sgn (v);
    mpz_clear (t);
    mpz_clear (v);
    return ret;
}

void debug_mpfr_value (mpfr_t foo, char *name)
{
    printf ("%s = ", name);
    mpfr_out_str (stdout, DBASE, 0, foo, GMP_RNDN);
    printf ("\n");
}

void abort_if_nan (mpfr_t foo, char *file, int line)
{
    if (mpfr_nan_p (foo)) {
	printf ("NaN in file %s, line %d.\n", file, line);
	abort ();
    }
}

/* Set ret to ulp(a), to zero if a = 0 */

void mpfr_setulp (mpfr_t ret, mpfr_t a) {
    mp_exp_t exp_a;

    if (mpfr_zero_p (a)) {
	mpfr_set_ui (ret, 0, GMP_RNDN);
	return;
    }
    exp_a = mpfr_get_exp (a);
    mpfr_set_ui_2exp (ret, 1, exp_a - mpfr_get_prec (a), GMP_RNDN);
}

void dump_poly (FILE *out, mpz_t *P, unsigned long d)
{
    unsigned long i;

    for (i = 0; i <= d; i++)
	gmp_fprintf (out, "%Zd\n", P[i]);
}

void load_poly (FILE *in, mpz_t *P, unsigned long d)
{
    unsigned long i;

    for (i = 0; i <= d; i++)
	gmp_fscanf (in, "%Zd\n", P[i]);
}

void dump_interval (FILE *out, interval I)
{
    gmp_fprintf (out, "%Zd %ld\n", I.c, I.k);
}

void load_interval (FILE *in, interval *I)
{
    gmp_fscanf (in, "%Zd %ld\n", I->c, &(I->k));
}

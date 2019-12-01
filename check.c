/* check.c -- check the library via sample integrals.

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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "crq.h"
#include "crq-impl.h"
#include "utils.h"

#define WANT_TIMINGS

extern mpfr_t max_dyi;
void compute_glc (mpfr_t *b, mpfr_t *roots, unsigned long n, mp_prec_t wp,
		  struct crq_gl_opts options, struct crq_stats *timings);
void compute_bnk (mpz_t * b, mpz_t d, unsigned long n);

enum { GL, NC };

void sin_wrapper (mpfr_t res, mpfr_t x)
{
    mpfr_sin (res, x, GMP_RNDN);
}

void square_wrapper (mpfr_t res, mpfr_t x)
{
    mpfr_mul (res, x, x, GMP_RNDN);
}

void sinsin (mpfr_t res, mpfr_t x)
{
    mpfr_sin (res, x, GMP_RNDN);
    mpfr_sin (res, res, GMP_RNDN);
}

void one_wrapper (mpfr_t res, mpfr_t x)
{
    mpfr_set_ui (res, 1, GMP_RNDN);
}

void sin_exp_cos (mpfr_t res, mpfr_t x)
{
    mpfr_t tmp;
    mpfr_init2 (tmp, mpfr_get_prec (res) + 5);
    
    mpfr_sin_cos (res, tmp, x, GMP_RNDN);
    mpfr_exp (tmp, tmp, GMP_RNDN);
    mpfr_mul (res, res, tmp, GMP_RNDN);
    mpfr_clear (tmp);
}

void darctan (mpfr_t res, mpfr_t x)
{
    mpfr_t tmp;
    mpfr_init2 (tmp, mpfr_get_prec (res) + 5);
    mpfr_sqr (tmp, x, GMP_RNDN);
    mpfr_add_ui (tmp, tmp, 1, GMP_RNDN);
    mpfr_ui_div (res, 1, tmp, GMP_RNDN);
    mpfr_clear (tmp);
}

void scaled_darctan (mpfr_t res, mpfr_t x)
{
    mpfr_t tmp;
    mpfr_init2 (tmp, mpfr_get_prec (res) + 5);
    mpfr_mul_ui (tmp, x, 100000, GMP_RNDN);
    mpfr_sqr (tmp, tmp, GMP_RNDN);
    mpfr_add_ui (tmp, tmp, 1, GMP_RNDN);
    mpfr_ui_div (res, 1, tmp, GMP_RNDN);
    mpfr_clear (tmp);
}

void x2sx3 (mpfr_t res, mpfr_t x)
{
    mpfr_t tmp;
    mpfr_init2 (tmp, mpfr_get_prec (res) + 5);
    
    mpfr_pow_ui (tmp, x, 3, GMP_RNDN);
    mpfr_sin (tmp, tmp, GMP_RNDN);
    mpfr_sqr (res, x, GMP_RNDN);
    mpfr_mul (res, res, tmp, GMP_RNDN);
    mpfr_clear (tmp);
}

void expx100 (mpfr_t res, mpfr_t x)
{
    mpfr_pow_ui (res, x, 100, GMP_RNDN);
    mpfr_neg (res, res, GMP_RNDN);
    mpfr_exp (res, res, GMP_RNDN);
}

void sqrtcomp (mpfr_t res, mpfr_t x)
{
    mpfr_sqr (res, x, GMP_RNDN);
    mpfr_ui_sub (res, 1, res, GMP_RNDN);
    mpfr_sqrt (res, res, GMP_RNDN);
}

void sqrt_wrapper (mpfr_t res, mpfr_t x)
{
    mpfr_sqrt (res, x, GMP_RNDN);
}

void max_sin_cos (mpfr_t res, mpfr_t x)
{
    mpfr_t tmp;
    mpfr_init2 (tmp, mpfr_get_prec (res));
    
    mpfr_sin_cos (res, tmp, x, GMP_RNDN);
    if (mpfr_cmp (res, tmp) < 0)
	mpfr_set (res, tmp, GMP_RNDN);
    mpfr_clear (tmp);
}

void expln (mpfr_t res, mpfr_t x)
{
    mpfr_t tmp;
    mpfr_init2 (tmp, mpfr_get_prec (res) + 5);
    mpfr_sqr (tmp, x, GMP_RNDN);
    mpfr_neg (tmp, tmp, GMP_RNDN);
    mpfr_exp (res, tmp, GMP_RNDN);
    mpfr_log (tmp, x, GMP_RNDN);
    mpfr_mul (res, res, tmp, GMP_RNDN);
    mpfr_clear (tmp);
}

void exp_wrapper (mpfr_t res, mpfr_t x)
{
#if VERBOSE >= 3
    printf ("On va calculer exp avec x = ");
    mpfr_out_str (stdout, DBASE, 0, x, GMP_RNDN);
    printf ("\n");
#endif
    mpfr_exp (res, x, GMP_RNDN);
    if (mpfr_nan_p (res))
	printf ("Quelque chose de pourri au royaume du Danemark.\n");
}

void nth_deriv_bound_exp (mpfr_t res, mpfr_t a, mpfr_t b, unsigned int n)
{
    mpfr_exp (res, b, GMP_RNDU);
}

/* t*log(1+t) */
void bailey1 (mpfr_t res, mpfr_t t)
{
    mpfr_log1p (res, t, GMP_RNDN);
    mpfr_mul (res, res, t, GMP_RNDN);
}

/* sqrt(1-t^2) */
void bailey2 (mpfr_t res, mpfr_t t)
{
    mpfr_sqr (res, t, GMP_RNDN);
    mpfr_ui_sub (res, 1, res, GMP_RNDN);
    mpfr_sqrt (res, res, GMP_RNDN);
}

struct gl_params {
    unsigned int n;
    mpfr_t *coeffs;
    mpfr_t *roots;
    mpfr_t **fxs;
};

struct nc_params {
    unsigned int n;
    mpz_t *coeffs;
    mpz_t *d;
    mpfr_t **fxs;
};

void gl_wrapper (mpfr_t res, mpfr_t a, mpfr_t b, aqfunc_t f, void *params)
{
    mpfr_t nan;
    struct gl_params *gparams;
    mpfr_init (nan);
    gparams = (struct gl_params *) params;

    gen_gl_closed (res, nan, a, b, gparams->n, gparams->coeffs, gparams->roots, f, nan,
		   mpfr_get_prec (res), gparams->fxs);
    mpfr_clear (nan);
}

void nc_wrapper (mpfr_t res, mpfr_t a, mpfr_t b, aqfunc_t f, void *params)
{
    mpfr_t nan;
    struct nc_params *ncparams;
    mpfr_init (nan);
    ncparams = (struct nc_params *) params;

    gen_nc_closed (nan, res, a, b, ncparams->n, ncparams->coeffs, *ncparams->d, f, nan, 
		   mpfr_get_prec (res) + mpz_sizeinbase (ncparams->coeffs[0], 2), ncparams->fxs);
    mpfr_clear (nan);
}

int main (int argc, char *argv[])
{
    unsigned int n, m, i;
    mp_prec_t prec, wp = 0;

    mpfr_t a, b;
    mpfr_t res;
    mpfr_t *coeffsgl, *roots, **fxs;
    mpz_t *coeffsnc;
    mpz_t nc_denom;
    mpfr_t exact_value;
    mpfr_t total_error, math_error, roundoff_error;
    mpfr_t nan;
    mpfr_t tmp;
    mpfr_t diffbound, ndiffbound;
    struct gl_params gparams;
    struct nc_params ncparams;
#ifdef WANT_TIMINGS
    int total_time = 0;
#endif
    struct crq_gl_opts options;
    int testcase = -1; /* default */
    int method = GL;
    int adaptive = 0;
    struct crq_stats stats;
    aqfunc_t *f;

    stats.usp_time = 0;
    stats.refine_time = 0;
    stats.coeff_time = 0;
    stats.newton_good = 0;
    stats.newton_bad = 0;
    stats.dich_prep_time = 0;
    stats.dicho_time = 0;
    stats.refinement_failures = 0;
 

    m = 1;
    n = 42;
    mpfr_init (nan);
    mpfr_set_nan (nan);
    options.iso = SIMPLE;
    options.ref = NEWTON;
    options.dicho = FULL;
    options.file = NULL;

    while (argc > 1 && argv[1][0] == '-') {
	if (argc >= 3 && strcmp (argv[1], "-p") == 0) {
	    prec = atoi (argv[2]);	/* target precision */
	    mpfr_set_default_prec (prec);
	    argv += 2;
	    argc -= 2;
	}
	else if (argc >= 3 && strcmp (argv[1], "-wp") == 0) {
	    wp = atoi (argv[2]);	/* working precision */
	    argv += 2;
	    argc -= 2;
	}
	else if (argc >= 3 && strcmp (argv[1], "-test") == 0) {
	    testcase = atoi (argv[2]);	/* which test */
	    argv += 2;
	    argc -= 2;
	}
	else if (argc >= 3 && strcmp (argv[1], "-n") == 0) {
	    n = atoi (argv[2]);	/* Method with n points */
	    argv += 2;
	    argc -= 2;
	}
	else if (argc >= 3 && strcmp (argv[1], "-m") == 0) {
	    m = atoi (argv[2]);	/* number of intervals */
	    argv += 2;
	    argc -= 2;
	}
	else if (argc >= 2 && strcmp (argv[1], "-newton") == 0) {
	    options.ref = NEWTON;
	    argv += 1;
	    argc -= 1;
	}
	else if (argc >= 2 && strcmp (argv[1], "-secante") == 0) {
	    options.ref = SECANTE;
	    argv += 1;
	    argc -= 1;
	}
	else if (argc >= 2 && strcmp (argv[1], "-nc") == 0) {
	    method = NC;
	    argv += 1;
	    argc -= 1;
	}
	else if (argc >= 2 && strcmp (argv[1], "-adapt") == 0) {
	    adaptive = 1;
	    argv += 1;
	    argc -= 1;
	}
	else if (argc >= 2 && strcmp (argv[1], "-gl") == 0) {
	    method = GL;
	    argv += 1;
	    argc -= 1;
	}
	else if (argc >= 2 && strcmp (argv[1], "-newton_int") == 0) {
	    options.ref = NEWTON_INT;
	    argv += 1;
	    argc -= 1;
	}
	else if (argc >= 2 && strcmp (argv[1], "-dicho") == 0) {
	    options.ref = DICHO;
	    argv += 1;
	    argc -= 1;
	}
	else if (argc >= 2 && strcmp (argv[1], "-dicho_trunc") == 0) {
	    options.ref = DICHO;
	    options.dicho = TRUNC;
	    argv += 1;
	    argc -= 1;
	}
	else if (argc >= 2 && strcmp (argv[1], "-usp2") == 0) {
	    options.iso = DOUBLE;
	    argv += 1;
	    argc -= 1;
	}
	else {
	    fprintf (stderr, "Invalid option: %s\n", argv[1]);
	    exit (1);
	}
    }

    if (wp == 0)
	wp = mpfr_get_default_prec ();

    mpfr_init2 (a, wp);
    mpfr_init2 (b, wp);
    mpfr_init2 (res, wp);
    mpfr_init2 (exact_value, wp);
    mpfr_init (math_error);
    mpfr_init2 (diffbound, wp);
    mpfr_init2 (ndiffbound, wp);

    mpfr_set_ui (a, 0, GMP_RNDN);
    mpfr_set_ui (b, 3, GMP_RNDN);

#if VERBOSE >= 1
    printf ("%d [n]\n%d [m]\n%lu [wp]\n%d [test]\n", n, m, wp, testcase);
    if (options.ref == DICHO)
	printf ("[dicho]\n");
    if (options.iso == DOUBLE)
	printf ("[usp_couple]\n");
    if (options.ref == NEWTON_INT)
	printf ("[newton_int]\n");
    if (options.ref == NEWTON)
	printf ("[newton]\n");
#endif
#ifdef WANT_TIMINGS
    total_time = cputime ();
#endif
    if (adaptive == 1) {
	/* We want to use an adaptive method, so compute
	 * weights once and for all */
	fxs = malloc (((n+1)/2) * sizeof (*fxs));
	for (i = 0; i < (n + 1) / 2; i++) {
	    fxs[i] = malloc (sizeof (mpfr_t));
	    mpfr_init2 (*fxs[i], wp + 42);
	}
	if (method == GL) {
	    coeffsgl = malloc (((n + 1) / 2) * sizeof(*coeffsgl));
	    roots = malloc ((n / 2) * sizeof(*coeffsgl));
	    for (i = 0; i < n / 2; i++)
		mpfr_init2 (roots[i], wp);
	    for (i = 0; i < (n+1) / 2; i++)
		mpfr_init2 (coeffsgl[i], wp);
	    compute_glc (coeffsgl, roots, n, wp, options, &stats);
	    gparams.n = n;
	    gparams.roots = roots;
	    gparams.coeffs = coeffsgl;
	    gparams.fxs = fxs;
	} else if (method == NC) {
	    coeffsnc = malloc (((n + 1) / 2) * sizeof (*coeffsnc));
	    for (i = 0; i < (n + 1) / 2; i++)
		mpz_init (coeffsnc[i]);
	    mpz_init (nc_denom);
	    compute_bnk (coeffsnc, nc_denom, n - 1);
	    ncparams.n = n;
	    ncparams.coeffs = coeffsnc;
	    ncparams.fxs = fxs;
	    ncparams.d = &nc_denom;
	}
    }

    switch (testcase) {
	case 0: /* f = 1, [a, b] = [0, 3] */
#if VERBOSE >= 1
	    printf ("f(x) = 1, [a, b] = [0, 3]\n");
#endif
	    f = one_wrapper;
	    mpfr_set_ui (exact_value, 3, GMP_RNDN);
	    mpfr_set_ui (diffbound, 0, GMP_RNDN);
	    mpfr_set_ui (ndiffbound, 0, GMP_RNDN);
	    break;

	case 1: /* f = exp, [a, b] = [0, 3] */
#if VERBOSE >= 1
	    printf ("f(x) = exp(x), [a, b] = [0, 3]\n");
#endif
	    f = exp_wrapper;
	    mpfr_set_ui (diffbound, 3, GMP_RNDN);
	    mpfr_exp (diffbound, diffbound, GMP_RNDU);
	    mpfr_set (ndiffbound, diffbound, GMP_RNDN); /* exact */
	    mpfr_set_prec (exact_value, wp + 50);
	    mpfr_set_ui (exact_value, 3, GMP_RNDN);
	    mpfr_exp (exact_value, exact_value, GMP_RNDN);
	    mpfr_sub_ui (exact_value, exact_value, 1, GMP_RNDN);
	    break;

	case 2: /* f(t) = t*log(1+t), [a, b] = [0, 1] */
#if VERBOSE >= 1
	    printf ("f(x) = x*log(1+x), [a, b] = [0, 1]\n");
#endif
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 1, GMP_RNDN);
	    f = bailey1;
	    mpfr_set_ui (exact_value, 1, GMP_RNDN);
	    mpfr_mul_2si(exact_value, exact_value, -2, GMP_RNDN);
	    break;

	case 3: /* f(t) = sqrt(1-t^2), [a, b] = [0, 1] */
#if VERBOSE >= 1
	    printf ("f(t) = sqrt(1-t^2), [a, b] = [0, 1]\n");
#endif
	    f = sqrtcomp;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 1, GMP_RNDN);
	    mpfr_const_pi (exact_value, GMP_RNDN);
	    mpfr_mul_2si (exact_value, exact_value, -2, GMP_RNDN);
	    break;

	case 4: /* f(t) = sin(sin(x)), [a, b] = [0, 1] */
#if VERBOSE >= 1
	    printf ("f(t) = sin(sin(x)), [a, b] = [0, 1]\n");
#endif
	    f = sinsin;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 1, GMP_RNDN);
	    break;
	
	case 5:
#if VERBOSE >= 1
	    printf ("f(t) = sin(t) exp(cos(t)), [a, b] = [0, pi]\n");
#endif
	    f = sin_exp_cos;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_const_pi (b, GMP_RNDN);
	    mpfr_set_ui (exact_value, 1, GMP_RNDN);
	    mpfr_exp (exact_value, exact_value, GMP_RNDN);
	    mpfr_init2 (tmp, wp);
	    mpfr_ui_div (tmp, 1, exact_value, GMP_RNDN);
	    mpfr_sub (exact_value, exact_value, tmp, GMP_RNDN);
	    mpfr_clear (tmp);
	    break;

	case 6:
#if VERBOSE >= 1
	    printf ("f(t) = 1/(1+t^2), [a, b] = [0, 1]\n");
#endif
	    f = darctan;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 1, GMP_RNDN);
	    mpfr_const_pi (exact_value, GMP_RNDN);
	    mpfr_div_2exp (exact_value, exact_value, 2, GMP_RNDN);
	    break;

	case 7:
#if VERBOSE >= 1
	    printf ("f(t) = 1/(1 + 10^10 * t^2), [a, b] = [0, 1]\n");
#endif
	    f = scaled_darctan;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 1, GMP_RNDN);
	    mpfr_set_ui (exact_value, 100000, GMP_RNDN);
	    mpfr_atan (exact_value, exact_value, GMP_RNDN);
	    mpfr_div_ui (exact_value, exact_value, 100000, GMP_RNDN);
	    break;

	case 8:
#if VERBOSE >= 1
	    printf ("f(t) = exp(-x^100), [a, b] = [0, 1.1]\n");
#endif
	    f = expx100;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 11, GMP_RNDN);
	    mpfr_div_ui (b, b, 10, GMP_RNDN);
	    break;

	case 9:
#if VERBOSE >= 1
	    printf ("f(t) = x^2 sin(x^3), [a, b] = [0, 10]\n");
#endif
	    f = x2sx3;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 10, GMP_RNDN);
	    mpfr_set_ui (exact_value, 1000, GMP_RNDN);
	    mpfr_cos (exact_value, exact_value, GMP_RNDN);
	    mpfr_ui_sub (exact_value, 1, exact_value, GMP_RNDN);
	    mpfr_div_ui (exact_value, exact_value, 3, GMP_RNDN);
	    break;

	case 10:
#if VERBOSE >= 1
	    printf ("f(t) = sqrt(x), [a, b] = [0, 1]\n");
#endif
	    f = sqrt_wrapper;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 1, GMP_RNDN);
	    mpfr_set_ui (exact_value, 2, GMP_RNDN);
	    mpfr_div_ui (exact_value, exact_value, 3, GMP_RNDN);
	    break;

	case 11:
#if VERBOSE >= 1
	    printf ("f(t) = max(sin(x), cos(x)), [a, b] = [0, 1]\n");
#endif
	    f = max_sin_cos;
	    mpfr_set_ui (a, 0, GMP_RNDN);
	    mpfr_set_ui (b, 1, GMP_RNDN);
	    break;
	
	case 12:
	    printf ("f(t) = exp(-t^2) * ln(t), [a, b] = [17, 42]\n");
	    f = expln;
	    mpfr_set_ui (a, 17, GMP_RNDN);
	    mpfr_set_ui (b, 42, GMP_RNDN);
	    break;
	
	default: /* only compute the coefficients of GL */
#ifdef WANT_TIMINGS
	    total_time = cputime ();
#endif
	    coeffsgl = malloc (((n + 1) / 2) * sizeof(*coeffsgl));
	    roots = malloc ((n / 2) * sizeof(*coeffsgl));
	    for (i = 0; i < n / 2; i++)
		mpfr_init2 (roots[i], wp);
	    for (i = 0; i < (n+1) / 2; i++)
		mpfr_init2 (coeffsgl[i], wp);
	    compute_glc (coeffsgl, roots, n, wp, options, &stats);
	    for (i = 0; i < n / 2; i++) {
		mpfr_out_str (stdout, DBASE, 0, roots[i], GMP_RNDN);
		printf (" [x%d]\n", i);
	    }
#ifdef WANT_TIMINGS
	    total_time = cputime () - total_time;
	    printf ("%d [Total CPU time]\n", total_time);
#if VERBOSE >= 1
	    printf ("%d [Uspensky time]\n", stats.usp_time);
	    printf ("%d [Root refinement time]\n", stats.refine_time);
	    printf ("%d [Coefficient computation time]\n", stats.coeff_time);
	    printf ("%d [Useful newton time]\n", stats.newton_good);
	    printf ("%d [Useless newton time]\n", stats.newton_bad);
	    printf ("%d [dicho preparation time]\n", stats.dich_prep_time);
	    printf ("%d [dicho time]\n", stats.dicho_time);
	    printf ("%d [refinement failures]\n", stats.refinement_failures);
#endif
#endif
	    return 0; /* XXX clean variables */
    }
    mpfr_init (roundoff_error);
    mpfr_init (total_error);
    /* Common code: compute the integral */
    if (method == GL) {
	if (adaptive == 1) {
	    adaptive_quad (res, a, b, f, wp, gl_wrapper, &gparams);
	} else {
	    compose_gl_stats (res, roundoff_error, a, b, m, n, f, diffbound, wp, options, &stats);
	    gl_error (math_error, a, b, m, n, ndiffbound);
	}
    } else if (method == NC) {
	if (adaptive == 1) {
	    adaptive_quad (res, a, b, f, wp, nc_wrapper, &ncparams);
	} else {
	    compose_nc (res, roundoff_error, a, b, m, n, f, diffbound, wp);
	    nc_closed_error (math_error, a, b, m, n, ndiffbound);
	}
    }
    mpfr_out_str (NULL, DBASE, 0, res, GMP_RNDN); printf (" [val]\n");
    mpfr_sub (exact_value, exact_value, res, GMP_RNDN);
    mpfr_abs (exact_value, exact_value, GMP_RNDN);
    mpfr_out_str (NULL, DBASE, 0, exact_value, GMP_RNDN);
    printf (" [measured error]\n");
    mpfr_out_str (NULL, DBASE, 0, math_error, GMP_RNDN);
    printf (" [method error]\n");
    mpfr_out_str (NULL, DBASE, 0, roundoff_error, GMP_RNDN);
    printf (" [roundoff error]\n");
    mpfr_add (total_error, roundoff_error, math_error, GMP_RNDU);
    mpfr_out_str (NULL, DBASE, 0, total_error, GMP_RNDN);
    printf (" [total error]\n");
    if (mpfr_cmp (math_error, roundoff_error) > 0)
	printf ("Insufficient order\n");
    else
	printf ("Insufficient prec\n");
#ifdef WANT_TIMINGS
    total_time = cputime () - total_time;
#endif

#ifdef WANT_TIMINGS
#if VERBOSE >= 1
    printf ("%d [Total CPU time]\n", total_time);
#endif
#endif
#if VERBOSE >= 1
    printf ("%d [Uspensky time]\n", stats.usp_time);
    printf ("%d [Root refinement time]\n", stats.refine_time);
    printf ("%d [Coefficient computation time]\n", stats.coeff_time);
    printf ("%d [Useful newton time]\n", stats.newton_good);
    printf ("%d [Useless newton time]\n", stats.newton_bad);
    printf ("%d [dicho preparation time]\n", stats.dich_prep_time);
    printf ("%d [dicho time]\n", stats.dicho_time);
    printf ("%d [refinement failures]\n", stats.refinement_failures);
#endif

    mpfr_clear (a);
    mpfr_clear (b);
    mpfr_clear (exact_value);
    mpfr_clear (res);
    return 0;
}

/* gl_gen.c -- generic gauss-legendre quadrature functions.

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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gmp.h>
#include <mpfr.h>
#include "usp.h"
#include "utils.h"
#include "crq.h"
#include "crq-impl.h"

#define TEST2

void affiche_interval (interval *I)
{
    mpz_out_str (stdout, DBASE, I->c);
    printf (" [I->c]\n%ld [I->k]\n", I->k);
}

void affiche_poly (mpz_t *P, unsigned long d)
{
    unsigned long j;
    printf ("gnuplot: p(x) = ");
    for (j = 0; j <= d; j++) {
	if (mpz_sgn (P[j]) > 0)
	    printf ("+");
	mpz_out_str (stdout, DBASE, P[j]);
	printf (" * x**%ld ", j);
    }
    printf ("\n");
    printf ("maple : ");
    for (j = 0; j <= d; j++) {
	if (mpz_sgn (P[j]) > 0)
	    printf ("+");
	mpz_out_str (stdout, DBASE, P[j]);
	printf (" * X^%ld ", j);
    }
    printf ("\n");
}

void LegendrePolynomial (mpz_t *p_odd, mpz_t *p_even, unsigned long *d,
			 unsigned long n)
{
    unsigned long alpha;
    unsigned long i, j;
    unsigned int parity;
    unsigned long val, minval;
    unsigned long u, v;
    mpz_t foo, bar;
    mpz_t *p[2];

    mpz_init (foo);
    mpz_init (bar);

    p[0] = p_even;
    p[1] = p_odd;
    
    mpz_set_ui (p_even[0], 1);
    d[0] = 0; /* p0 = 1 */
    mpz_set_ui (p_odd[0], 1);
    d[1] = 0; /* p1 = x */

    parity = 0;
    for (i = 2; i <= n; i++) { /* We already have p0 and p1 */
	alpha = d[1 - parity] - d[parity];
	/* Special case for p[parity][0] */
	if (parity == 0) {
	    mpz_mul_ui (foo, p[parity][0], i - 1);
	    mpz_mul_2exp (foo, foo, alpha);
	    mpz_neg (p[parity][0], foo);
	}
	else {
	    mpz_mul_ui (foo, p[parity][0], i - 1);
	    mpz_mul_2exp (foo, foo, alpha);
	    mpz_mul_ui (bar, p[1 - parity][0], 2 * i - 1);
	    mpz_sub (p[parity][0], bar, foo);
	}

	minval = mpz_scan1 (p[parity][0], 0);

	for (j = 1; j <= i/2; j++) {
	    mpz_mul_ui (foo, p[parity][j], i - 1);
	    mpz_mul_2exp (foo, foo, alpha);

	    mpz_mul_ui (bar, p[1 - parity][j + parity - 1], 2 * i - 1);
	    mpz_sub (p[parity][j], bar, foo);
	    val = mpz_scan1 (p[parity][j], 0);
	    if (val < minval)
		minval = val;
	}
	/* We need val_2(i) */
	v = i;
	u = 0;
	while ((v & 1) == 0) {
	    v = v >> 1;
	    u++;
	}
	d[parity] = d[1 - parity] + u - minval;
	for (j = 0; j <= i/2; j++) {
	    mpz_divexact_ui (p[parity][j], p[parity][j], v);
	    mpz_tdiv_q_2exp (p[parity][j], p[parity][j], minval);
	}
	parity = 1 - parity;
    }
    mpz_clear (foo);
    mpz_clear (bar);
}

void refine_root (mpfr_t res, interval *I, mpz_t *p, mpz_t *pp,
		  mpfr_t *pl, mpfr_t *pr, unsigned long n, mp_prec_t wp,
		  struct crq_gl_opts options, struct crq_stats *timings)
{
    if (timings != NULL)
	timings->refine_time -= cputime ();
#if VERBOSE >= 3
    printf ("REFINE ROOT\n");
#endif

    switch (options.ref) {
	case DICHO:
#if VERBOSE >= 3
	    printf ("CAS DICHO\n");
#endif
	    if (options.dicho == FULL) {
		/* We use *only* naive dichotomy for root refinement */
		root_dicho (res, I, p, n, wp, timings);
	    } else
		root_dicho_trunc (res, I, p, pl, pr, n, wp);
	    break;
	
	case NEWTON_INT:
#if VERBOSE >= 3
	    printf ("CAS NEWTONINT\n");
#endif
	case NEWTON:
#if VERBOSE >= 3
	    printf ("CAS NEWTON\n");
#endif
	    root_newton (res, I, p, pp, n, wp, options, timings);
	    break;
	
	case SECANTE:
#if VERBOSE >= 3
	    printf ("CAS SECANTE\n");
#endif
	    root_secante (res, I, p, n, wp, options, timings);
	    break;
    }
    if (timings != NULL)
	timings->refine_time += cputime ();
}

/* Computes coefficients and roots for n-points Gauss-Legendre integration
 * method. Assumes b[0]...b[n/2] roots[0]...roots[n/2] are initialized.
 * roots[k] is the (k+1)th positive root of Ln(x), and b[k] is the coefficient
 * associated with the (k+1)th positive root of Ln(x). If n is odd, b[n/2+1] is
 * the coefficient associated with the root 0.
 * options is a bit mask that selects the isolation and refinement method
 * and options:
 * 
 */

void compute_glc (mpfr_t *b, mpfr_t *roots, unsigned long n, mp_prec_t wp,
		  struct crq_gl_opts options, struct crq_stats *timings)
{
    /* First compute the Legendre polynomials. Notice that for n odd,
     * only the odd coefficients of p_n are non zero. Same for n even.
     * We have P_odd(X) = 2^(-d_odd) (p_odd[0].X + p_odd[1].X^3 + ...)
     * and     P_even(X) = 2^(-d_even) (p_even[0] + p_even[1].X^2 + ...)
     * P_n has floor(n/2) + 1 coefficients
     */
    mpz_t *p_even, *p_odd, *p_copy, *pp_copy;
    mpfr_t *pl, *pr;
    unsigned long d[2];
    mpz_t *p[2];
    unsigned long i;
    unsigned int parity;
    unsigned int nbroot;
    interval *rootsI;
    mpz_t foo, bar;
    mpfr_t baz, qux, gee, gorp;
    mpfr_t err1, err2, err3, err4, dbaz;
    mp_prec_t extra_prec = 4;
    int loaded = 0;

    mpfr_init (err1);
    mpfr_init (err2);
    mpfr_init (err3);
    mpfr_init (err4);
    mpfr_init (dbaz);
    mpfr_init (gorp);
    
    p_even = malloc ((n / 2 + 1) * sizeof(*p_even));
    p_odd = malloc ((n / 2 + 1) * sizeof(*p_even));
    p_copy = malloc ((n / 2 + 1) * sizeof(*p_even));
    pp_copy = malloc ((n / 2) * sizeof(*p_even));

    for (i = 0; i <= n/2; i++) {
	mpz_init (p_even[i]);
	mpz_init (p_odd[i]);
	mpz_init (p_copy[i]);
	mpz_set_ui (p_even[i], 0);
	mpz_set_ui (p_odd[i], 0);
    }
    for (i = 0; i < n/2; i++)
	mpz_init (pp_copy[i]);

    mpz_init (foo);
    mpz_init (bar);

    p[0] = p_even;
    p[1] = p_odd;

    /* Init p_even and p_odd */

#if VERBOSE >= 1
    fprintf (stderr, "Computing polynomials...\n");
#endif
    parity = (n+1) % 2;
    if ((options.iso == FROMFILE) && (options.file != NULL)) {
	fscanf (options.file, "%lu\n", d + (1 - parity));
	load_poly (options.file, p[1 - parity], n / 2);
	loaded = 1;
	rootsI = malloc ((n/2) * sizeof (*rootsI));
	for (i = 0; i < n/2; i++) {
	    mpz_init (rootsI[i].c);
	    load_interval (options.file, rootsI + i);
	}
    } else  { /* Compute LegendrePolynomial */
	if (timings != NULL)
	    timings->pol_time -= cputime ();
	LegendrePolynomial (p_odd, p_even, d, n);
	if (options.file != NULL) {
	    fprintf (options.file, "%lu\n", d[1 - parity]);
	    dump_poly (options.file, p[1 - parity], n / 2);
	}
	if (timings != NULL)
	    timings->pol_time += cputime ();
    }

#if VERBOSE >= 2
    affiche_poly (p[1 - parity], n/2);
#endif

    /* We need the first derivative of p_n in any case. */
    deriv (p[parity], p[1 - parity], n/2);


    if (options.iso != FROMFILE) {
	/* Uspensky modifies polynomial argument */
	for (i = 0; i <= n/2; i++)
	    mpz_set (p_copy[i], p[1 - parity][i]);
	for (i = 0; i < n/2; i++)
	    mpz_set (pp_copy[i], p[parity][i]);
	if (timings != NULL)
	    timings->usp_time = cputime ();
	if (options.iso == DOUBLE)
	    rootsI = Uspensky_couple (p_copy, pp_copy, n/2, n/2 - 1, &nbroot);
	else
	    rootsI = Uspensky (p_copy, n/2, &nbroot);
	if (timings != NULL)
	    timings->usp_time = cputime () - timings->usp_time;

	/* If we have only one root in [0,1], which can happen with
	 * degree 3 or 2, Uspensky can return k < 0, which makes dicho fail.
	 */
	if (nbroot == 1) {
	    mpz_set_ui (rootsI[0].c, 0);
	    rootsI[0].k = 0;
	}
    }

    mpfr_init2 (baz, wp);
    mpfr_init2 (qux, wp);
    mpfr_init2 (gee, wp);

    if (options.dicho == TRUNC) {
	/* TODO : use mpfi for intervals */
	pl = malloc ((n/2 + 1) * sizeof (*pl));
	pr = malloc ((n/2 + 1) * sizeof (*pr));
	for (i = 0; i <= n/2; i++) {
	    mpfr_init2 (pl[i], wp + 42);
	    mpfr_init2 (pr[i], wp + 42);
	    mpfr_set_z (pl[i], p[1 - parity][i], GMP_RNDD);
	    mpfr_set_z (pr[i], p[1 - parity][i], GMP_RNDU);
	}
    }

#if VERBOSE >= 1
    fprintf (stderr, "Refining roots and computing coefficients...\n");
#endif
#if VERBOSE >= 3
    for (i = 0; i < n/2; i++) {
	double foo;
	if (options.iso == DOUBLE)
	    foo = mpz_get_ui (rootsI[2 * i].c) / ((float) (1 << rootsI[2 * i].k));
	else
	    foo = mpz_get_ui (rootsI[i].c) / ((float) (1 << rootsI[i].k));
	printf ("gnuplot: set arrow %ld from %lf, graph 0 to %lf, graph 1 nohead lt -1\n",
		i + 1, foo, foo);
    }
#endif

    for (i = 0; i < n/2; i++) {
#if VERBOSE >= 1
	fprintf (stderr, "%ld of %ld\n", i+1, n/2);
#endif
	mpfr_set_nan (roots[i]);
	/* Refine root */
	do {
	    if (options.iso == DOUBLE)
		refine_root (roots[i], rootsI + 2*i, p[1 - parity], p[parity], pl, pr,
		     n/2, wp + extra_prec, options, timings);
	    else
		refine_root (roots[i], rootsI + i, p[1 - parity], p[parity], pl, pr,
		     n/2, wp + extra_prec, options, timings);
#if VERBOSE >= 3
	mpfr_out_str (stdout, DBASE, 0, roots[i], GMP_RNDN);
	printf (" [root2 %ld]\n", i);
#endif
	/* We have the square of the root x_i. Compute the associated
	 * coefficient */
	/* Error on x_i^2 is less than one ulp, so error on \sqrt(x_i) is
	 * less than 7/2 ulp(\sqrt(xi)) */
	if (timings != NULL)
	    timings->coeff_time -= cputime ();
	mpfr_sub_ui (baz, roots[i], 1, GMP_RNDN);
	mpfr_neg (baz, baz, GMP_RNDN); /* 1 - xi^2 */
	mpfr_setulp (dbaz, baz);
	mpfr_setulp (gorp, roots[i]);
	mpfr_mul_2si (gorp, gorp, 1, GMP_RNDN);
	mpfr_add (dbaz, dbaz, gorp, GMP_RNDU);
	mpfr_mul_2si (dbaz, dbaz, 1 - mpfr_get_exp (baz), GMP_RNDU);
	/* dbaz such that |\widehat{baz} - baz| <= 1/2 dbaz \ulp(\widehat{baz})
	 * or |\widehat{baz} - baz | <= dbaz (baz)
	 */
#if 0
	debug_mpfr_value (dbaz, "dbaz");
#endif
	horner_eval (qux, roots[i], p[parity], n/2 - 1, dbaz, err1);
#if 0
	debug_mpfr_value (err1, "erreur sur Q'(x^2)");
	debug_mpfr_value (qux, "Q'(x^2)");
#endif
	/* qux = Q'n(x^2), with absolute error at most err1 */
	mpfr_mul_2si (dbaz, dbaz, -1, GMP_RNDU);  /* back to absolute error */

	if (n & 1) {
	    mpfr_setulp (err2, roots[i]);
	    mpfr_mul_2si (err2, err2, 1 - mpfr_get_exp (roots[i]), GMP_RNDU);
	    horner_eval (gee, roots[i], p[1 - parity], n/2, err2, err3);
	    /* gee = Q(x^2) with error at most err3 */
	    mpfr_setulp (err2, roots[i]);
	    mpfr_mul (err4, err2, err3, GMP_RNDU);
	    mpfr_mul (err2, err2, gee, GMP_RNDN); /* no rounding */
	    mpfr_abs (err2, err2, GMP_RNDN);
	    mpfr_add (err4, err4, err2, GMP_RNDU);
	    mpfr_mul (err2, roots[i], err3, GMP_RNDU);
	    mpfr_add (err4, err4, err2, GMP_RNDU);
	    mpfr_mul (qux, qux, roots[i], GMP_RNDN);
	    mpfr_mul_2si (qux, qux, 1, GMP_RNDN);
	    mpfr_mul_2si (err4, err4, 1, GMP_RNDN);
	    /* qux = 2x^2 Q'(x^2) with error at most err4 */
	    mpfr_add (qux, qux, gee, GMP_RNDN);
	    mpfr_add (err4, err3, err4, GMP_RNDU);
	    mpfr_setulp (err2, qux);
	    mpfr_add (err4, err4, err2, GMP_RNDU);
	    /* qux = P'(x^2) with error at most err4 */
	    mpfr_sqrt (roots[i], roots[i], GMP_RNDN);
	}
	else {
	    mpfr_sqrt (roots[i], roots[i], GMP_RNDN);
	    mpfr_setulp (err2, roots[i]);
	    mpfr_mul_ui (err2, err2, 4, GMP_RNDU);
	    /* roots[i] is x_i with error at most err2 */
	    mpfr_mul (err4, err1, err2, GMP_RNDU);
	    mpfr_mul (err3, err1, roots[i], GMP_RNDU);
	    mpfr_add (err4, err4, err3, GMP_RNDU);
	    mpfr_mul (err3, err2, qux, GMP_RNDU);
	    mpfr_add (err4, err4, err3, GMP_RNDU);
	    mpfr_mul (qux, qux, roots[i], GMP_RNDN);
	    mpfr_mul_2si (qux, qux, 1, GMP_RNDN);
	    mpfr_mul_2si (err4, err4, 1, GMP_RNDN);
	    /* qux = P'(x^2) with error at most err4 */
	} /* pn'(xi) */
	mpfr_sqr (err1, err4, GMP_RNDU);
	mpfr_mul_2si (err2, err4, 1, GMP_RNDU);
	mpfr_mul (err2, err2, qux, GMP_RNDU);
	mpfr_abs (err2, err2, GMP_RNDU);
	mpfr_add (err4, err1, err2, GMP_RNDU);
	mpfr_sqr (qux, qux, GMP_RNDN);
	mpfr_setulp (err1, qux);
	mpfr_mul_2si (err1, err1, -1, GMP_RNDN);
	mpfr_add (err4, err4, err1, GMP_RNDU);
	/* qux is P'^2(x_i) with error at most err4 */
	mpfr_mul (err2, err4, baz, GMP_RNDU);
	mpfr_mul (err3, dbaz, qux, GMP_RNDU);
	mpfr_add (err1, err2, err3, GMP_RNDU);
	mpfr_mul (err2, dbaz, err4, GMP_RNDU);
	mpfr_add (err1, err1, err2, GMP_RNDU);
	mpfr_mul (baz, baz, qux, GMP_RNDN);
	mpfr_setulp (err2, baz);
	mpfr_mul_2si (err2, err2, -1, GMP_RNDN);
	mpfr_add (err1, err1, err2, GMP_RNDN);
	/* baz is P'^2(x) (1-x^2) with error at most err1 */
	/* translate err1 as ulp(baz) : */
	mpfr_mul_2si (err1, err1, mpfr_get_prec (baz) - mpfr_get_exp (baz),
		      GMP_RNDU);
	mpfr_set_ui (qux, 2, GMP_RNDN);
	mpfr_div (b[i], qux, baz, GMP_RNDN);
	mpfr_mul_2si (b[i], b[i], 2 * d[1 - parity], GMP_RNDN);
	mpfr_mul_2si (err1, err1, 3, GMP_RNDU);
	mpfr_add_ui (err1, err1, 1, GMP_RNDU);
	mpfr_mul_2si (err1, err1, -1, GMP_RNDU);
	/* Hopefully, b[i] is now w_i with error at most err1 ulp(b[i]).
	 * See if we can round */
	mpfr_setulp (err2, roots[i]);
	mpfr_mul_ui (err2, err2, 4, GMP_RNDU);
	mpfr_add_ui (roots[i], roots[i], 1, GMP_RNDN);
	mpfr_setulp (err3, roots[i]);
	mpfr_add (err2, err2, err3, GMP_RNDU);
	mpfr_mul_2si (roots[i], roots[i], -1, GMP_RNDN);
	mpfr_mul_2si (err2, err2, -1 + mpfr_get_prec (roots[i]) - mpfr_get_exp (roots[i]), GMP_RNDN);
	/* roots[i] is (1+x_i) / 2  with error at most err2 ulp(roots[i]) */

	if (timings != NULL)
	    timings->coeff_time += cputime ();
	if ((! mpfr_can_round (b[i], mpfr_get_prec (b[i])
				  - mpfr_get_exp (err1) - 1, GMP_RNDN,
		 			     GMP_RNDN, wp)) ||
	    (! mpfr_can_round (roots[i], mpfr_get_prec (roots[i])
					  - mpfr_get_exp (err2) - 1, GMP_RNDN,
					  GMP_RNDN, wp))) {
	    extra_prec *= 2;
	    mpfr_set_prec (qux, wp + extra_prec);
	    mpfr_set_prec (baz, wp + extra_prec);
	    mpfr_set_prec (gee, wp + extra_prec);
	    mpfr_set_prec (b[i], wp + extra_prec);
	    mpfr_prec_round (roots[i], wp + extra_prec, GMP_RNDN);
	    continue;
	} else
	    break;
	} while (1);

	mpfr_prec_round (b[i], wp, GMP_RNDN);
	mpfr_prec_round (roots[i], wp, GMP_RNDN);
#if 0
	printf ("%ld [extra_prec]\n", extra_prec);
#endif
	abort_if_nan (b[i], __FILE__, __LINE__);
	
#if VERBOSE >= 3
	mpfr_out_str (stdout, DBASE, 0, roots[i], GMP_RNDN);
	printf (" [root %ld]\n", i);
	mpfr_out_str (stdout, DBASE, 0, b[i], GMP_RNDN);
	printf (" [b %ld]\n", i);
#endif
	extra_prec /= 2;
    }

    /* Compute coefficient for zero if needed */
    if (n & 1) {
	if (timings != NULL)
	    timings->coeff_time -= cputime ();
	mpfr_set_z (baz, p[1 - parity][0], GMP_RNDN);
	mpfr_sqr (baz, baz, GMP_RNDN);
	mpfr_set_ui (qux, 2, GMP_RNDN);
	mpfr_div (baz, qux, baz, GMP_RNDN);
	mpfr_mul_2si (b[n/2], baz, 2 * d[1 - parity], GMP_RNDN);
	if (timings != NULL)
	    timings->coeff_time += cputime ();
#if VERBOSE >= 3
	mpfr_out_str (stdout, DBASE, 0, b[n/2], GMP_RNDN);
	printf (" [b zero]\n");
#endif
	if (timings != NULL)
	    timings->coeff_time += cputime ();
    }
    /* Save rootsI if requested */

    if ((loaded == 0) && (options.file != NULL)) {
	for (i = 0; i < n/2; i++) {
	    if (options.iso == DOUBLE)
		dump_interval (options.file, rootsI[2 * i]);
	    else
		dump_interval (options.file, rootsI[i]);
	}
    }

    /* Free allocated memory */
    
    if (options.dicho == TRUNC) {
	for (i = 0; i <= n/2; i++) {
	    mpfr_clear (pl[i]);
	    mpfr_clear (pr[i]);
	}
	free (pl);
	free (pr);
    }

    for (i = 0; i <= n/2; i++) {
	mpz_clear (p_even[i]);
	mpz_clear (p_odd[i]);
	mpz_clear (p_copy[i]);
    }
    for (i = 0; i < n/2; i++)
	mpz_clear (pp_copy[i]);
    if (options.iso == DOUBLE) {
	for (i = 0; i < (n - 1 - (n & 1)); i++)
	    mpz_clear (rootsI[i].c);
    } else {
	for (i = 0; i < n/2; i++)
	    mpz_clear (rootsI[i].c);
    }
    free (p_even);
    free (p_odd);
    free (p_copy);
    free (pp_copy);
    free (rootsI);
    mpz_clear (foo);
    mpz_clear (bar);
    mpfr_clear (baz);
    mpfr_clear (qux);
    mpfr_clear (gee);

    mpfr_clear (err1);
    mpfr_clear (err2);
    mpfr_clear (err3);
    mpfr_clear (err4);
    mpfr_clear (dbaz);
    mpfr_clear (gorp);
}

/* Error bound for m-composed GL generic rule for n points on [a, b].
 * This accoutns only for the method error.
 * M is a bound of |f^(2n)| on [a, b]
 */

void gl_error (mpfr_t res, mpfr_t a, mpfr_t b, unsigned int m,
	       unsigned int n, mpfr_t M)
{
    mpfr_t foo;

    mpfr_init2 (foo, mpfr_get_prec (res));

    mpfr_sub (res, b, a, GMP_RNDU);
    mpfr_div_ui (res, res, m, GMP_RNDU);
    mpfr_pow_ui (res, res, 2 * n + 1, GMP_RNDU);
    mpfr_fac_ui (foo, n, GMP_RNDU);
    mpfr_pow_ui (foo, foo, 4, GMP_RNDU);
    mpfr_mul (res, res, foo, GMP_RNDU);
    mpfr_fac_ui (foo, 2 * n, GMP_RNDD);
    mpfr_pow_ui (foo, foo, 3, GMP_RNDD);
    mpfr_mul_ui (foo, foo, 2 * n + 1, GMP_RNDD);
    mpfr_div (res, res, foo, GMP_RNDU);
    mpfr_mul (res, res, M, GMP_RNDU);
    mpfr_mul_ui (res, res, m, GMP_RNDU);
    mpfr_clear (foo);
}

/* Compute the Horner evaluation of the polynomial p of degree n at
 * x. Assume x is given as ~x where ~x = (1+foo)x, |foo| <= inerr.
 * Return outerr such that |p(x) -~p(~x) | <= outerr.
 */

void horner_eval (mpfr_t res, mpfr_t x, mpz_t *p, unsigned long n,
		  mpfr_t inerr, mpfr_t outerr)
{
    mpfr_t foo, bar;
    mpfr_t delta;
    mpfr_set_z (res, p[n], GMP_RNDN);
    mp_prec_t prec_x, prec_res;

    mpfr_setulp (outerr, res);
    if (n == 0)
	return;

    prec_res = mpfr_get_prec (res);
    prec_x = mpfr_get_prec (x);

    mpfr_init (delta);
    mpfr_init (foo);
    mpfr_init (bar);
    mpfr_set_ui (foo, 1, GMP_RNDN);
    mpfr_mul_2si (foo, foo, - prec_res, GMP_RNDN);
    mpfr_add (delta, foo, inerr, GMP_RNDU);
    mpfr_mul (foo, foo, inerr, GMP_RNDU);
    mpfr_add (delta, delta, foo, GMP_RNDU);
    mpfr_set_ui (foo, 1, GMP_RNDN);
    mpfr_mul_2si (foo, foo, - prec_res, GMP_RNDN);
    if (mpfr_cmp (delta, foo) < 0)
	mpfr_set (delta, foo, GMP_RNDU);
    mpfr_set_ui (outerr, 0, GMP_RNDN);

    mpfr_prec_round (res, prec_res, GMP_RNDN);
    do {
	n--;
#if 0
	/* Prec is enough, rounding mode doesn't matter */
	mpfr_prec_round (res, prec_res + prec_x, GMP_RNDD);
#endif
	mpfr_mul (res, res, x, GMP_RNDN);
	mpfr_abs (foo, res, GMP_RNDU);
	mpfr_add_z (res, res, p[n], GMP_RNDN);
	mpfr_abs (bar, res, GMP_RNDU);
	mpfr_add (foo, foo, bar, GMP_RNDU);
	mpfr_mul (foo, foo, delta, GMP_RNDU);
	mpfr_abs (bar, x, GMP_RNDU);
	mpfr_mul (outerr, outerr, bar, GMP_RNDU);
	mpfr_add (outerr, outerr, foo, GMP_RNDU);
#if 0
	mpfr_prec_round (res, prec_res, GMP_RNDN);
#endif
    } while (n > 0);
    mpfr_clear (foo);
    mpfr_clear (bar);
    mpfr_clear (delta);
}

/* Compute p''(x) */
void horner_evalp2 (mpfr_t res, mpfr_t x, mpz_t *p, unsigned long n)
	
{
    mpfr_t foo;
    mpfr_set_z (res, p[n], GMP_RNDN);

    if (n < 2) {
	mpfr_set_ui (res, 0, GMP_RNDN);
	return;
    }
    mpfr_set_z (res, p[n], GMP_RNDN);
    mpfr_mul_ui (res, res, n, GMP_RNDN);
    mpfr_mul_ui (res, res, n - 1, GMP_RNDN);
    mpfr_init2 (foo, mpfr_get_prec (res));
    printf ("Debut debuggage de mongolien\n");
    mpfr_out_str (stdout, 10, 0, res, GMP_RNDN);
    printf ("\n");
    
    n--;
    while (n > 1) {
	mpfr_mul (res, res, x, GMP_RNDN);
	mpfr_set_z (foo, p[n], GMP_RNDN);
	mpfr_mul_ui (foo, foo, n, GMP_RNDN);
	mpfr_mul_ui (foo, foo, n - 1, GMP_RNDN);
	mpfr_add (res, res, foo, GMP_RNDN);
	n--;
	mpfr_out_str (stdout, 10, 0, res, GMP_RNDN);
	printf ("\n");
    }
    mpfr_clear (foo);
    printf ("Fin debuggage de mongolien\n");
}

void horner_evalb (mpfr_t res, mpfr_t x, mpz_t *p, unsigned long n)
{
    mpfr_set_z (res, p[n], GMP_RNDN);

    if (n == 0)
	return;

    do {
	n--;
	mpfr_mul (res, res, x, GMP_RNDN);
	mpfr_add_z (res, res, p[n], GMP_RNDN);
    } while (n > 0);
}

void deriv (mpz_t *ret, mpz_t *p, unsigned long n) {
    unsigned long i;

    if (n == 0) {
	mpz_set_ui (ret[0], 0);
	return;
    }

    for (i = 0; i < n; i++)
	mpz_mul_ui (ret[i], p[i + 1], i + 1);
}

void error_on_yi (mpfr_t res, mpfr_t yi, mpfr_t xi, mpfr_t wi, mpfr_t m,
		  unsigned long wp)
{
    mpfr_t foo;
    mpfr_t bar;

    mpfr_init (foo);
    mpfr_init (bar);

#if 0
    printf ("m = ");
    mpfr_out_str (stdout, DBASE, 0, m, GMP_RNDN); printf ("\n");
    printf ("xi = ");
    mpfr_out_str (stdout, DBASE, 0, xi, GMP_RNDN); printf ("\n");
    printf ("wi = ");
    mpfr_out_str (stdout, DBASE, 0, wi, GMP_RNDN); printf ("\n");
    printf ("yi = ");
    mpfr_out_str (stdout, DBASE, 0, yi, GMP_RNDN); printf ("\n");
#endif

    mpfr_set_ui (foo, 5, GMP_RNDU);
    mpfr_set_ui (bar, 1, GMP_RNDU);
    mpfr_mul_2si (foo, foo, -1, GMP_RNDU);
    mpfr_mul_2si (bar, bar, -wp, GMP_RNDU);
    mpfr_add (foo, foo, bar, GMP_RNDU);
    mpfr_setulp (bar, yi);
    mpfr_mul (res, foo, bar, GMP_RNDU);

    mpfr_set_ui (foo, 17, GMP_RNDU);
    mpfr_mul_2si (foo, foo, -2, GMP_RNDU);
    mpfr_set_ui (bar, 17, GMP_RNDU);
    mpfr_mul_2si (bar, bar, -wp -2, GMP_RNDU);
    mpfr_add (foo, foo, bar, GMP_RNDU);
    mpfr_mul (foo, foo, m, GMP_RNDU);
    mpfr_mul (foo, foo, wi, GMP_RNDU);
    mpfr_setulp (bar, xi);
    mpfr_mul (foo, foo, bar, GMP_RNDU);
    mpfr_add (res, res, foo, GMP_RNDU);

#if 0
    printf ("delta_yi = ");
    mpfr_out_str (stdout, DBASE, 0, res, GMP_RNDN); printf ("\n");
#endif

    mpfr_clear (foo);
    mpfr_clear (bar);
}
    

/* Closed GL generic rule for n points on [a, b] and integrand f.
 * coeffs and roots are the associated parameters computed with
 * compute_glc, m is a bound of |f'| on [a, b].
 * fxs is an auxiliary (n+1)/2 table of mpfr_t.
 * wp is the working precision.
 */

void gen_gl_closed (mpfr_t res, mpfr_t errb, mpfr_t a, mpfr_t b, unsigned int n,
		    mpfr_t *coeffs, mpfr_t *roots, aqfunc_t f, mpfr_t m,
		    mp_prec_t wp, mpfr_t *fxs[])
{
    unsigned int i;
    mpfr_t acc;
    mpfr_t apb, bma, foo, bar;
    mpfr_t maxdy, currdy;

    assert (n >= 2);

    mpfr_init2 (acc, wp);
    mpfr_init (maxdy);
    mpfr_init (currdy);
    mpfr_set_ui (acc, 0, GMP_RNDN);

    mpfr_init2 (apb, wp);
    mpfr_init2 (bma, wp);
    mpfr_init2 (foo, wp);
    mpfr_init2 (bar, wp);

    mpfr_add (apb, a, b, GMP_RNDN);
    mpfr_sub (bma, b, a, GMP_RNDN);

    for (i = 0; i < n/2; i++) {
#if 0
	debug_mpfr_value (roots[i], "V_i");
#endif
	mpfr_mul (foo, roots[i], bma, GMP_RNDN);
	mpfr_add (foo, foo, a, GMP_RNDN);
#if 0
	debug_mpfr_value (foo, "X_i");
#endif
	/*mpfr_mul_2si (foo, foo, -1, GMP_RNDN);*/
	f (*fxs[i], foo);
	mpfr_mul (*fxs[i], *fxs[i], coeffs[i], GMP_RNDN);
	
	error_on_yi (currdy, *fxs[i], foo, coeffs[i], m, wp);
	if (i == 0)
	    mpfr_set (maxdy, currdy, GMP_RNDU);
	else {
	    if (mpfr_cmp (currdy, maxdy) > 0)
		mpfr_set (maxdy, currdy, GMP_RNDU);
	}
	mpfr_ui_sub (foo, 1, roots[i], GMP_RNDN);
	mpfr_mul (foo, foo, bma, GMP_RNDN);
	mpfr_add (foo, a, foo, GMP_RNDN);
	/*mpfr_mul_2si (foo, foo, -1, GMP_RNDN);*/
	f (bar, foo);
	mpfr_mul (bar, bar, coeffs[i], GMP_RNDN);
	error_on_yi (currdy, bar, foo, coeffs[i], m, wp);
	if (mpfr_cmp (currdy, maxdy) > 0)
	    mpfr_set (maxdy, currdy, GMP_RNDU);
	mpfr_add (*fxs[i], *fxs[i], bar, GMP_RNDN);
#if VERBOSE >= 3
	mpfr_out_str (stdout, DBASE, 0, *fxs[i], GMP_RNDN); printf ("\n");
#endif
    }
    if (n & 1) {
	mpfr_mul_2si (foo, apb, -1, GMP_RNDN);
	f (*fxs[n/2], foo);
	mpfr_mul (*fxs[n/2], *fxs[n/2], coeffs[n/2], GMP_RNDN);
    }
    /* TODO : directly calling mpfr_sum is not suitable here */
#if 0
    mpfr_set (max_yi, *fxs[0], GMP_RNDU);
    for (i = 1; i < n/2; i++) {
	if (mpfr_cmp (max_yi, *fxs[i]) < 0)
	    mpfr_set (max_yi, *fxs[i], GMP_RNDU);
    }
    if (n & 1) {
	if (mpfr_cmp (max_yi, *fxs[n/2]) < 0)
	    mpfr_set (max_yi, *fxs[n/2], GMP_RNDU);
    }
#endif
    mpfr_sum (res, (mpfr_ptr *) fxs, (n+1) / 2, GMP_RNDN);
    mpfr_sub (foo, b, a, GMP_RNDN);
    mpfr_mul_2si (res, res, -1, GMP_RNDN);
    mpfr_mul (res, res, foo, GMP_RNDN);

    mpfr_sub (bma, b, a, GMP_RNDU);
    mpfr_mul (errb, maxdy, bma, GMP_RNDU);
    mpfr_mul_ui (errb, errb, n, GMP_RNDU);
    mpfr_set_ui (bar, 1, GMP_RNDU);
    mpfr_mul_2si (bar, bar, -wp, GMP_RNDU);
    mpfr_add_ui (bar, bar, 1, GMP_RNDU);
    mpfr_mul (errb, errb, bar, GMP_RNDU);

    mpfr_set_ui (bar, 1, GMP_RNDU);
    mpfr_mul_2si (bar, bar, -wp + 1, GMP_RNDU);
    mpfr_mul_ui (bar, bar, 3, GMP_RNDU);
    mpfr_add_ui (bar, bar, 9, GMP_RNDU);
    mpfr_mul_2si (bar, bar, mpfr_get_exp (res) - wp, GMP_RNDU);
    mpfr_add (errb, errb, bar, GMP_RNDU);

    mpfr_clear (acc);
    mpfr_clear (foo);
    mpfr_clear (bar);
    mpfr_clear (apb);
    mpfr_clear (bma);
}

/* Compose m times the closed n points GL method on [a, b].
 * Returns a bound on the evaluation error in errb.
 * Fills the timings struct with timings data. */
void compose_gl_stats (mpfr_t res, mpfr_t errb, mpfr_t a, mpfr_t b,
		      unsigned int m, unsigned int n, aqfunc_t f,
		      mpfr_t mp, mp_prec_t wp, struct crq_gl_opts options,
		      struct crq_stats *timings)
{
    mpfr_t *roots, *coeffs, **fxs;
    mpfr_t **gls;
    unsigned int i, step;
    mpfr_t acc, bmin, bmax, foo, bar;
    mpfr_t curr_errb;
    
    coeffs = malloc (((n + 1) / 2) * sizeof(*coeffs));
    roots = malloc ((n / 2) * sizeof(*coeffs));

    for (i = 0; i < n / 2; i++)
	mpfr_init2 (roots[i], wp);

    for (i = 0; i < (n+1) / 2; i++)
	mpfr_init2 (coeffs[i], wp);

    mpfr_init2 (acc, wp);
    mpfr_init2 (bmin, wp);
    mpfr_init2 (bmax, wp);
    mpfr_init2 (foo, wp);
    mpfr_init2 (bar, wp);
    mpfr_init (curr_errb);

    gls = malloc (m * sizeof (*gls));
    for (i = 0; i < m; i++) {
        gls[i] = malloc (sizeof (mpfr_t));
        mpfr_init2 (*gls[i], wp);
    }

    fxs = malloc (((n+1)/2) * sizeof (*fxs));
    for (i = 0; i < (n + 1) / 2; i++) {
        fxs[i] = malloc (sizeof (mpfr_t));
        mpfr_init2 (*fxs[i], wp);
    }

#if VERBOSE >= 1
    fprintf (stderr, "Computing roots and coefficients...\n");
#endif
    compute_glc (coeffs, roots, n, wp, options, timings);
#if VERBOSE >= 1
    fprintf (stderr, "Composing GL\n");
#endif

    mpfr_set (bmin, a, GMP_RNDN);
    mpfr_set_ui (acc, 0, GMP_RNDN);
    mpfr_set_ui (errb, 0, GMP_RNDN);

    for (step = 0; step < m; step++) {
	/* Compute new bmax */
	mpfr_mul_ui (foo, a, m - step - 1, GMP_RNDN);
	mpfr_mul_ui (bar, b, step + 1, GMP_RNDN);
	mpfr_add (bmax, foo, bar, GMP_RNDN);
	mpfr_div_ui (bmax, bmax, m, GMP_RNDN);
	/* Compute GL-closed n point on [bmin, bmax] */
	gen_gl_closed (*gls[step], curr_errb, bmin, bmax, n, coeffs, roots, f,
		       mp, wp, fxs);
	mpfr_add (errb, errb, curr_errb, GMP_RNDU);

	/* Set bmin for next step */
	mpfr_set (bmin, bmax, GMP_RNDN);
    }
    mpfr_sum (res, (mpfr_ptr *) gls, m, GMP_RNDN);
    /* TODO : mpfr_sum guaranties error <= 0.5 ulp but we should do it
     * by hand */
    mpfr_setulp (foo, res);
    mpfr_mul_2si (foo, foo, -1, GMP_RNDN);
    mpfr_add (errb, errb, foo, GMP_RNDN);

    for (i = 0; i < n / 2; i++)
	mpfr_clear (roots[i]);

    for (i = 0; i < (n+1) / 2; i++)
	mpfr_clear (coeffs[i]);

    for (i = 0; i < m; i++) {
        mpfr_clear (*gls[i]);
        free (gls[i]);
    }
    for (i = 0; i < (n + 1) / 2; i++) {
        mpfr_clear (*fxs[i]);
        free (fxs[i]);
    }
    free (coeffs);
    free (roots);
    free (gls);
    free (fxs);
    mpfr_clear (acc);
    mpfr_clear (bmin);
    mpfr_clear (bmax);
    mpfr_clear (foo);
    mpfr_clear (bar);
    mpfr_clear (curr_errb);
}

void compose_gl (mpfr_t res, mpfr_t errb, mpfr_t a, mpfr_t b, unsigned int m,
		 unsigned int n, aqfunc_t f, mpfr_t mp, mp_prec_t wp,
		 struct crq_gl_opts options)
{
    compose_gl_stats (res, errb, a, b, m, n, f, mp, wp, options,
			NULL);
}

#ifdef TEST1

int main (int argc, char *argv[]) {
    unsigned long i, n;
    mpfr_t *b;
    mpfr_t *roots;

    n = atoi (argv[1]);

    roots = malloc (n * sizeof(*roots));
    b = malloc (n * sizeof(*b));

    for (i = 0; i < n; i++) {
	mpfr_init (roots[i]);
	mpfr_init (b[i]);
    }

    compute_glc (b, roots, n, atoi (argv[2]));
    return 0;
}

#endif

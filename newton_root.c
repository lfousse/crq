/* newton_root.c -- Newton root refinement.

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

#include "crq-impl.h"
#include "crq.h"
#include <stdio.h>
#include <mpfr.h>
#include "usp.h"
#include "interval.h"
#include "utils.h"
#include <stdlib.h>

/* TODO : try the formula 
 *
 * x_{m,n} = (1 - 1/(8n^2) + 1/(8n^3)) cos(Pi (4m-1)/(4n+2)) + O(n^{-4})
 * x_{n,1} > x_{n,2} > ... > x_{n, n}
 */

extern void dicho (interval *I, mpz_t *P, unsigned long n, unsigned long p);

void set_middle (mpfr_t root, interval *I) {
    mp_prec_t cur_prec;

    cur_prec = mpz_sizeinbase (I->c, 2);
    mpfr_set_prec (root, cur_prec * 2 + 3);
    mpfr_set_z (root, I->c, GMP_RNDN);
    mpfr_mul_2ui (root, root, 1, GMP_RNDN);
    mpfr_add_ui (root, root, 1, GMP_RNDN);
    mpfr_div_2ui (root, root, I->k + 1, GMP_RNDN);
}

/* Aberth correction */

void aberth (mpfr_t res, mpfr_t *roots, unsigned int i, unsigned int n)
{
    unsigned int j;
    mpfr_t foo;

    mpfr_set_ui (res, 0, GMP_RNDN);
    mpfr_init2 (foo, mpfr_get_prec (res));

    for (j = 0; j < n; j++) {
	if (j == i)
	    continue;
	mpfr_sub (foo, roots[i], roots[j], GMP_RNDN);
	mpfr_ui_div (foo, 1, foo, GMP_RNDN);
	mpfr_add (res, res, foo, GMP_RNDN);
    }
    mpfr_clear (foo);
}

/* Parallel newton root finding */
void newton_para (mpfr_t *roots, mpz_t *P, mpz_t *Pp, unsigned int n,
		  mp_prec_t wp)
{
    mpfr_t acorr;
    mpfr_t px, ppx, nan;
    long good;
    unsigned int done = 0;
    unsigned int i;
    mp_prec_t cur_prec;

    mpfr_init (px);
    mpfr_init (ppx);
    mpfr_init (acorr);
    mpfr_init (nan);

    do {
	done = 0;
	for (i = 0; i < n; i++) {
	    cur_prec = mpfr_get_prec (roots[i]);
	    if (cur_prec > 2 * wp + 53)
		continue;
	    mpfr_set_prec (px, 2 * cur_prec + 2);
	    mpfr_set_prec (ppx, 2 * cur_prec + 2);
	    mpfr_set_prec (acorr, 2 * cur_prec + 2);
	    aberth (acorr, roots, i, n);
	    horner_eval (px, roots[i], P, n, nan, nan);
	    horner_eval (ppx, roots[i], Pp, n - 1, nan, nan);
	    mpfr_div (px, px, ppx, GMP_RNDN);
	    mpfr_mul (acorr, acorr, px, GMP_RNDN);
	    mpfr_ui_sub (acorr, 1, acorr, GMP_RNDN);
	    mpfr_div (px, px, acorr, GMP_RNDN);
	    if (mpfr_zero_p (px)) {
		if (mpfr_zero_p (roots[i]))
		    good = 0;
		else
		    good = mpfr_get_prec (roots[i]);
	    }
	    else {
		if (mpfr_zero_p (roots[i]))
		    good = 0;
		else
		    good = mpfr_get_exp (roots[i]) - mpfr_get_exp (px);
	    }
	    if (good < 0)
		good = 0;
	    mpfr_sub (roots[i], roots[i], px, GMP_RNDN);
	    mpfr_prec_round (roots[i], 2 * good + 53, GMP_RNDN);
	    done++;
	}
    } while (done != 0);
}

void newton_prepare (mpfr_t *roots, interval *I, unsigned int n)
{
    unsigned int i;
    mp_prec_t cur_prec;

    for (i = 0; i < n; i++) {
	cur_prec = mpz_sizeinbase (I[i].c, 2);
	mpfr_set_prec (roots[i], cur_prec * 2 + 3);
	mpfr_set_z (roots[i], I[i].c, GMP_RNDN);
	mpfr_mul_2ui (roots[i], roots[i], 1, GMP_RNDN);
	mpfr_add_ui (roots[i], roots[i], 1, GMP_RNDN);
	mpfr_div_2ui (roots[i], roots[i], I[i].k + 1, GMP_RNDN);
    }
}

#define max(a,b) ((a) < (b) ? (b) : (a))

int newton (mpfr_t root, mpz_t *p, mpz_t *pp, unsigned long d, interval *I,
	     mp_prec_t wp, int *nbiter)
{
    mpfr_t px, ppx, nan;
    mp_prec_t cur_prec;
    long good = 0;
    unsigned int i;
    mp_prec_t coeffprec;
#if VERBOSE >= 3
    mpfr_t orig_point;
    mpfr_t p2x;
#endif

    cur_prec = mpfr_get_prec (root);
    mpfr_init2 (px, cur_prec);
    mpfr_init2 (ppx, cur_prec);
    mpfr_init (nan);

    coeffprec = mpz_sizeinbase (p[0], 2);
    for (i = 1; i <= d; i++) {
	mp_prec_t tmp = mpz_sizeinbase (p[i], 2);
	coeffprec = max (coeffprec, tmp);
    }
#if VERBOSE >= 3
    printf ("Valeur initiale pour Newton: ");
    mpfr_out_str (stdout, DBASE, 0, root, GMP_RNDN);
    printf ("\n");
#endif
#if VERBOSE >= 3
    if (is_in_interval (root, I) == 0) {
	abort ();
    }
#endif

#if VERBOSE >= 3
    mpfr_init2 (orig_point, mpfr_get_prec (root));
    mpfr_set (orig_point, root, GMP_RNDN);
    mpfr_init (p2x);
#endif

    *nbiter = 0;
    while (cur_prec < 2 * wp) {
#if 0
	mpfr_set_prec (px, 2 * max(cur_prec, coeffprec));
	mpfr_set_prec (ppx, max(cur_prec, coeffprec));
#else
	mpfr_set_prec (px, 2 * cur_prec);
	mpfr_set_prec (ppx, cur_prec);
#endif
	horner_eval (px, root, p, d, nan, nan);
	horner_eval (ppx, root, pp, d - 1, nan, nan);
#if VERBOSE >= 3
	affiche_poly (pp, d - 1);
	printf ("cur_prec = %ld\n", cur_prec);
	printf ("Tangente calculee :\n");
	printf ("g(x) = (x - ");
	mpfr_out_str (stdout, 10, 0, root, GMP_RNDN);
	printf (") * ");
	mpfr_out_str (stdout, 10, 0, ppx, GMP_RNDN);
	printf (" + ");
	mpfr_out_str (stdout, 10, 0, px, GMP_RNDN);
	printf ("\n");
#endif
	mpfr_div (px, px, ppx, GMP_RNDN);
	if (mpfr_zero_p (px)) {
	    if (mpfr_zero_p (root))
		good = 0;
	    else
		good = mpfr_get_prec (root);
	}
	else {
	    if (mpfr_zero_p (root))
		good = 0;
	    else
		good = mpfr_get_exp (root) - mpfr_get_exp (px);
	}
#if VERBOSE >= 3
	printf ("Bits corrects: %ld\n", good);
#endif
	if (good < 0)
	    good = wp - 1;
	cur_prec += good;
	mpfr_prec_round (root, cur_prec, GMP_RNDN);
	mpfr_sub (root, root, px, GMP_RNDN);
	if (is_in_interval (root, I) == 0) {
#if VERBOSE >= 3
	    mpz_t mantissa;
	    unsigned int k;
	    k = - mpfr_get_z_exp (mantissa, orig_point);
	    if (sign_eval (mantissa, p, k, d) * 
		sign_eval_p2 (mantissa, p, k, d) < 0)
		printf ("P(x) et P''(x) de signe different\n");
	    else
		printf ("P(x) et P''(x) de meme signe.\n");
#endif
	    mpfr_clear (px);
	    mpfr_clear (ppx);
	    return 1;
	}
	(*nbiter)++;
    }
    /* mpfr_prec_round (root, wp, GMP_RNDN); */
    mpfr_clear (px);
    mpfr_clear (ppx);
    return 0;
}

/*#define NEWTON_BENCH*/
#ifdef NEWTON_BENCH
void root_newton (mpfr_t ret, interval *I, mpz_t *p, mpz_t *pp,
			unsigned int d, mp_prec_t wp, struct crq_gl_opts options)
{
    mp_prec_t extra_prec;
    mpz_t foo;
    mpfr_t root;
    int nret;
    int nbiter;
    interval Itmp;
    long kmin, kmax;
    int curr_failures = 0;

    mpz_init (foo);
    mpfr_init (root);

    sign_set (&(I->polysigns), SPA,
	      sign_eval (I->c, p, I->k, d));
    mpz_add_ui (foo, I->c, 1);
    sign_set (&(I->polysigns), SPB,
	      sign_eval (foo, p, I->k, d));

    extra_prec = 42;
    mpz_init (Itmp.c);
    kmin = I->k;
    Itmp.polysigns = I->polysigns;
    kmax = 100;
    /* See of kmax is sufficient */
    do {
	mpz_set (Itmp.c, I->c);
	Itmp.k = I->k;
	dicho (&Itmp, p, d, kmax);
	set_middle (root, &Itmp);
	if (options.ref == NEWTON_INT)
	    nret = newton_int (root, p, pp, d, I, wp + extra_prec, &nbiter);
	else
	    nret = newton (root, p, pp, d, I, wp + extra_prec, &nbiter);
	if (nret == 0)
	    break;
	kmax *= 2;
    } while (nret != 0);
    /* Now, dichotomy between kmin and kmax to seek the lowest
     * acceptable value for which newton converge */
    do {
	mpz_set (Itmp.c, I->c);
	Itmp.k = I->k;
	dicho (&Itmp, p, d, (kmin + kmax) / 2);
	set_middle (root, &Itmp);
	if (options.ref == NEWTON_INT)
	    nret = newton_int (root, p, pp, d, I, wp + extra_prec, &nbiter);
	else
	    nret = newton (root, p, pp, d, I, wp + extra_prec, &nbiter);
	if (nret == 0)
	    kmax = (kmax + kmin) / 2;
	else
	    kmin = (kmax + kmin) / 2;
    } while (kmax - kmin > 1);
    printf ("%ld [mink]\n", kmax);

    mpfr_set (ret, root, GMP_RNDN);
    mpfr_clear (root);
    mpz_clear (foo);
    newton_failures += curr_failures;
#if VERBOSE >= 3
    printf ("%d [curr_failures]\n", curr_failures);
#endif
}

#else
void root_newton (mpfr_t ret, interval *I, mpz_t *p, mpz_t *pp, unsigned int d,
		  mp_prec_t wp, struct crq_gl_opts options,
		  struct crq_stats *timings)
{
    mp_prec_t cur_prec, extra_prec;
    mpz_t foo;
    mpfr_t root;
    int nret;
    int nbiter;
    int curr_failures = 0;
    int time_tmp = 0;

    mpz_init (foo);
    mpfr_init (root);

    sign_set (&(I->polysigns), SPA,
	      sign_eval (I->c, p, I->k, d));
    mpz_add_ui (foo, I->c, 1);
    sign_set (&(I->polysigns), SPB,
	      sign_eval (foo, p, I->k, d));
    extra_prec = 10;
    cur_prec = mpz_sizeinbase (I->c, 2);
    do {
	/* Take the middle of I as starting point */
	    set_middle (root, I);
	if (timings != NULL)
	    time_tmp = cputime ();
	if (options.ref == NEWTON_INT)
	    nret = newton_int (root, p, pp, d, I, wp + extra_prec, &nbiter);
	else
	    nret = newton (root, p, pp, d, I, wp + extra_prec, &nbiter);
	if (timings != NULL)
	    time_tmp = cputime () - time_tmp;
	if (nret == 0)
	    break;
	cur_prec += 2; /* TODO : better? */
	dicho (I, p, d, cur_prec);
	curr_failures++;
	if (timings != NULL)
	    timings->newton_bad += time_tmp;
#if 0
	printf ("Newton s'est vautré après %d itérations.\n", nbiter);
#endif
    } while (nret != 0);
#if 0
    /* Newton converged. See if we can round the sqrt */
    do {
	exp = mpfr_get_z_exp (u, root);
	diff = mpfr_get_prec (root) - wp - extra_prec;
	mpz_fdiv_q_2exp (u, u, diff);
	sa = sign_eval (u, p, - diff - exp, d);
	mpz_add_ui (u, u, 1);
	sb = sign_eval (u, p, - diff - exp, d);
	if (sa * sb > 0) {
	    extra_prec *= 2;
	    sa = newton (root, p, pp, d, I, wp + extra_prec);
	    printf ("Help! Help! I'm being repressed! %i\n", sa);
	    continue;
#if 0
	    mpz_out_str (stdout, 2, u);
	    printf (" %ld %ld\nsa = %d, sb = %d\n", diff, exp, sa, sb);
	    printf ("diff = %ld\npika = %ld\n", diff, mpfr_get_prec (root));
	    exit (1);
	}
	mpfr_set_prec (ret, wp + extra_prec);
	mpfr_prec_round (root, wp + extra_prec, GMP_RNDD);
	mpfr_sqrt (ret, root, GMP_RNDD);
	if (mpfr_can_round (ret, wp + extra_prec - 3, GMP_RNDD, GMP_RNDN,
			    wp))
	    break;
	else
	    extra_prec *= 2;
	newton (root, p, pp, d, I, wp + extra_prec);
	} while (1);
#endif
#endif
    if (timings != NULL)
	timings->newton_good += time_tmp;
    mpfr_set (ret, root, GMP_RNDN);
    mpfr_clear (root);
    mpz_clear (foo);
    if (timings != NULL)
	timings->refinement_failures += curr_failures;
#if VERBOSE >= 3
    printf ("%d [curr_failures]\n", curr_failures);
#endif
}
#endif /* NEWTON_BENCH */

int newton_int (mpfr_t root, mpz_t *p, mpz_t *pp, unsigned long d, interval *I,
		mp_prec_t wp, int *nbiter)
{
    mpfr_t pxl, pxr, ppxl, ppxr;
    mpfr_t newxl, newxr;
    mpfr_t adjl, adjr;
    mpfr_t curxl, curxr;
    mpfr_t tmp;
    mp_prec_t cur_prec, expected_prec;

    mpfr_init2 (curxl, mpz_sizeinbase (I->c, 2) + 2);
    mpfr_init2 (curxr, mpfr_get_prec (curxl));
    mpfr_set_z (curxl, I->c, GMP_RNDN);
    mpfr_set_z (curxr, I->c, GMP_RNDN);
    mpfr_add_ui (curxr, curxr, 1, GMP_RNDN);
    mpfr_div_2ui (curxl, curxl, I->k, GMP_RNDN);
    mpfr_div_2ui (curxr, curxr, I->k, GMP_RNDN);

    mpfr_init (pxl);
    mpfr_init (pxr);
    mpfr_init (ppxl);
    mpfr_init (ppxr);
    mpfr_init (newxl);
    mpfr_init (newxr);
    mpfr_init (adjl);
    mpfr_init (adjr);
    mpfr_init (tmp);

    *nbiter = 0;
    do {
	mpfr_sub (tmp, curxr, curxl, GMP_RNDN);
	
	cur_prec = mpfr_get_exp (curxr) - mpfr_get_exp (tmp);
	expected_prec = (cur_prec + 5) * 2;
#if VERBOSE >= 3
	mpfr_dump (curxl);
	mpfr_dump (curxr);
	printf ("Precision courante: %ld\n", cur_prec);
#endif
	/* Set appropriate precision */
	mpfr_set_prec (root, expected_prec);
	mpfr_set_prec (pxl, expected_prec);
	mpfr_set_prec (pxr, expected_prec);
	mpfr_set_prec (ppxl, expected_prec);
	mpfr_set_prec (ppxr, expected_prec);
	mpfr_set_prec (adjl, expected_prec);
	mpfr_set_prec (adjr, expected_prec);
	mpfr_set_prec (newxl, expected_prec);
	mpfr_set_prec (newxr, expected_prec);

	mpfr_add (root, curxl, curxr, GMP_RNDN);
	mpfr_mul_2si (root, root, -1, GMP_RNDN);

	horner_eval_int (pxl, pxr, root, root, p, d);
#if VERBOSE >= 3
	printf ("pxl = "); mpfr_dump (pxl);
	printf ("pxr = "); mpfr_dump (pxr);
#endif
	horner_eval_int (ppxl, ppxr, curxl, curxr, pp, d - 1);
	INTERVAL(ppxl, ppxr);
#if VERBOSE >= 3
	printf ("ppxl = "); mpfr_dump (ppxl);
	printf ("ppxr = "); mpfr_dump (ppxr);
#endif
	myint_div (adjl, adjr, pxl, pxr, ppxl, ppxr);
#if VERBOSE >= 3
	printf ("adjl = "); mpfr_dump (adjl);
	printf ("adjr = "); mpfr_dump (adjr);
	INTERVAL(adjl, adjr);
#endif
	if (mpfr_nan_p (adjl)) {
	    if (mpfr_sgn (pxl) * mpfr_sgn (pxr) > 0) {
		if (mpfr_sgn (pxl) == sign_get (I->polysigns, SPB)) {
		    mpfr_set (newxl, curxl, GMP_RNDD);
		    mpfr_set (newxr, root, GMP_RNDU);
		} else {
		    mpfr_set (newxl, root, GMP_RNDD);
		    mpfr_set (newxr, curxr, GMP_RNDU);
		}
	    } else
		return 1;
	} else 
	    myint_sub (newxl, newxr, root, root, adjl, adjr);
#if VERBOSE >= 3
	printf ("root = "); mpfr_dump (root);
	printf ("newxl = "); mpfr_dump (newxl);
	printf ("newxr = "); mpfr_dump (newxr);
	INTERVAL(newxl, newxr);
#endif
	mpfr_prec_round (curxl, expected_prec, GMP_RNDD);
	mpfr_prec_round (curxr, expected_prec, GMP_RNDU);
	/* Intersect with current interval */
	if ((mpfr_cmp (newxl, curxr) > 0) || (mpfr_cmp (newxr, curxl) < 0)) {
	    printf ("Intersection vide\n");
	    mpfr_dump (newxl);
	    mpfr_dump (newxr);
	    abort ();
	}
	if (mpfr_cmp (newxl, curxl) > 0)
	    mpfr_set (curxl, newxl, GMP_RNDD);
	if (mpfr_cmp (newxr, curxr) < 0)
	    mpfr_set (curxr, newxr, GMP_RNDU);
	INTERVAL(curxl, curxr);
    } while (cur_prec < 2 * wp);
    return 0;
}

/* End of Interval arithmetic */

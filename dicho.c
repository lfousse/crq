/* dicho.c -- dichotomy functions.

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
#include <stdlib.h>
#include <mpfr.h>
#include "crq-impl.h"
#include "usp.h"
#include "utils.h"

extern void affiche_interval (interval *);
extern void affiche_poly (mpz_t *, unsigned long);

int sign_eval_int (mpfr_t x, mpfr_t *pl, mpfr_t *pr, unsigned long n)
{
    mpfr_t resl, resr;
    int ret;

    mpfr_init2 (resl, mpfr_get_prec (x));
    mpfr_init2 (resr, mpfr_get_prec (x));

    horner_eval_int2 (resl, resr, x, x, pl, pr, n);
    if (mpfr_sgn (resl) >= 0)
	ret = 1;
    else if (mpfr_sgn (resr) <= 0)
	ret = -1;
    else
	ret = 42;
    mpfr_clear (resl);
    mpfr_clear (resr);
    return ret;
}

void dicho_prep (interval *I, mpz_t *P, unsigned long n)
{
    unsigned int i;
    mpz_t foo;
#if VERBOSE >= 3
    printf ("Dans dicho_prep, l'intervalle : ");
    affiche_interval (I);
    printf ("Dans dicho_prep, le poly d'origine: ");
    affiche_poly (P, n);
#endif
    Homoth (P, - I->k, n);
#if VERBOSE >= 3
    printf ("Dans dicho_prep, apres homothetie: ");
    affiche_poly (P, n);
#endif
    X2XPC (P, I->c, n);
#if VERBOSE >= 3
    printf ("Dans dicho_prep, apres translation: ");
    affiche_poly (P, n);
#endif
    sign_set (&(I->polysigns), SPA, mpz_sgn (P[0]));
    mpz_init (foo);
    mpz_set (foo, P[0]);
    for (i = 1; i <= n; i++)
	mpz_add (foo, foo, P[i]);
    sign_set (&(I->polysigns), SPB, mpz_sgn (foo));
    mpz_clear (foo);
}

void dicho_mieux (interval *I, mpz_t *P, unsigned long n, unsigned long p,
		  int *miroir)
{
    int sa, sb;
    sa = sign_get (I->polysigns, SPA);
    sb = sign_get (I->polysigns, SPB);
#if VERBOSE >= 3
    printf ("Dicho mieux, sa = %d, sb = %d.\n", sa, sb);
#endif
    while (mpz_sizeinbase (I->c, 2) < p) {
#if VERBOSE >= 3
	printf ("D'ebut de boucle, miroir [%d], poly : ", *miroir);
	affiche_poly (P, n);
#endif
	if (*miroir == 0) {
	    /* Update poly */
	    Homoth (P, -1, n);
	    X2XP1 (P, n);
	    if (mpz_sgn (P[0]) == sb) {
		*miroir = 1;
		/* Update interval */
		mpz_mul_2exp (I->c, I->c, 1);
		I->k++;
	    } else {
		/* Update interval */
		mpz_mul_2exp (I->c, I->c, 1);
		mpz_add_ui (I->c, I->c, 1);
		I->k++;
	    }
	} else {
	    /* Update poly */
	    Homoth (P, -1, n);
	    X2XM1 (P, n);
	    if (mpz_sgn (P[0]) == sa) {
		*miroir = 0;
		/* Update interval */
		mpz_mul_2exp (I->c, I->c, 1);
		mpz_add_ui (I->c, I->c, 1);
		I->k++;
	    } else {
		/* Update interval */
		mpz_mul_2exp (I->c, I->c, 1);
		I->k++;
	    }
	}
#if VERBOSE >= 3
	affiche_poly (P, n);
#endif
    }
}
	    
void dicho (interval *I, mpz_t *P, unsigned long n, unsigned long p)
{
    int sa, sb, sm;
    mpz_t a, b, m;
    unsigned long k;

    mpz_init (m);
    mpz_init (a);
    mpz_init (b);

    mpz_set (a, I->c);
    mpz_set (b, I->c);
    mpz_add_ui (b, b, 1);
    k = I->k;

    sa = sign_get (I->polysigns, SPA);
    sb = sign_get (I->polysigns, SPB);
    
    if (sa == sb) {
	printf ("Erreur, pas de racine dans l'intervalle.\n");
	mpz_out_str (stdout, DBASE, a);
	printf ("\n%ld sa = %d, sb = %d\n", I->k, sa, sb);
	exit (1);
    }

    while (mpz_sizeinbase (a, 2) < p) {
	mpz_add (m, a, b); /* TODO : inefficace */
	k = k + 1;
	sm = sign_eval (m, P, k, n);
	if (sm == sa) {
	    mpz_swap (a, m);
	    mpz_mul_2exp (b, b, 1);
	}
	else {
	    mpz_swap (b, m);
	    mpz_mul_2exp (a, a, 1);
	}
    }
    mpz_set (I->c, a);
    I->k = k;
    mpz_clear (a);
    mpz_clear (b);
    mpz_clear (m);
}

void root_dicho (mpfr_t ret, interval *I, mpz_t *p, unsigned int d,
		 mp_prec_t wp, struct crq_stats *timings)
{
    mpz_t foo;
    mpz_init (foo);
    mp_prec_t extra_prec = 4;
#ifdef IMP_DICHO
    int miroir;
    mpz_t *pcopy;
    unsigned int i;

    pcopy = malloc ((d+1) * sizeof (*pcopy));
    for (i = 0; i <= d; i++) {
	mpz_init (pcopy[i]);
	mpz_set (pcopy[i], p[i]);
    }
    if (timings != NULL)
	timings->dich_prep_time -= cputime ();
    dicho_prep (I, pcopy, d);
    if (timings != NULL)
	timings->dich_prep_time += cputime ();
#else
    sign_set (&(I->polysigns), SPA,
	      sign_eval (I->c, p, I->k, d));
    mpz_add_ui (foo, I->c, 1);
    sign_set (&(I->polysigns), SPB,
	      sign_eval (foo, p, I->k, d));
#endif
    do {
	if (timings != NULL)
	    timings->dicho_time -= cputime ();
#ifdef IMP_DICHO
	miroir = 0;
	dicho_mieux (I, pcopy, d, wp + extra_prec, &miroir);
#else
	dicho (I, p, d, wp + extra_prec);
#endif
	if (timings != NULL)
	    timings->dicho_time += cputime ();
	mpfr_set_prec (ret, wp + extra_prec);
	mpfr_set_z (ret, I->c, GMP_RNDN);
	mpfr_mul_2si (ret, ret, - I->k, GMP_RNDN);
#if 0
	mpfr_sqrt (ret, ret, GMP_RNDD);
	if (mpfr_can_round (ret, wp + extra_prec - 2, GMP_RNDD, GMP_RNDN, wp))
	    break;
	extra_prec *= 2;
#endif
    } while (0);
    mpfr_prec_round (ret, wp, GMP_RNDN);
    mpz_clear (foo);
}

void truncate_poly (mpfr_t *pl, mpfr_t *pr, mpz_t *P, unsigned int d)
{
    unsigned int i;

    for (i = 0; i <= d; i++) {
	mpfr_set_z (pl[i], P[i], GMP_RNDD);
	mpfr_set_z (pr[i], P[i], GMP_RNDU);
    }
}

/* Dicho with truncated interval coefficients for the polynomial */
int root_dicho_trunc (mpfr_t ret, interval *I, mpz_t *p, mpfr_t *pl, mpfr_t *pr,
		      unsigned int d, mp_prec_t np)
{
    mpfr_t a, b, t;
    int st, sa, sb;
    unsigned int j;
    mp_prec_t cur_prec, wp;
    mp_prec_t extra_prec = 4;

    cur_prec = mpz_sizeinbase (I->c, 2);
    wp = mpfr_get_prec (pl[0]);
    mpfr_init2 (a, wp);
    mpfr_set_z (a, I->c, GMP_RNDN);
    mpfr_mul_2si (a, a, - I->k, GMP_RNDN);
    mpfr_init2 (b, wp);
    mpfr_set_z (b, I->c, GMP_RNDN);
    mpfr_add_ui (b, b, 1, GMP_RNDN);
    mpfr_mul_2si (b, b, - I->k, GMP_RNDN);

    sa = sign_eval_int (a, pl, pr, d);
    sb = sign_eval_int (b, pl, pr, d);
    if ((sa == 42) || (sb == 42))
	abort ();
    mpfr_init2 (t, wp + 42);

    do {
	mpfr_add (t, a, b, GMP_RNDN);
	mpfr_mul_2si (t, t, -1, GMP_RNDN);
	st = sign_eval_int (t, pl, pr, d);
	if (st == 42) {
	    extra_prec *= 2;
	    /* Truncate polynomial to a better precision */
	    for (j = 0; j <= d; j++) {
		mpfr_set_prec (pl[j], wp + extra_prec);
	    	mpfr_set_prec (pr[j], wp + extra_prec);
		mpfr_set_z (pl[j], p[j], GMP_RNDD);
		mpfr_set_z (pr[j], p[j], GMP_RNDU);
	    }
	    mpfr_set_prec (t, wp + extra_prec);
	}
	if (st * sa >= 0)
	    mpfr_set (a, t, GMP_RNDD); /* no rounding here */
	else
	    mpfr_set (b, t, GMP_RNDU); /* no rounding here */
	cur_prec++;
	wp++;
    } while (cur_prec < np);
    mpfr_set (ret, t, GMP_RNDN);
    mpfr_prec_round (ret, np, GMP_RNDN);
    mpfr_clear (a);
    mpfr_clear (b);
    mpfr_clear (t);
    return 0;
}

/* secante_root.c -- Secante iteration root refinement.

Copyright 2005,2006 Laurent Fousse <laurent@komite.net>

This file is part of the CRQ Library.

The CRQ Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 2.1
as published by the Free Software Foundation.

The CRQ Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include "crq-impl.h"
#include <assert.h>
#include <stdio.h>
#include <mpfr.h>
#ifndef NDEBUG
#include "utils.h"
#endif

/* TODO : this function might be added to MPFR
 * in some form or another */
int mpfr_eqbits (mpfr_t a, mpfr_t b)
{
    int minb, maxb, m;

    minb = 0;
    maxb = mpfr_get_prec (a);

    while (maxb - minb > 1) {
	m = (maxb + minb) / 2;
	if (mpfr_eq (a, b, m))
	    minb = m;
	else
	    maxb = m;
    }
    return minb;
}

int secante (struct crq_gl_root root, mpz_t *p, unsigned long deg,
	      mp_prec_t wp)
{
    mp_prec_t guarde;
    mp_prec_t prec[2];
    long C, t;
    int signs[3];
    int soffset;
    int boucles;

    mpfr_t x[2], f[2], d, n, herr, zero, bleft, bright;
    int par, smin, smax;

    assert (mpfr_cmp (*root.xmin, *root.xmax) < 0);

    prec[0] = mpfr_get_prec (*root.xmin);
    prec[1] = mpfr_get_prec (*root.xmax);
    C = 0;
    guarde = 1;

    mpfr_init2 (x[0], prec[0] + guarde);
    mpfr_init2 (x[1], prec[1] + guarde);
    mpfr_init2 (bleft, prec[0]);
    mpfr_set (bleft, *root.xmin, GMP_RNDN);
    mpfr_init2 (bright, prec[1]);
    mpfr_set (bright, *root.xmax, GMP_RNDN);
    mpfr_set (x[0], *root.xmin, GMP_RNDN);
    mpfr_set (x[1], *root.xmax, GMP_RNDN);
    mpfr_init2 (f[0], prec[0] + guarde);
    mpfr_init2 (f[1], prec[1] + guarde);
    mpfr_init (d);
    mpfr_init (n);
    mpfr_init (zero);
    mpfr_init (herr);

    mpfr_set_ui (zero, 0, GMP_RNDN);

    assert (mpfr_cmp (x[0], x[1]) != 0);
#if VERBOSE >= 3
    printf ("\nxmin = ");
    mpfr_out_str (stdout, DBASE, 0, *root.xmin, GMP_RNDN);
    printf ("\nxmax = ");
    mpfr_out_str (stdout, DBASE, 0, *root.xmax, GMP_RNDN);
#endif
#ifndef NDEBUG
    /* Check that sign(f(xmin)).sign(f(xmax)) <= 0 */
    {
	int smin, smax;
	if (mpfr_zero_p (*root.xmin))
	    smin = mpz_sgn (p[0]);
	else {
	    mpz_t z;
	    mp_exp_t k;
	    mpz_init (z);
	    k = mpfr_get_z_exp (z, *root.xmin);
	    smin = sign_eval (z, p, -k, deg);
	    mpz_clear (z);
	}
	if (mpfr_zero_p (*root.xmax))
	    smax = mpz_sgn (p[0]);
	else {
	    mpz_t z;
	    mp_exp_t k;
	    mpz_init (z);
	    k = mpfr_get_z_exp (z, *root.xmax);
	    smax = sign_eval (z, p, -k, deg);
	    mpz_clear (z);
	}
	assert (smax * smin <= 0);
    }
#endif
    if (prec[0] > prec[1])
	par = 0;
    else
	par = 1;

    do {
	horner_eval (f[1 - par], x[1 - par], p, deg, zero, herr);
	if (mpfr_cmpabs (f[1 - par], herr) > 0)
	    break;
	guarde *= 2;
	mpfr_set_prec (f[1 - par], prec[1 - par] + guarde);
    } while (1);
    
    if (par == 0) {
	smax = mpfr_sgn (f[1 - par]);
	smin = - smax;
    } else {
	smin = mpfr_sgn (f[1 - par]);
	smax = - smin;
    }
    
    boucles = 0;
    soffset = 0;
    for (;;) {
#if VERBOSE >= 3
	printf ("<secante>\n");
	printf ("x_{n+1} = ");
	mpfr_dump (x[par]);
	printf ("prec[n+1] = %ld\n", prec[par]);
	printf ("wp = %ld\n", wp);
	printf ("prec[n] = %ld\n", prec[1 - par]);
	printf ("C = %ld\n", C);
#endif
	t = prec[par] + C + guarde;
	if (t < 2)
	    t = 53;
	mpfr_set_prec (f[par], t);
	mpfr_set_prec (d, t);
	mpfr_set_prec (n, t);

	do {
	    horner_eval (f[par], x[par], p, deg, zero, herr);
	    if (mpfr_cmpabs (f[par], herr) > 0)
		break;
	    t += guarde;
	    guarde *= 2;
	    mpfr_set_prec (f[par], t);
	} while (1);
	signs[soffset] = mpfr_sgn (f[par]);
	if (signs[soffset] * smin > 0) {
	    if (mpfr_cmp (*root.xmin, x[par]) < 0) {
		mpfr_set_prec (*root.xmin, mpfr_get_prec (x[par]));
		mpfr_set (*root.xmin, x[par], GMP_RNDN);
#if VERBOSE >= 3
		printf ("Update xmin.\n");
#endif
		assert (mpfr_cmp (*root.xmin, *root.xmax) < 0);
	    }
	} else {
	    if (mpfr_cmp (*root.xmax, x[par]) > 0) {
		mpfr_set_prec (*root.xmax, mpfr_get_prec (x[par]));
		mpfr_set (*root.xmax, x[par], GMP_RNDN);
#if VERBOSE >= 3
		printf ("Update xmax.\n");
#endif
		assert (mpfr_cmp (*root.xmin, *root.xmax) < 0);
	    }
	}
	/* See if we got three times the same sign in a row */
	if ((boucles > 2) && (signs[0] == signs[1]) && (signs[1] == signs[2])) {
#if VERBOSE >= 3
	    printf ("signs[:] uniformes\n");
#endif
	    return 0;
	}
	mpfr_sub (n, x[par], x[1 - par], GMP_RNDN);
	mpfr_sub (d, f[par], f[1 - par], GMP_RNDN);
	if (mpfr_zero_p (d)) {
#if VERBOSE >= 3
	    printf ("f[n+1] - f[n] == 0\n");
#endif
	    return 0;
	}
	if (mpfr_zero_p (n)) {
#if VERBOSE >= 3
	    printf ("x[n+1] - x[n] == 0\n");
#endif
	    return 0;
	}
	mpfr_div (n, n, d, GMP_RNDN);
	mpfr_mul (n, n, f[par], GMP_RNDN);
#if VERBOSE >= 3
	printf ("r = ");
	mpfr_dump (n);
#endif
	t = mpfr_get_exp (x[par]) - mpfr_get_exp (n) - prec[par];
#if VERBOSE >= 3
	printf ("c = %ld\n", t);
#endif
	t = 0.9 * C + 0.1 * (C + t);
	C = t;
#if VERBOSE >= 3
	printf ("C' = %ld\n", C);
#endif
	t = prec[par] + t;
	if (t < 2)
	    t = 53;
	prec[par] = t;
#if VERBOSE >= 3
	printf ("prec[n+1]' = %ld\n", prec[par]);
#endif
	t = prec[1 - par] + prec[par] + C;
	if (t < 2)
	    t = 53;
	prec[1 - par] = t;
#if VERBOSE >= 3
	printf ("prec[n+2] = %ld\n", prec[1 - par]);
#endif
	mpfr_prec_round (x[1 - par], prec[1 - par], GMP_RNDN);
#if VERBOSE >= 3
	printf ("E(r) = %ld\n", mpfr_get_exp (n));
#endif
	mpfr_sub (x[1 - par], x[par], n, GMP_RNDN);
	if (mpfr_cmp (x[1 - par], bleft) < 0) {
#if VERBOSE >= 3
	    printf ("sortie de piste par la gauche\n");
#endif
	    return 0;
	}
	if (mpfr_cmp (x[1 - par], bright) > 0) {
#if VERBOSE >= 3
	    printf ("sortie de piste par la droite\n");
#endif
	    return 0;
	}

	par = 1 - par;
#if VERBOSE >= 3
	printf ("</secante>\n");
#endif
	boucles++;
	soffset = (soffset + 1) % 3;
	mpfr_sub (d, *root.xmax, *root.xmin, GMP_RNDN);
	if (mpfr_zero_p (d)) {
#if VERBOSE >= 3
	    printf ("xmin == xmax.\n");
#endif
	    return (mpfr_get_prec (*root.xmax) >= wp);
	}
#if VERBOSE >= 3
	printf ("exp(xmax) = %ld, wp = %ld, exp(d) = %ld\n",
		mpfr_get_exp (*root.xmax), wp, mpfr_get_exp (d));
#endif
	if (mpfr_get_exp (*root.xmax) - ((signed long) wp) > mpfr_get_exp (d))
	    return 1;
    }
}

int safe_horner_sign (mpz_t *p, unsigned long d, mpfr_t x, unsigned long base,
		      unsigned long *guarde)
{
    mpfr_t zero, herr, fx;
    int ret;

    mpfr_init (herr);
    mpfr_init (zero);
    mpfr_init2 (fx, base + *guarde);
    mpfr_set_ui (zero, 0, GMP_RNDN);

    do {
#if VERBOSE >= 3
	printf ("pon\n");
#endif
	horner_eval (fx, x, p, d, zero, herr);
#if VERBOSE >= 3
	mpfr_dump (fx);
#endif
#if VERBOSE >= 3
	mpfr_dump (herr);
#endif
	if (mpfr_cmpabs (fx, herr) > 0)
	    break;
	*guarde *= 2;
	mpfr_set_prec (fx, base + *guarde);
    } while (1);
    ret = mpfr_sgn (fx);
    mpfr_clear (fx);
#ifndef NDEBUG
    /* How safe is safe_horner_sign ? */
    {
	if (mpfr_zero_p (x)) {
	    assert (mpz_sgn (p[0]) == ret);
	    return ret;
	}
	mpz_t z;
	mp_exp_t k;
	int proper_sign;
	mpz_init (z);
	k = mpfr_get_z_exp (z, x);
	proper_sign = sign_eval (z, p, -k, d);
	assert (proper_sign == ret);
	mpz_clear (z);
    }
#endif
    return ret;
}

void root_secante (mpfr_t root, interval *I, mpz_t *p, unsigned long d,
		   mp_prec_t wp, struct crq_gl_opts options, struct crq_stats *timings)
{
    struct crq_gl_root inter;
    int ret;
    int t;
    mp_prec_t sp;
    unsigned long guarde;
    int smin, smoy;

    mpfr_t xmin, xmax, fm, zero, herr, moy;

    sp = mpz_sizeinbase (I->c, 2);
    if (sp <= 53) {
	sp = 53;
    }
    mpfr_init2 (xmin, sp);
    mpfr_init2 (xmax, sp);
    mpfr_init2 (moy, sp);
    mpfr_init2 (fm, sp);
    mpfr_init (zero);
    mpfr_init (herr);
    guarde = 1;

    mpfr_set_ui (zero, 0, GMP_RNDN);

    mpfr_set_z (xmin, I->c, GMP_RNDN);
    mpfr_div_2ui (xmin, xmin, I->k, GMP_RNDN);
    mpfr_set_z (xmax, I->c, GMP_RNDN);
    mpfr_add_ui (xmax, xmax, 1, GMP_RNDN);
    mpfr_div_2ui (xmax, xmax, I->k, GMP_RNDN);

    assert (mpfr_cmp (xmin, xmax) < 0);
#if VERBOSE >= 3
    mpfr_dump (xmin);
    mpfr_dump (xmax);
#endif
    inter.xmin = &xmin;
    inter.xmax = &xmax;
    inter.approx = &moy;
    smin = safe_horner_sign (p, d, xmin, sp, &guarde);

    do {
	ret = secante (inter, p, d, wp);
	if (ret)
	    break;
	assert (mpfr_cmp (xmin, xmax) < 0);

#ifndef NDEBUG
	/* Check that sign(f(xmin)).sign(f(xmax)) <= 0 */
	{
	    int smin, smax;
	    if (mpfr_zero_p (xmin))
		smin = mpz_sgn (p[0]);
	    else {
		mpz_t z;
		mp_exp_t k;
		mpz_init (z);
		k = mpfr_get_z_exp (z, xmin);
		smin = sign_eval (z, p, -k, d);
		mpz_clear (z);
	    }
	    if (mpfr_zero_p (xmax))
		smax = mpz_sgn (p[0]);
	    else {
		mpz_t z;
		mp_exp_t k;
		mpz_init (z);
		k = mpfr_get_z_exp (z, xmax);
		smax = sign_eval (z, p, -k, d);
		mpz_clear (z);
	    }
	    assert (smax * smin <= 0);
	}
#endif
	if (mpfr_get_prec (*inter.xmin) > mpfr_get_prec (*inter.xmax))
	    t = mpfr_get_prec (*inter.xmin);
	else
	    t = mpfr_get_prec (*inter.xmax);
#if VERBOSE >= 3
	printf ("t = %d\n", t);
	printf ("prec(xmin) = %d, prec(xmax) = %d\n", mpfr_get_prec (*inter.xmin), mpfr_get_prec (*inter.xmax));
#endif
	mpfr_set_prec (*inter.approx, t + 1);
	mpfr_add (*inter.approx, *inter.xmin, *inter.xmax, GMP_RNDN);
	mpfr_div_2ui (*inter.approx, *inter.approx, 1, GMP_RNDN);
	smoy = safe_horner_sign (p, d, *inter.approx, sp, &guarde);
#if VERBOSE >= 3
	printf ("moy = ");
	mpfr_out_str (stdout, DBASE, 0, *inter.approx, GMP_RNDN);
	printf ("\nsmoy = %d\n", smoy);
	printf ("\nsmin = %d\n", smin);
#endif
	if (smoy * smin > 0) {
	    mpfr_prec_round (*inter.xmin, t + 1, GMP_RNDN);
	    mpfr_set (*inter.xmin, *inter.approx, GMP_RNDN);
	} else {
	    mpfr_prec_round (*inter.xmax, t + 1, GMP_RNDN);
	    mpfr_set (*inter.xmax, *inter.approx, GMP_RNDN);
	}
    } while (! ret);
    /* convert [xmin, xmax] to an interval of the
     * form [c/2^k, (c+1)/2^k] */
    t = mpfr_eqbits (*inter.xmin, *inter.xmax);
    mpfr_set_prec (moy, t);
    mpfr_set (moy, *inter.xmin, GMP_RNDD);
    t = mpfr_get_z_exp (I->c, moy);
    I->k = - t;
    mpfr_set (root, *inter.xmin, GMP_RNDN);
}

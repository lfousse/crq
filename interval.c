/* interval.c -- custom interval methods.

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

/* Interval arithmetic */

void myint_add (mpfr_t retl, mpfr_t retr, mpfr_t al, mpfr_t ar, mpfr_t bl,
	       mpfr_t br)
{
    mpfr_add (retl, al, bl, GMP_RNDD);
    mpfr_add (retr, ar, br, GMP_RNDU);
}

void myint_sub (mpfr_t retl, mpfr_t retr, mpfr_t al, mpfr_t ar, mpfr_t bl,
	       mpfr_t br)
{
    int i;
    i = mpfr_sub (retl, al, br, GMP_RNDD);
    i = mpfr_sub (retr, ar, bl, GMP_RNDU);
}

void myint_mul (mpfr_t retl, mpfr_t retr, mpfr_t al, mpfr_t ar, mpfr_t bl,
	       mpfr_t br, mpfr_t *tmp)
{
    if (mpfr_sgn (al) >= 0) {
	if (mpfr_sgn (bl) >= 0) {
	    mpfr_mul (retl, al, bl, GMP_RNDD);
	    mpfr_mul (retr, ar, br, GMP_RNDU);
	    return;
	}
	else {
	    mpfr_mul (retl, ar, bl, GMP_RNDD);
	    if (mpfr_sgn (br) >= 0)
		mpfr_mul (retr, br, ar, GMP_RNDU);
	    else
		mpfr_mul (retr, br, al, GMP_RNDU);
	    return;
	}
    }
    else {
	if (mpfr_sgn (bl) >= 0) {
	    mpfr_mul (retl, br, al, GMP_RNDD);
	    if (mpfr_sgn (ar) >= 0)
		mpfr_mul (retr, ar, br, GMP_RNDU);
	    else
		mpfr_mul (retr, ar, bl, GMP_RNDU);
	    return;
	}
	else { /* [a] and [b] contain 0 */
	    mpfr_mul (tmp[0], ar, br, GMP_RNDU);
	    mpfr_mul (tmp[1], al, bl, GMP_RNDU);
	    if (mpfr_cmp (tmp[0], tmp[1]) > 0)
		mpfr_set (retr, tmp[0], GMP_RNDU);
	    else
		mpfr_set (retr, tmp[1], GMP_RNDU);
	    mpfr_mul (tmp[0], ar, bl, GMP_RNDD);
	    mpfr_mul (tmp[1], al, br, GMP_RNDD);
	    if (mpfr_cmp (tmp[0], tmp[1]) < 0)
		mpfr_set (retl, tmp[0], GMP_RNDD);
	    else
		mpfr_set (retl, tmp[1], GMP_RNDD);
	    return;
	}
    }
}

void myint_div (mpfr_t retl, mpfr_t retr, mpfr_t al, mpfr_t ar, mpfr_t bl,
	       mpfr_t br)
{
    /* mpfr_dump (al);
    mpfr_dump (ar);
    mpfr_dump (bl);
    mpfr_dump (br); */
    /* NaN case */
    if (mpfr_sgn (bl) * mpfr_sgn (br) <= 0) {
	mpfr_set_nan (retl);
	mpfr_set_nan (retr);
	return;
    }
    if (mpfr_sgn (al) * mpfr_sgn (ar) <= 0) {
	/* 0 \in [al, ar] */
	if (mpfr_sgn (bl) >= 0) {
	    mpfr_div (retl, al, br, GMP_RNDD);
	    mpfr_div (retr, ar, bl, GMP_RNDU);
	    return;
	} else {
	    mpfr_div (retl, ar, bl, GMP_RNDD);
	    mpfr_div (retr, al, br, GMP_RNDU);
	    return;
	}
    } else {
	if (mpfr_sgn (al) >= 0) {
	    if (mpfr_sgn (bl) >= 0) {
		mpfr_div (retl, al, br, GMP_RNDD);
		mpfr_div (retr, ar, bl, GMP_RNDU);
	    } else {
		mpfr_div (retl, ar, br, GMP_RNDD);
		mpfr_div (retr, al, bl, GMP_RNDU);
	    }
	} else {
	    if (mpfr_sgn (br) >= 0) {
		mpfr_div (retl, al, bl, GMP_RNDD);
		mpfr_div (retr, ar, br, GMP_RNDU);
	    } else {
		mpfr_div (retl, ar, bl, GMP_RNDD);
		mpfr_div (retr, al, br, GMP_RNDU);
	    }
	}
    }
}

void myint_add_z (mpfr_t resl, mpfr_t resr, mpfr_t xl, mpfr_t xr, mpz_t b)
{
    mpfr_add_z (resl, xl, b, GMP_RNDD);
    mpfr_add_z (resr, xr, b, GMP_RNDU);
}

void horner_eval_int (mpfr_t resl, mpfr_t resr, mpfr_t xl, mpfr_t xr, mpz_t *p,
		      unsigned long n)
{
    mpfr_t foo[2];
    mpfr_set_z (resl, p[n], GMP_RNDD);
    mpfr_set_z (resr, p[n], GMP_RNDU);
    mpfr_init2 (foo[0], mpfr_get_prec (resl));
    mpfr_init2 (foo[1], mpfr_get_prec (resl));

    if (n == 0)
	return;

    do {
	n--;
	myint_mul (resl, resr, resl, resr, xl, xr, foo);
	INTERVAL(resl, resr);
	myint_add_z (resl, resr, resl, resr, p[n]);
	INTERVAL(resl, resr);
    } while (n > 0);
    mpfr_clear (foo[0]);
    mpfr_clear (foo[1]);
}

void horner_eval_int2 (mpfr_t resl, mpfr_t resr, mpfr_t xl, mpfr_t xr,
		       mpfr_t *pl, mpfr_t *pr, unsigned long n)
{
    mpfr_t foo[2];
    mpfr_set (resl, pl[n], GMP_RNDD);
    mpfr_set (resr, pr[n], GMP_RNDU);
    mpfr_init2 (foo[0], mpfr_get_prec (resl));
    mpfr_init2 (foo[1], mpfr_get_prec (resl));

    if (n == 0) {
	mpfr_clear (foo[0]);
	mpfr_clear (foo[1]);
	return;
    }

    do {
	n--;
	myint_mul (resl, resr, resl, resr, xl, xr, foo);
	INTERVAL(resl, resr);
	myint_add (resl, resr, resl, resr, pl[n], pr[n]);
	INTERVAL(resl, resr);
    } while (n > 0);
    mpfr_clear (foo[0]);
    mpfr_clear (foo[1]);
}

void intersect (mpfr_t resl, mpfr_t resr, mpfr_t al, mpfr_t ar, mpfr_t bl,
		mpfr_t br)
{
    if (mpfr_cmp (al, bl) < 0)
	mpfr_set (resl, bl, GMP_RNDD);
    else
	mpfr_set (resl, al, GMP_RNDD);
    if (mpfr_cmp (ar, br) < 0)
	mpfr_set (resr, ar, GMP_RNDD);
    else
	mpfr_set (resr, br, GMP_RNDD);
}

int is_in_interval (mpfr_t root, interval *I)
{
    mpfr_t a, b;
    mp_prec_t p;

    p = mpz_sizeinbase (I->c, 2) + 3;
    mpfr_init2 (a, p);
    mpfr_init2 (b, p);
    mpfr_set_z (a, I->c, GMP_RNDN);
    mpfr_set (b, a, GMP_RNDN);
    mpfr_add_ui (b, b, 1, GMP_RNDN);
    mpfr_mul_2si (a, a, - I->k, GMP_RNDN);
    mpfr_mul_2si (b, b, - I->k, GMP_RNDN);

    if (mpfr_cmp (root, a) < 0) {
#if VERBOSE >= 3
	printf ("Dans is_in_interval, trop petit:\n");
	mpfr_out_str (stdout, DBASE, 0, a, GMP_RNDN);
	printf (" [a]\n");
	mpfr_out_str (stdout, DBASE, 0, root, GMP_RNDN);
	printf (" [r00t]\n");
	mpfr_out_str (stdout, DBASE, 0, b, GMP_RNDN);
	printf (" [b]\n");
	mpz_out_str (stdout, DBASE, I->c);
	printf (" [I->c]\n%ld\n", I->k);
#endif
	mpfr_clear (a);
	mpfr_clear (b);
	return 0;
    }

    if (mpfr_cmp (root, b) > 0) {
#if VERBOSE >= 3
	printf ("Dans is_in_interval, trop grand:\n");
	mpfr_out_str (stdout, DBASE, 0, a, GMP_RNDN);
	printf (" [a]\n");
	mpfr_out_str (stdout, DBASE, 0, root, GMP_RNDN);
	printf (" [r00t]\n");
	mpfr_out_str (stdout, DBASE, 0, b, GMP_RNDN);
	printf (" [b]\n");
#endif
	mpfr_clear (a);
	mpfr_clear (b);
	return 0;
    }
    mpfr_clear (a);
    mpfr_clear (b);
    return 1;
}


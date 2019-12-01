/* nc_gen.c -- generic method for Newton-Cotes.

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
#include <string.h>
#include <gmp.h>
#include "crq.h"
#include "crq-impl.h"
#include <assert.h>
#include <mpfr.h>

#include "utils.h"
#include "listz.h"

#if 0
#define NDEBUG
#endif

void
list_mul (listz_t a, listz_t b, unsigned int k, int monic_b,
	          listz_t c, unsigned int l, int monic_c, listz_t t);

void ProductTree (listz_t *tree, unsigned int n0, unsigned int n1, listz_t T, mpz_t n);
int list_mul_mem (unsigned int len);
int polyeval_tellegen (listz_t b, unsigned int k, listz_t *Tree, listz_t tmp,
                   unsigned int sizeT, listz_t invF, mpz_t n);
void frac_rec (mpz_t a, mpz_t b, mpz_t u, mpz_t n, mpz_t B);

#define FUNCTION exp_wrapper
#define LEFT_BOUND 0
#define RIGHT_BOUND 3
/* #define DEBUG_EXACT */
/* #define TEST2 */
/*#define DTIS*/

#define ALLOW_CHEAT

#ifdef DEBUG_EXACT

void mpfr_add_exact (mpfr_t ret, mpfr_t a, mpfr_t b)
{
    mp_prec_t needed;

    if (mpfr_nan_p (a) || mpfr_nan_p (b) || mpfr_nan_p (ret)) {
	printf ("Duh. NAN.\n");
	abort ();
    }
    if mpfr_zero_p (a) {
	mpfr_set_prec (ret, mpfr_get_prec (b));
	mpfr_set (ret, b, GMP_RNDN);
	return;
    }
    if mpfr_zero_p (b) {
	mpfr_set_prec (ret, mpfr_get_prec (a));
	mpfr_set (ret, a, GMP_RNDN);
	return;
    }

    /* Grossly overestimate needed precision */
    if (mpfr_get_exp (a) > mpfr_get_exp (b)) {
	needed = mpfr_get_exp (a) - mpfr_get_exp (b);
	if (mpfr_get_prec (a) > mpfr_get_prec (b))
	    needed += mpfr_get_prec (a);
	else
	    needed += mpfr_get_prec (b);
    }
    else {
	needed = mpfr_get_exp (b) - mpfr_get_exp (a);
	if (mpfr_get_prec (a) > mpfr_get_prec (b))
	    needed += mpfr_get_prec (a);
	else
	    needed += mpfr_get_prec (b);
    }
    needed++;

    mpfr_prec_round (ret, needed, GMP_RNDN);
    mpfr_add (ret, a, b, GMP_RNDN);
    if (mpfr_nan_p (ret)) {
	printf ("There is a pestilence upon this land. Nothing is sacred\n");
	abort ();
    }
}

#endif

/* unused */
#if 0
void mpfr_sum_no_loop (mpfr_t res, mpfr_ptr tab, unsigned int n) 
{
    mpfr_srcptr *perm;

    perm = malloc (n * sizeof (*perm));
    mpfr_count_sort (tab, n, perm);

    free (perm);
}
#endif


/* assume b[0]...b[n/2] and d are initialized
   Return b[] and d such that b[k]/d is the k-th Newton-Cotes coefficient
   of order n, i.e. (-1)^(n-k)/k!/(n-k)!*int(t*(t-1)*...*(t-n)/(t-k),t=0..n).
   
   Remark: we don't need b[n/2+1]...b[n] since b[k]=b[n-k].
*/
void compute_bnk (mpz_t * b, mpz_t d, unsigned long n)
{
    unsigned long k, g, j;
    mpz_t num, den, pow_n, *p, fac_n;

    mpz_init (num);
    mpz_init (den);
    mpz_init (pow_n);
    mpz_init (fac_n);

    /* compute (multiple of) common denominator */
    mpz_fac_ui (fac_n, n);
    mpz_set_ui (num, 1);
    for (k = 2; k <= n; k++)
	mpz_lcm_ui (num, num, k);
    mpz_divexact (num, fac_n, num);	/* n! / lcm(1, 2, ..., n) */
    mpz_mul_ui (num, num, n);
    mpz_set_ui (pow_n, n);
    mpz_set_ui (den, 1);
    for (k = 1; k <= n; k++) {
	mpz_mul_ui (pow_n, pow_n, n);	/* pow_n = n^(k+1) */
	g = mpz_gcd_ui (NULL, pow_n, k + 1);
	mpz_divexact_ui (b[0], pow_n, g);
	mpz_gcd (num, num, b[0]);
	mpz_lcm_ui (den, den, (k + 1) / g);
    }
    mpz_gcd (num, num, den);
    mpz_divexact (den, den, num);

    p = (mpz_t *) malloc ((n + 1) * sizeof (mpz_t));
    for (k = 0; k <= n; k++)
	mpz_init (p[k]);
    /* put coefficients of (t-1)*(t-2)*...*(t-n) in p */
    mpz_set_ui (p[0], 1);
    for (k = 1; k <= n; k++) {
	/* p[0..k-1] contains the coefficients of (t-1)*...*(t-(k-1)) */
	mpz_set_ui (p[k], 0);
	for (j = k; j > 0; j--) {
	    mpz_mul_ui (p[j], p[j], k);
	    mpz_sub (p[j], p[j - 1], p[j]);
	}
	mpz_mul_si (p[0], p[0], -k);
    }
    /* compute b[0] */
    mpz_set (pow_n, den);
    mpz_set_ui (b[0], 0);
    for (j = 0; j <= n; j++) {
	mpz_mul_ui (pow_n, pow_n, n);	/* pow_n = den * n^(j+1) */
	mpz_mul (num, p[j], pow_n);
	mpz_divexact_ui (num, num, j + 1);
	mpz_add (b[0], b[0], num);
    }
    if (n % 2)
	mpz_neg (b[0], b[0]);

    /* compute b[1]...b[n/2] */
    mpz_set_ui (d, 1);		/* d = n!/k!/(n-k)! */
    for (k = 1; k <= n / 2; k++) {
	mpz_mul_ui (d, d, n - (k - 1));
	mpz_divexact_ui (d, d, k);
	/* multiply by t-k */
	for (j = n - 1; j > 0; j--)
	    mpz_addmul_ui (p[j], p[j + 1], k);
	/* divide by t-(k-1) */
	mpz_set_ui (p[0], 0);
	for (j = 0; j < n; j++)
	    mpz_submul_ui (p[j], p[j + 1], k - 1);
	/* set b[k] */
	mpz_set (pow_n, den);
	mpz_set_ui (b[k], 0);
	for (j = 0; j <= n; j++) {
	    mpz_mul_ui (pow_n, pow_n, n);	/* pow_n = den * n^(j+1) */
	    mpz_mul (num, p[j], pow_n);
	    mpz_divexact_ui (num, num, j + 1);
	    mpz_add (b[k], b[k], num);
	}
	if ((n - k) % 2)
	    mpz_neg (b[k], b[k]);
	mpz_mul (b[k], b[k], d);
    }

    for (k = 0; k <= n; k++)
	mpz_clear (p[k]);
    free (p);

    mpz_mul (d, fac_n, den);

#if 0
    /* compute smallest denominator */
    mpz_gcd (den, d, b[0]);
    for (k = 1; k <= n / 2; k++)
	mpz_gcd (den, den, b[k]);
    mpz_divexact (d, d, den);
    for (k = 0; k <= n / 2; k++)
	mpz_divexact (b[k], b[k], den);
#endif

    mpz_clear (num);
    mpz_clear (den);
    mpz_clear (pow_n);
    mpz_clear (fac_n);
}

void compute_bnk2 (mpz_t * b, mpz_t d, unsigned long n)
{
    listz_t *tree;
    listz_t T, S, R, invTop;
    int i, logn;
    unsigned long k, g;
    int c;
    mpz_t P, prodtmp, inv, zi, adj;
    mpfr_t acc;
    mpz_t *denum, B;
    mpfr_t u, lognf;
    unsigned long sqrtnp1;
    int tempspace = 100 * list_mul_mem (2 * n);
    int time_adj, time_bound, time_calcinv, time_calcp, time_calcr;
    int time_calcs, time_lcm, time_polyeval, time_ptree, time_refrac;

    mpz_t num, den, pow_n, fac_n;
#ifdef ALLOW_CHEAT
    int optnums[] = { 4, 4, 15 /* 8 */, 34, 133, 279, 831, 1606 /* 256 */, 4143, 8765, 21105};
    int optdenums[] = { 4, 3, 15 /* 8 */, 30, 119, 232, 727, 1387 /* 256 */, 3674, 7796, 19135 };
#endif

    mpz_init (num);
    mpz_init (den);
    mpz_init (pow_n);
    mpz_init (fac_n);


    T = malloc (tempspace * sizeof (*T));
    for (i = 0; i < tempspace; i++) {
	mpz_init (T[i]);
    }

    time_bound = cputime ();
    /* Chose P prime with 2, ... n and big enough
     * to recover the coefficients */
    mpz_init (P);
    mpfr_init (u);
    mpfr_init (acc);
    mpfr_init (lognf);
    mpfr_set_ui (lognf, 0, GMP_RNDN);
    mpfr_set_ui (acc, 0, GMP_RNDN);
    /* Compute an upper bound of log2 ((n-1)!) */
    for (i = 2; i < n; i++) {
	mpfr_set_ui (u, i, GMP_RNDU);
	mpfr_log2 (u, u, GMP_RNDU);
	mpfr_add (lognf, lognf, u, GMP_RNDU);
    }
    mpfr_set_ui (u, n + 1, GMP_RNDN);
    mpfr_sqrt (u, u, GMP_RNDN);
    sqrtnp1 = mpfr_get_ui (u, GMP_RNDU);
    for (i = 2; i < sqrtnp1 - 1; i++) {
	mpfr_set_ui (u, i, GMP_RNDD);
	mpfr_log2 (u, u, GMP_RNDD);
	mpfr_mul_ui (u, u, 3, GMP_RNDD);
	mpfr_sub (acc, acc, u, GMP_RNDU);
    }
    mpfr_set_ui (u, sqrtnp1 - 1, GMP_RNDD);
    mpfr_log2 (u, u, GMP_RNDD);
    mpfr_mul_2ui (u, u, 1, GMP_RNDD);
    mpfr_sub (acc, acc, u, GMP_RNDU);
    for (i = sqrtnp1; i <= n + 1 - sqrtnp1; i++) {
	mpfr_set_ui (u, i, GMP_RNDD);
	mpfr_log2 (u, u, GMP_RNDD);
	mpfr_sub (acc, acc, u, GMP_RNDU);
    }
    mpfr_set_ui (u, n - 1, GMP_RNDU);
    mpfr_log2 (u, u, GMP_RNDU);
    mpfr_mul_ui (u, u, n - 1, GMP_RNDU);
    mpfr_add (acc, acc, u, GMP_RNDU);
    mpfr_set_ui (u, n, GMP_RNDU);
    mpfr_log2 (u, u, GMP_RNDU);
    mpfr_add (acc, acc, u, GMP_RNDU);
    mpfr_mul_2ui (u, lognf, 1, GMP_RNDU);
    /* u is now a (rough) upper bound of log2 (B_denum) */
    mpfr_add (acc, acc, u, GMP_RNDU);
    /* acc is now an upper bound of log2 (B_num) */
#ifdef ALLOW_CHEAT
    {
	int U = n;
	int log2n = 0;
	while (U > 1) {
	    log2n++;
	    U = U / 2;
	}
	printf ("log2n = %d\n", log2n);
	mpfr_set_ui (acc, optnums[log2n - 1], GMP_RNDU);
	mpfr_set_ui (u, optdenums[log2n - 1], GMP_RNDU);
    }
#endif
    printf ("%ld [log2 num]\n%ld [log2 denum]\n", mpfr_get_ui (acc, GMP_RNDU), mpfr_get_ui (u, GMP_RNDU));
    time_bound = cputime () - time_bound;
    time_calcp = cputime ();
    mpz_init (B);
    mpz_set_ui (B, 1);
    mpz_mul_2exp (B, B, mpfr_get_ui (acc, GMP_RNDU));

    mpz_set_ui (P, 1);
    mpz_mul_2exp (P, P, mpfr_get_ui (acc, GMP_RNDU) + mpfr_get_ui (u, GMP_RNDU));
    c = 1;
    mpz_add_ui (P, P, 1);
    for (;;) {
check:
	for (i = 3; i <= n; i++) {
	    if (mpz_fdiv_ui (P, i) == 0) {
		c += 2;
		mpz_add_ui (P, P, 2);
		goto check;
	    }
	}
	break;
    }
    time_calcp = cputime () - time_calcp;
#ifndef NDEBUG
    printf ("P = ");
    mpz_out_str (stdout, 10, P);
    printf ("\n");
#endif

    time_ptree = cputime ();
    /* Allocate memory for the tree */
    logn = 0;
    while ((n >> logn) > 0)
	logn++;
    logn++;
    tree = malloc (logn * sizeof (*tree));
    for (i = 0; i < logn; i++) {
	tree[i] = malloc ((n+1) * sizeof (mpz_t));
	/* +1 for root */
    }
    ProductTree (tree, 0, n - 1, T, P);
    time_ptree = cputime () - time_ptree;
#ifndef NDEBUG
    printf ("Top = ");
    pari_poly (tree[0], n, 1);
    printf ("\n");
#endif
    /* tree[0] is now x(x-1)...(x-(n-1)), of degree n, leading
     * coefficient implicit */

    /* S of degree n-1, with coefficients s_i = (n-1)^(n-i) / (n-i) mod P..
     * not monic, so no implicit coefficient. */
    time_calcs = cputime ();
    S = malloc (n * sizeof (*S));
    mpz_init (prodtmp);
    mpz_init (inv);
    mpz_init (zi);
    mpz_set_ui (prodtmp, n - 1);
    for (i = 0; i < n; i++) {
	mpz_init (S[n - 1 - i]);
	mpz_set (S[n - 1 - i], prodtmp);
	mpz_mul_ui (prodtmp, prodtmp, n - 1);
	mpz_set_ui (zi, i + 1);
	mpz_invert (inv, zi, P);
	mpz_mul (S[n - 1 - i], S[n - 1 - i], inv);
    }
    list_mod (S, S, n, P);
    time_calcs = cputime () - time_calcs;
#ifndef NDEBUG
    printf ("S = ");
    pari_poly (S, n - 1, 0);
    printf ("\n");
#endif

    /* We need the high part of the product tree[0] and S */
    time_calcr = cputime ();
    R = malloc (2 * n * sizeof (*R));
    for (i = 0; i < 2 * n; i++)
	mpz_init (R[i]);
    /* TODO : should use list_mul_high */
    list_mul (R, tree[0], n, 1, S, n, 0, T);
    R = R + n; /* "divide by x^n */
    /* R has degree n-1, and we want R(0), R(1), ... R(n-1) */
    list_mod (R, R, n, P);
    time_calcr = cputime () - time_calcr;
#if 0
    for (i = 0; i < n; i++)
	mpz_set_ui (R[i], 0);
    mpz_set_ui (R[0], 1);
    mpz_set_ui (R[1], 1);
#endif
#ifndef NDEBUG
    printf ("R = ");
    pari_poly (R, n - 1, 0);
    printf ("\n");
#endif
    time_calcinv = cputime ();
    invTop = malloc (n * sizeof (*invTop));
    for (i = 0; i < n; i++)
	mpz_init (invTop[i]);
    mpz_init (tree[0][n]);
    mpz_set_ui (tree[0][n], 1);
    PolyInvert (invTop, tree[0] + 1, n, T, P);
    time_calcinv = cputime () - time_calcinv;
#ifndef NDEBUG
    printf ("Inv = ");
    pari_poly (invTop, n - 1, 0);
    printf ("\n");
#endif
    time_polyeval = cputime ();
    polyeval_tellegen (R, n, tree + 1, T, tempspace, invTop, P);
    time_polyeval = cputime () - time_polyeval;
    time_adj = cputime ();
    mpz_init (adj);
    mpz_set_ui (adj, n - 1);
    for (i = n - 2; i > 1; i--) {
	mpz_mul_ui (adj, adj, i);
	mpz_mod (adj, adj, P);
    }
    if (n % 2 == 0) {
	mpz_neg (adj, adj);
    }
    for (i = 0; i < n; i++) {
	mpz_invert (inv, adj, P);
	mpz_mul (R[i], R[i], inv);
	mpz_mod (R[i], R[i], P);
	mpz_mul_ui (adj, adj, i + 1);
	if (i < n - 2) {
	    mpz_set_ui (inv, n - 1 - i);
	    mpz_invert (inv, inv, P);
	    mpz_mul (adj, adj, inv);
	}
	mpz_neg (adj, adj);
	mpz_mod (adj, adj, P);
    }
    time_adj = cputime () - time_adj;
    denum = malloc ((n+1) / 2 * sizeof (*denum));
    /* Now, reconstruct the fractions */
    time_refrac = cputime ();
    for (i = 0; i < (n+1)/2; i++) {
	mpz_init (denum[i]);
	frac_rec (b[i], denum[i], R[i], P, B);
    }
    time_refrac = cputime () - time_refrac;
	
#if 0
    printf ("Computing lcm %d\n", cputime () - usage);
    mpz_set (d, denum[0]);
    for (i = 1; i < (n+1) / 2; i++) {
	printf ("i = %d\n", i);
	mpz_lcm (d, d, denum[i]);
    }
#else
    time_lcm = cputime ();
    /* compute (multiple of) common denominator */
    mpz_fac_ui (fac_n, n);
    mpz_set_ui (num, 1);
    for (k = 2; k <= n; k++)
	mpz_lcm_ui (num, num, k);
    mpz_divexact (num, fac_n, num);	/* n! / lcm(1, 2, ..., n) */
    mpz_mul_ui (num, num, n);
    mpz_set_ui (pow_n, n);
    mpz_set_ui (den, 1);
    for (k = 1; k <= n; k++) {
	mpz_mul_ui (pow_n, pow_n, n);	/* pow_n = n^(k+1) */
	g = mpz_gcd_ui (NULL, pow_n, k + 1);
	mpz_divexact_ui (b[0], pow_n, g);
	mpz_gcd (num, num, b[0]);
	mpz_lcm_ui (den, den, (k + 1) / g);
    }
    mpz_gcd (num, num, den);
    mpz_divexact (den, den, num);
    mpz_mul (d, den, fac_n);
    time_lcm = cputime () - time_lcm;
#endif
    for (i = 0; i < (n+1) / 2; i++) {
	printf ("i = %d\n", i);
	mpz_div (denum[i], d, denum[i]);
	mpz_mul (b[i], b[i], denum[i]);
    }
    mpz_clear (P);
    mpz_clear (zi);
    mpz_clear (inv);
    printf ("%d [Time adjusting]\n", time_adj);
    printf ("%d [Time bound]\n", time_bound);
    printf ("%d [Time inv]\n", time_calcinv);
    printf ("%d [Time computing P]\n", time_calcp);
    printf ("%d [Time computing R]\n", time_calcr);
    printf ("%d [Time computing S]\n", time_calcs);
    printf ("%d [Time computing lcm]\n", time_lcm);
    printf ("%d [Time polyeval]\n", time_polyeval);
    printf ("%d [Time producttree]\n", time_ptree);
    printf ("%d [Time frac_rec]\n", time_refrac);
}

/* Compute the polynomial product (x-n0)(x-(n0+1))...(x-n1)
 * in tree[0][0], tree[0][1]
 * and all intermediate products in tree[1][0 -> (n1-n0)/2]...
 * Reduced mod n.
 */

void ProductTree (listz_t *tree, unsigned int n0, unsigned int n1, listz_t T, mpz_t n)
{
    unsigned int l, r, i;

#if VERBOSE >= 3
    printf ("Appel a ProductTree n0 = %d n1 = %d\n", n0, n1);
#endif
    if (n0 == n1) {
	mpz_init (tree[0][n0]);
	mpz_set_ui (tree[0][n0], n0);
	mpz_neg (tree[0][n0], tree[0][n0]);
	mpz_mod (tree[0][n0], tree[0][n0], n);
	return;
    }

    l = (n1 - n0 + 1) / 2;
    r = n1 - n0 + 1 - l;

    /* Compute left and right subtrees */
    ProductTree (tree + 1, n0, n0 + l - 1, T, n);
#if 0
    printf ("Sous-arbre de gauche: ");
    for (i = n0; i < n0 + l; i++) {
	mpz_out_str (stdout, 10, tree[1][i]);
	printf (" ");
    }
    printf ("\n");
#endif
    ProductTree (tree + 1, n0 + l, n1, T, n);
#if 0
    printf ("Sous-arbre de droite: ");
    for (i = n0 + l; i < n1; i++) {
	mpz_out_str (stdout, 10, tree[1][i]);
	printf (" ");
    }
    printf ("\n");
#endif
    /* This level's product */
#if VERBOSE >= 3
    printf ("(");
    pari_poly (tree[1] + n0, l);
    printf (") * (");
    pari_poly (tree[1] + n0 + l, r);
    printf (") - (");
#endif
    for (i = n0; i <= n1; i++)
	mpz_init (tree[0][i]);
    list_mul (tree[0] + n0, tree[1] + n0 + l, r, 1, tree[1] + n0, l, 1, T);
    list_mod (tree[0] + n0, tree[0] + n0, n1 - n0 + 1, n);
#if VERBOSE >= 3
    pari_poly (tree[0] + n0, l + r);
    printf (")\n");
#endif
#if 0
    printf ("Notre arbre : ");
    for (i = n0; i < n1; i++) {
	mpz_out_str (stdout, 10, tree[0][i]);
	printf (" ");
    }
    printf ("\n");
#endif
}

void pari_poly (mpz_t *p, int n, int monic)
{
    int i;

    for (i = 0; i < n; i++) {
	mpz_out_str (stdout, 10, p[i]);
	printf (" * x^%d", i);
	printf (" + ");
    }
    if (monic)
	printf ("x^%d", n);
    else {
	mpz_out_str (stdout, 10, p[n]);
	printf (" * x^%d", n);
    }
}

/* Error bound for composed NC generic rule for n points on [a, b].
 * This accounts only for the method error.
 * M is a bound of |f^(n)| for even n and of |f^(n+1)| for odd n
 * on [a, b].
 * m is the order of the composition.
 */

void nc_closed_error (mpfr_t res, mpfr_t a, mpfr_t b, unsigned int m, unsigned int n, mpfr_t M)
{
    mpfr_t h;
    unsigned int power;

    mpfr_init (h);

    mpfr_sub (h, b, a, GMP_RNDU);
    mpfr_div_ui (h, h, n - 1, GMP_RNDU);
    mpfr_div_ui (h, h, m, GMP_RNDU);

    if (n & 1)
	power = n + 2;
    else
	power = n + 1;

    mpfr_pow_ui (res, h, power, GMP_RNDU);
    mpfr_mul (res, res, M, GMP_RNDU);
    /* with odd n, |c_n| <= 1/8
     * with even n, |c_n| <= 1/4
     */

    if (n & 1)
	mpfr_mul_2si (res, res, -3, GMP_RNDU);
    else
	mpfr_mul_2si (res, res, -2, GMP_RNDU);
    mpfr_mul_ui (res, res, m, GMP_RNDU);
    mpfr_clear (h);
}
    
/* Closed NC generic rule for n points on [a, b] 
   wp is the working precision
 */
void
gen_nc_closed (mpfr_t errb, mpfr_t res, mpfr_t a, mpfr_t b, unsigned int n,
	       mpz_t *coeffs, mpz_t d, aqfunc_t f, mpfr_t m, mp_prec_t wp,
	       mpfr_t **fxs)
{
    mpfr_t D, U;
    mpfr_t fx;
    mpfr_t x, foo, bar, apb;
    mpfr_t acc;
    mpfr_t t1, t2, err_stat, err_diff, err_eval;
    mpz_t gee;
    unsigned int i;
    mpfr_t curr_deltati, max_deltati;
#ifdef DTIS
    mpfr_t dtis;
#endif

    assert (n >= 2);

    mpfr_init2 (acc, wp);
    mpfr_set_ui (acc, 0, GMP_RNDN);

    mpfr_init2 (x, wp);
    mpfr_init2 (foo, wp);
    mpfr_init2 (bar, wp);
    mpfr_init2 (curr_deltati, wp);
    mpfr_init2 (max_deltati, wp);
#ifdef DTIS
    mpfr_init2 (dtis, wp);
#endif
    mpfr_set (x, a, GMP_RNDN);

    mpfr_init2 (fx, wp);
    mpfr_init2 (apb, wp);
    mpfr_add (apb, a, b, GMP_RNDN);

    for (i = 0; i < n/2; i++) {
	f (*fxs[i], x);
	/* Compute x' such that x + x' = a + b */
	mpfr_sub (x, apb, x, GMP_RNDN);
	f (foo, x);
	mpfr_add (*fxs[i], *fxs[i], foo, GMP_RNDN);
	mpfr_mul_z (*fxs[i], *fxs[i], coeffs[i], GMP_RNDN);
	/* We need to keep track of the max error bound on
	 * ti = f(xi) * ni. It's not merely given by ulp(ti) */
	mpfr_setulp (curr_deltati, *fxs[i]);
#if VERBOSE >= 3
	printf ("\nLigne %d curr_deltati = ", __LINE__);
	mpfr_out_str (stdout, DBASE, 5, curr_deltati, GMP_RNDN);
#endif
	mpfr_mul_ui (curr_deltati, curr_deltati, 3, GMP_RNDN);
#if VERBOSE >= 3
	printf ("\nLigne %d curr_deltati = ", __LINE__);
	mpfr_out_str (stdout, DBASE, 5, curr_deltati, GMP_RNDN);
#endif
	mpfr_mul_2si (curr_deltati, curr_deltati, -1, GMP_RNDN);
#if VERBOSE >= 3
	printf ("\nLigne %d curr_deltati = ", __LINE__);
	mpfr_out_str (stdout, DBASE, 5, curr_deltati, GMP_RNDN);
#endif
	/* 3/2 ulp(ti) */

	/* Duh. ulp(0) fails of course */
	if (mpfr_zero_p (x))
	    mpfr_set_ui (foo, 0, GMP_RNDN);
	else
	    mpfr_setulp (foo, x);
	mpfr_mul_ui (foo, foo, 6, GMP_RNDU);
	if (mpz_sgn (coeffs[i]) > 0) /* abs */
	    mpfr_mul_z (foo, foo, coeffs[i], GMP_RNDU);
	else {
	    mpfr_mul_z (foo, foo, coeffs[i], GMP_RNDD);
	    mpfr_abs (foo, foo, GMP_RNDU);
	}
	
	mpfr_mul (foo, foo, m, GMP_RNDU);
	/* foo = 6 ni m ulp(xi) */
	
	mpfr_add (curr_deltati, curr_deltati, foo, GMP_RNDU);
#if VERBOSE >= 3
	printf ("\nLigne %d curr_deltati = ", __LINE__);
	mpfr_out_str (stdout, DBASE, 5, curr_deltati, GMP_RNDN);
	printf ("\nm = ");
	mpfr_out_str (stdout, DBASE, 5, m, GMP_RNDN);
#endif

#ifdef DTIS
	if (i == 0)
	    mpfr_set (dtis, curr_deltati, GMP_RNDN);
	else
	    mpfr_add (dtis, dtis, curr_deltati, GMP_RNDU);
#else
	if (i == 0)
	    mpfr_set (max_deltati, curr_deltati, GMP_RNDN);
	if (mpfr_cmp (curr_deltati, max_deltati) > 0) {
	    mpfr_set (max_deltati, curr_deltati, GMP_RNDN);
	}
#endif
	/* Compute new x */
	mpfr_mul_ui (foo, a, n - i - 2, GMP_RNDN);
	mpfr_mul_ui (bar, b, i + 1, GMP_RNDN);
	mpfr_add (x, foo, bar, GMP_RNDN);
	mpfr_div_ui (x, x, n - 1, GMP_RNDN);
    }
    if (n & 1) {
	mpfr_mul_2si (x, apb, -1, GMP_RNDN);
	f (*fxs[n/2], x);
	mpfr_mul_z (*fxs[n/2], *fxs[n/2], coeffs[n/2], GMP_RNDN);
	/* TODO : deltati... */
    }

    /* TODO : directly calling mpfr_sum is not suitable here */
    mpfr_sum (acc, (mpfr_ptr *) fxs, (n+1)/2, GMP_RNDN);

    mpfr_init2 (D, wp);
    mpfr_init2 (U, wp);

    /* TODO: Look for lost significant bits here. */
    mpfr_sub (D, b, a, GMP_RNDN);

    mpz_init (gee);
    mpz_mul_ui (gee, d, n - 1);
    mpfr_div_z (U, acc, gee, GMP_RNDN);
    mpz_clear (gee);
    mpfr_mul (res, U, D, GMP_RNDN);

    /* Error bound */
    mpfr_init (t1);
    mpfr_init (t2);
    mpfr_init (err_eval);
    mpfr_init (err_stat);
    mpfr_init (err_diff);

    /* The formerly known as "static" error */
    mpfr_set_ui (t1, 19, GMP_RNDU);
    mpfr_set_ui (t2, 35, GMP_RNDU);
    mpfr_mul_2si (t2, t2, -wp - 1, GMP_RNDU);
    mpfr_add (err_stat, t1, t2, GMP_RNDU);
    mpfr_mul_2si (err_stat, err_stat, mpfr_get_exp (res) - mpfr_get_prec (res) + 1, GMP_RNDU);

#ifdef DTIS
    mpfr_set_ui (err_eval, 1, GMP_RNDU);
#else
    mpfr_set_ui (err_eval, n, GMP_RNDU);
#endif
    mpfr_div_ui (err_eval, err_eval, n - 1, GMP_RNDU);
    mpfr_mul_ui (err_eval, err_eval, 5, GMP_RNDU);
    mpfr_div_z (err_eval, err_eval, d, GMP_RNDU);
    mpfr_mul (err_eval, err_eval, D, GMP_RNDU);
    mpfr_mul_2si (err_eval, err_eval, -1, GMP_RNDU);
    mpfr_set_ui (t1, 1, GMP_RNDU);
    mpfr_mul_2si (t1, t1, -wp, GMP_RNDU);
    mpfr_add_ui (t1, t1, 1, GMP_RNDU);
    mpfr_mul (err_eval, err_eval, t1, GMP_RNDU);
#ifdef DTIS
    mpfr_mul (err_eval, dtis, err_eval, GMP_RNDU);
#else
    mpfr_mul (err_eval, max_deltati, err_eval, GMP_RNDU);
#endif

    mpfr_setulp (t1, b);
    mpfr_setulp (t2, a);
    mpfr_add (err_diff, t1, t2, GMP_RNDU);
    mpfr_abs (U, U, GMP_RNDU);
    mpfr_mul (err_diff, err_diff, U, GMP_RNDU);
    mpfr_mul_2si (err_diff, err_diff, -1, GMP_RNDU);

    mpfr_add (errb, err_diff, err_eval, GMP_RNDU);
    mpfr_add (errb, errb, err_stat, GMP_RNDU);

    mpfr_clear (acc);
    mpfr_clear (fx);
    mpfr_clear (D);
    mpfr_clear (U);
    mpfr_clear (x);
    mpfr_clear (foo);
    mpfr_clear (bar);
    mpfr_clear (apb);
    mpfr_clear (t1);
    mpfr_clear (t2);
    mpfr_clear (err_eval);
    mpfr_clear (err_stat);
    mpfr_clear (err_diff);
    mpfr_clear (curr_deltati);
    mpfr_clear (max_deltati);
}

double ulp_erreur (mpfr_t ret, mpfr_t exact, mpfr_t approche)
{
    mpfr_t diff;
    double r;
    mpfr_init2 (diff, mpfr_get_prec (exact));

    /* BUGFIX: ajout du test somme nulle */

    if ((mpfr_cmp_ui (approche, 0) == 0) && (mpfr_cmp_ui (exact, 0) != 0)) {
	mpfr_set_inf (ret, 1);
	mpfr_clear (diff);
	return 0;
    }
    if (mpfr_sgn (exact) >= 0) {
	if (mpfr_cmp (exact, approche) >= 0)
	    mpfr_sub (diff, exact, approche, GMP_RNDU);
	else
	    mpfr_sub (diff, approche, exact, GMP_RNDU);
    }
    else {
	if (mpfr_cmp (exact, approche) >= 0)
	    mpfr_sub (diff, approche, exact, GMP_RNDD);
	else
	    mpfr_sub (diff, exact, approche, GMP_RNDD);
    }
    mpfr_abs (diff, diff, GMP_RNDU);
    mpfr_mul_2si (diff, diff, mpfr_get_prec (approche) - approche->_mpfr_exp,
		  GMP_RNDN);

    r = mpfr_get_d1 (diff);
    mpfr_swap (diff, ret);
    mpfr_clear (diff);
    return r;
}

/* Compose m times the closed n points NC method on [a, b]
 * M is a bound (absolute value) on the n-th (or n+1-th, depends on the
 * parity of n) derivative of f
 * m1 a bound (absolute value) on the derivative of f.
 * */
void
compose_nc (mpfr_t res, mpfr_t gerrb, mpfr_t a, mpfr_t b, unsigned int m,
	    unsigned int n, aqfunc_t f, mpfr_t m1, mp_prec_t wp)
{
    mpfr_t bmin, bmax;
    mpfr_t acc;
    mpfr_t foo, bar, errb;
    unsigned int step, i;
    mpz_t *coeffs;
    mpz_t d;
    int st;
    mpfr_t nc_error;
    mpfr_t **ncs;
    mpfr_t **fxs;
    assert (m > 0);

    /* Initialize and compute the Bnk */
    coeffs = malloc (((n + 1) / 2) * sizeof (*coeffs));
    for (i = 0; i < (n + 1) / 2; i++)
	mpz_init (coeffs[i]);
    mpz_init (d);

    mpfr_init (errb);

    st = cputime ();
    compute_bnk (coeffs, d, n - 1);
    printf ("compute_bnk took %dms\n", cputime () - st);

    mpfr_init2 (bmin, wp);
    mpfr_init2 (bmax, wp);
    mpfr_init2 (foo, wp);
    mpfr_init2 (bar, wp);
    mpfr_init2 (acc, wp);
    mpfr_init2 (nc_error, wp);

    /* Initialise tab for nc_gen calls */
    fxs = malloc (((n+1)/2) * sizeof (*fxs));
    for (i = 0; i < (n+1)/2; i++) {
        fxs[i] = malloc (sizeof (mpfr_t));
        mpfr_init2 (*fxs[i], wp);
    }
    /* Initialise tab for nc_gen return values */
    ncs = malloc (m * sizeof (*ncs));
    for (i = 0; i < m; i++) {
        ncs[i] = malloc (sizeof (mpfr_t));
        mpfr_init2 (*ncs[i], wp);
    }

    mpfr_set (bmin, a, GMP_RNDN);
    mpfr_set_ui (acc, 0, GMP_RNDN);
    mpfr_set_ui (gerrb, 0, GMP_RNDN);

    for (step = 0; step < m; step++) {
	/* Compute new bmax */
	mpfr_mul_ui (foo, a, m - step - 1, GMP_RNDN);
	mpfr_mul_ui (bar, b, step + 1, GMP_RNDN);
	mpfr_add (bmax, foo, bar, GMP_RNDN);
	mpfr_div_ui (bmax, bmax, m, GMP_RNDN);
	/* Compute NC-closed n point on [bmin, bmax] */
	gen_nc_closed (errb, *ncs[step], bmin, bmax, n, coeffs, d, f, m1, wp,
		       fxs);
	mpfr_add (gerrb, gerrb, errb, GMP_RNDU);
	/* Set bmin for next step */
	mpfr_set (bmin, bmax, GMP_RNDN);
#ifndef NDEBUG
	if (mpfr_nan_p (acc)) {
	    printf ("Nan. step = %d.\n", step);
	    abort ();
	}
#endif
    }

    mpfr_sum (res, (mpfr_ptr *) ncs, m, GMP_RNDN);
    /* mpfr_set (res, acc, GMP_RNDN); */
    for (i = 0; i < (n + 1) / 2; i++)
	mpz_clear (coeffs[i]);
    mpz_clear (d);
    free (coeffs);
    mpfr_clear (acc);
    mpfr_clear (bar);
    mpfr_clear (bmax);
    mpfr_clear (bmin);
    mpfr_clear (foo);
    for (i = 0; i < (n+1)/2; i++) {
	mpfr_clear (*fxs[i]);
	free (fxs[i]);
    }
    free (fxs);
    for (i = 0; i < m; i++) {
	mpfr_clear (*ncs[i]);
	free (ncs[i]);
    }
    free (ncs);
}

#ifdef TEST2
void sin_wrapper (mpfr_t res, mpfr_t x, mp_prec_t rndm)
{
    mpfr_sin (res, x, rndm);
}

void square_wrapper (mpfr_t res, mpfr_t x, mp_prec_t rndm)
{
    mpfr_mul (res, x, x, rndm);
}

void exp_wrapper (mpfr_t res, mpfr_t x, mp_prec_t rndm)
{
    mpfr_exp (res, x, rndm);
}

void nth_deriv_bound_exp (mpfr_t res, mpfr_t a, mpfr_t b, unsigned int n)
{
    mpfr_exp (res, b, GMP_RNDU);
}

/* t*log(1+t) */
void bailey1 (mpfr_t res, mpfr_t t, mp_prec_t rndm)
{
    mpfr_log1p (res, t, rndm);
    mpfr_mul (res, res, t, rndm);
}

int main (int argc, char *argv[])
{
    unsigned int n, m;
    mp_prec_t prec, wp = 0;

    mpfr_t a, b;
    mpfr_t res;
    mpfr_t M, m1;
    mpfr_t exact_value;
    mpfr_t erreur;
    aqderiv_bound_t *nth_deriv_bound;


    m = 1;
    n = 42;

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
	else if (argc >= 3 && strcmp (argv[1], "-n") == 0) {
	    n = atoi (argv[2]);	/* Order of NC method */
	    argv += 2;
	    argc -= 2;
	}
	else if (argc >= 3 && strcmp (argv[1], "-m") == 0) {
	    m = atoi (argv[2]);	/* number of intervals */
	    argv += 2;
	    argc -= 2;
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
    mpfr_init2 (erreur, wp);

    mpfr_set_ui (a, LEFT_BOUND, GMP_RNDN);
    mpfr_set_ui (b, RIGHT_BOUND, GMP_RNDN);

    printf ("%d [n]\n%d [m]\n%lu [wp]\n", n, m, wp);

    mpfr_init (M);
    mpfr_init (m1);
#if (FUNCTION == exp_wrapper)
    nth_deriv_bound = nth_deriv_bound_exp;
#endif
    nth_deriv_bound (m1, a, b, 1);
    if (n & 1)
	nth_deriv_bound (M, a, b, n + 1);
    else
	nth_deriv_bound (M, a, b, n);
    compose_nc (res, a, b, m, n, FUNCTION, M, m1, wp);

    printf ("val = ");
    mpfr_out_str (NULL, DBASE, 0, res, GMP_RNDN); printf ("\n");
#if (FUNCTION == exp_wrapper)
    mpfr_set_ui (exact_value, RIGHT_BOUND, GMP_RNDN);
    mpfr_exp (exact_value, exact_value, GMP_RNDN);
    mpfr_set_ui (a, LEFT_BOUND, GMP_RNDN);
    mpfr_exp (a, a, GMP_RNDN);
    mpfr_sub (exact_value, exact_value, a, GMP_RNDN);
    printf ("Exact val: ");
    mpfr_out_str (NULL, DBASE, 0, exact_value, GMP_RNDN); printf ("\n");
    /* derreur = ulp_erreur (erreur, exact_value, res);
    printf ("\nErreur (inexact): %f ulp\n", derreur); */
    mpfr_sub (exact_value, exact_value, res, GMP_RNDN);
    mpfr_abs (exact_value, exact_value, GMP_RNDN);
    mpfr_out_str (NULL, DBASE, 0, exact_value, GMP_RNDN);
    printf (" [measured error]\n");
#endif
    mpfr_clear (a);
    mpfr_clear (b);
    mpfr_clear (erreur);
    mpfr_clear (exact_value);
    mpfr_clear (res);
    mpfr_clear (M);
    mpfr_clear (m1);
    return 0;
}

#endif

#ifdef TEST3

int main (void) {
    mpfr_t x, ulpx;

    mpfr_init2 (x, 4);
    mpfr_init (ulpx);

    mpfr_set_ui (x, 3, GMP_RNDN);
    mpfr_setulp (ulpx, x);
    mpfr_out_str (stdout, DBASE, 0, x, GMP_RNDN);
    printf (" [x]\n");
    mpfr_out_str (stdout, DBASE, 0, ulpx, GMP_RNDN);
    printf (" [ulpx]\n");

}

#endif

#ifdef TEST1

int main (int argc, char *argv[])
{
    unsigned long n = atoi (argv[1]);
    unsigned long k;
    mpz_t *b, d;

    mpz_init (d);
    b = (mpz_t *) malloc ((n / 2 + 1) * sizeof (mpz_t));
    for (k = 0; k <= n / 2; k++)
	mpz_init (b[k]);
    compute_bnk (b, d, n);
    for (k = 0; k <= n; k++) {
	printf ("b[%lu]=", k);
	mpz_out_str (stdout, 10, b[(k <= n / 2) ? k : n - k]);
	printf ("/");
	mpz_out_str (stdout, 10, d);
	printf ("\n");
    }
    for (k = 0; k <= n / 2; k++)
	mpz_clear (b[k]);
    free (b);
    mpz_clear (d);

    return 0;
}

#endif

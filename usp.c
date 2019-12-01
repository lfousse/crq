/* usp.c -- Uspensky's root isolation.

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
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA.

Copyright Guillaume Hanrot 1999, 2002,
with some improvements by Fabrice Rouillier.
based on Descartes' rule, following the procedure described by
Rouillier & Zimmermann.

This file was taken from the Quadrics Intersection (QI) software :
http://www.loria.fr/equipes/vegas/qi/
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <math.h>
#include "usp.h"

#define vali(x) mpz_scan1((x), 0)
#define ilog2(a) mpz_sizeinbase(a,2)
#define TOT_POS -1

/* Probably out of memory long since on nontrivial examples, but anyway. */
#define DEPTH_MAX 1024

static unsigned long usign;
static int b1 = 0, b2 = 0;

interval *clean_output(interval * roots, unsigned long deg,
		       unsigned int nbroot);

int RemoveContent(mpz_t * P, unsigned long deg)
{
    unsigned long cont, i, z;

    i = 0;
    while (mpz_sgn(P[i]) == 0)
	i++;
    cont = vali(P[i]);

    for (; (i <= deg) && cont; i++) {
	if (mpz_sgn(P[i]) != 0) {
	    z = vali(P[i]);
	    if (z < cont)
		cont = z;
	}
    }
    if (cont == 0)
	return 0;

    for (i = 0; i <= deg; i++)
	mpz_fdiv_q_2exp(P[i], P[i], cont);

    return cont;
}

/* Johnson's bound : 2*max(abs(-ai/an)^(1/(n-i)),
   the maximum being taken over the i with sgn(a_i) != sgn(an) */
long bound_roots(mpz_t * t, unsigned long deg)
{
    unsigned long i;
    long maxpow, currpow, currpow2, lan, tpos = 1;

    currpow = currpow2 = 0;
    lan = ilog2(t[deg]) - 1;	/* puiss de 2 < an */

    maxpow = -lan;

    for (i = 0; i < deg; i++) {
	if (mpz_sgn(t[deg]) != mpz_sgn(t[i])) {
	    tpos = 0;
	    currpow = ilog2(t[i]);
	    currpow -= lan;	/* 2^currpow >= abs(-ai/an) */

	    if (currpow > 0)
		currpow2 = currpow / (deg - i);
	    else
		currpow2 = ((-currpow) / (deg - i));

	    /* Bug fix Sylvain in parentheses */
	    if (currpow2 * ((long) (deg - i)) != currpow)
		currpow2++;
	    /* 2^currpow2 >= abs(-ai/an)^(1/(n-i)) */

	    if (currpow2 > maxpow)
		maxpow = currpow2;
	}
    }
    if (tpos == 1)
	return -1;

    /* here 2^maxpow > max(abs(-ai/an)^(1/(n-i)), add one to get the bound */
    maxpow++;
    return maxpow;
}

/* From the polynomial P, of degree deg, and the bound b such that
        all positive roots of P are <= 2^b, compute a polynomial Q0
        which has all its real roots in ]0, 1[, namely P(X*2^b). */
int change_pol(mpz_t * P, unsigned long b, unsigned long deg)
{
    unsigned long i, j = b;

    for (i = 1; i <= deg; i++, j += b) {
	mpz_mul_2exp(P[i], P[i], j);
    }
    return RemoveContent(P, deg);
}

/* Computes P(X*2^k) */
int Homoth(mpz_t * P, long k, unsigned long deg)
{
    unsigned long i;
    long j;

    if (k > 0) {
	j = k;
	for (i = 1; i <= deg; i++, j += k)
	    mpz_mul_2exp(P[i], P[i], j);
    } else {
	j = deg * (-k);
	for (i = 0; i < deg; i++, j += k)
	    mpz_mul_2exp(P[i], P[i], j);
    }

    /* Remove possible large power of 2 in content */
    return RemoveContent(P, deg);
}

/* Replaces P by the polynomial P(X+1) */
void X2XP1(mpz_t * P, unsigned long deg)
{
    long i, j;
    for (i = 0; i <= (long) deg - 1; i++)
	for (j = deg - 1; j >= i; j--)
	    mpz_add(P[j], P[j], P[j + 1]);
}

/* Replaces P by the polynomial P(X-1) */
void X2XM1(mpz_t * P, unsigned long deg)
{
    long i, j;
    for (i = 0; i <= (long) deg - 1; i++)
	for (j = deg - 1; j >= i; j--)
	    mpz_sub(P[j], P[j], P[j + 1]);
}

/* Replaces P by the polynomial P(X+C) */
void X2XPC(mpz_t * P, mpz_t c, unsigned long deg)
{
    long i, j;
    mpz_t foo;
    mpz_init(foo);
    for (i = 0; i <= (long) deg - 1; i++) {
	for (j = deg - 1; j >= i; j--) {
	    mpz_mul(foo, P[j + 1], c);
	    mpz_add(P[j], P[j], foo);
	}
    }
    mpz_clear(foo);

    return;
}

/* Number of sign changes in the coefficients of P(1/(X+1)) 
   (Descartes' rule) */
unsigned long Descartes(mpz_t * P, unsigned long deg, long sigh,
			long *flag)
{
    unsigned long nb = 0;
    long s, t;
    unsigned long i, j;
    mpz_t *Q;

    if (deg == 0)
	return 0;

    /* 
       Prune the computation if all the coefficients are of the sign of P[deg] 
       In that case any subsequent interval shall have the same property, 
       we put *flag at 1 to point this to Uspensky_rec.
     */
    j = deg;
    t = mpz_sgn(P[j]);
    while (j >= 0 && mpz_sgn(P[j]) == t) {
	if (j == 0)
	    break;
	j--;
    }

    if (j < 0) {
	*flag = -1;
	return nb;
    }

    Q = (mpz_t *) malloc((deg + 1) * sizeof(mpz_t));

    /*check_alloc((!Q),"mpz_t *","Descartes"); */

    for (i = 0; i <= deg; i++)
	mpz_init_set(Q[i], P[i]);

    for (j = 0; j <= deg - 1; j++)
	mpz_add(Q[j + 1], Q[j + 1], Q[j]);

    s = mpz_sgn(Q[deg]);

    *flag = s && (s == mpz_sgn(P[0])) && (s == -sigh);

    for (i = 1; i <= deg - 1; i++) {
	/* 
	   Prune the computation if all further coefficients are of the sign of 
	   Q[deg-i] 
	 */
	j = deg - i;
	t = s;
	while (j >= 0 && t == 0) {
	    t = mpz_sgn(Q[j]);
	    j--;
	}
	while (j >= 0 && mpz_sgn(Q[j]) == t) {
	    if (j == 0)
		break;

	    j--;
	}
	if (j < 0) {
	    for (i = 0; i <= deg; i++)
		mpz_clear(Q[i]);
	    free(Q);
	    return nb;
	}

	for (j = 0; j <= deg - i - 1; j++)
	    mpz_add(Q[j + 1], Q[j + 1], Q[j]);

	if (s == 0) {
	    s = mpz_sgn(Q[deg - i]);
	} else if (s == -mpz_sgn(Q[deg - i])) {
	    if ((nb == 1 && !*flag) || nb == 2) {
		for (i = 0; i <= deg; i++)
		    mpz_clear(Q[i]);
		free(Q);
		return (nb + 1);
	    }
	    nb++;
	    s = -s;
	}
    }

    if (s == -mpz_sgn(Q[0]))
	nb++;
    for (i = 0; i <= deg; i++)
	mpz_clear(Q[i]);
    free(Q);

    return nb;
}

/* Returns the sign of P(1/2) */
long evalhalf(mpz_t * P, unsigned long deg)
{
    long j;
    int ret;
    mpz_t x, y;

    mpz_init_set(x, P[deg]);
    mpz_init(y);

    for (j = deg - 1; j >= 0; j--) {
	mpz_mul_2exp(y, P[j], deg - j);
	mpz_add(x, x, y);
    }

    mpz_clear(y);
    ret = mpz_sgn(x);
    mpz_clear(x);
    return ret;
}

void add_root(interval * roots, mpz_t c, int k, unsigned int flag,
	      unsigned int nbroot)
{
    int b = (usign ? b1 : b2);

    mpz_init(roots[nbroot].c);

    if (k <= b) {
	if (usign) {
	    mpz_neg(roots[nbroot].c, c);
	    /* Bug fix Sylvain when root is exact */
	    if (!flag)
		mpz_sub_ui(roots[nbroot].c, roots[nbroot].c, 1);
	    mpz_mul_2exp(roots[nbroot].c, roots[nbroot].c, b - k);
	} else
	    mpz_mul_2exp(roots[nbroot].c, c, b - k);

	roots[nbroot].k = k - b;
	roots[nbroot].isexact = flag;

	return;
    } else {
	if (usign) {
	    mpz_neg(roots[nbroot].c, c);
	    /* Bug fix Sylvain when root is exact */
	    if (!flag)
		mpz_sub_ui(roots[nbroot].c, roots[nbroot].c, 1);
	} else
	    mpz_set(roots[nbroot].c, c);

	roots[nbroot].k = k - b;
	roots[nbroot].isexact = flag;
    }

    return;
}

/* 
   Check interval [c/2^k, (c+1)/2^k]. The value of k is returned, this
   is necessary to know from where we come [i.e., what exactly is the 
   current polynomial P] when several recursive calls return in a row. 
   In practice, this is used to update the polynomial at HERE */
long Uspensky_rec(mpz_t * P, mpz_t c, unsigned long k, unsigned long *Deg,
		  interval * roots, unsigned int *nbroot)
{
    unsigned long i, j, nb;
    long oldk;
    long shalf, flag;
    mpz_t tmp;

    if (k > DEPTH_MAX) {
	fprintf(stderr, "Maximal depth reached.\n");
	exit(-1);
    }

    mpz_init(tmp);

    /* Check whether c/2^k is a root */
    if (mpz_cmp_ui(P[0], 0) == 0) {
	i = 1;
	while (mpz_cmp_ui(P[i], 0) == 0) {
	    i++;
	}

	for (j = 0; j < i; j++) {
	    add_root(roots, c, k, 1, *nbroot);
	    (*nbroot)++;
	}

	*Deg -= i;		/* Update the polynomial */
	for (j = 0; j <= *Deg; j++, i++)
	    mpz_set(P[j], P[i]);
    }

    /* 
       Compute the sign of P(1/2) ; thus if Descartes bound is 2, 
       whereas sign(P(0)) = sign(P(1)) = -sign(P(1/2)) we have
       found two roots. 
     */
    shalf = evalhalf(P, *Deg);

    /* Compute Descartes' bound */
    nb = Descartes(P, *Deg, shalf, &flag);
    if (flag == TOT_POS) {
	mpz_clear(tmp);
	return TOT_POS;
    }

    switch (nb) {
    case 0:			/* no root */
	mpz_clear(tmp);
	return k;

    case 1:			/* exactly one root */
	add_root(roots, c, k, 0, *nbroot);
	(*nbroot)++;
	mpz_clear(tmp);
	return k;

    case 2:			/* if flag!=0, one root in each half of the current interval */
	if (flag) {
	    mpz_set(tmp, c);

	    mpz_mul_2exp(tmp, tmp, 1);
	    add_root(roots, tmp, k + 1, 0, *nbroot);
	    (*nbroot)++;

	    mpz_add_ui(tmp, tmp, 1);
	    add_root(roots, tmp, k + 1, 0, *nbroot);
	    (*nbroot)++;

	    mpz_clear(tmp);
	    return k;
	}

    default:	/* recursive call on each half of the interval */
	mpz_set(tmp, c);

	mpz_mul_2exp(tmp, tmp, 1);
	Homoth(P, -1, *Deg);
	oldk = Uspensky_rec(P, tmp, k + 1, Deg, roots, nbroot);
	if (oldk == TOT_POS) {
	    mpz_clear(tmp);
	    return TOT_POS;
	}

	mpz_add_ui(tmp, tmp, 1);
	X2XP1(P, *Deg);

	if (oldk > (long) k + 1)
	    Homoth(P, oldk - (k + 1), *Deg);

	oldk = Uspensky_rec(P, tmp, k + 1, Deg, roots, nbroot);
	if (oldk == TOT_POS) {
	    mpz_clear(tmp);
	    return TOT_POS;
	}

	mpz_clear(tmp);
	return oldk;
    }
}

interval *Uspensky(mpz_t * Q, unsigned long deg, unsigned int *nbroot)
{
    unsigned long deg1 = deg;
    long i, j;
    interval *roots = (interval *) malloc(deg1 * sizeof(interval));

    /*check_alloc((!Q),"interval *","Uspensky"); */

    mpz_t *P, e;
    int nb_z;

    mpz_init_set_ui(e, 0);
    *nbroot = 0;

    /* Remove 0 roots in all cases */
    nb_z = 0;
    while (mpz_sgn(Q[nb_z]) == 0)
	nb_z++;
    for (j = 0; j < nb_z; j++) {
	add_root(roots, e, 0, 1, *nbroot);
	(*nbroot)++;
    }

    deg1 = deg - nb_z;		/* Update the polynomial */
    P = (mpz_t *) malloc((deg1 + 1) * sizeof(mpz_t));

    /*check_alloc((!Q),"mpz_t *","Uspensky"); */

    for (j = 0; j <= (long) deg1; j++)
	mpz_init_set(P[j], Q[nb_z + j]);

    /* First work on the positive roots. */
    b2 = bound_roots(P, deg1);

    if (b2 < 0)
	goto NEGATIVE;

    change_pol(P, b2, deg1);

    usign = 0;
    Uspensky_rec(P, e, 0, &deg1, roots, nbroot);

    /* Change P into P(-X) to look for negative roots */
  NEGATIVE:
    deg1 = deg - nb_z;
    for (i = deg1; i >= 0; i--) {
	if (i % 2 == 1)
	    mpz_neg(P[i], Q[nb_z + i]);
	else
	    mpz_set(P[i], Q[nb_z + i]);
    }

    b1 = bound_roots(P, deg1);

    if (b1 >= 0) {
	change_pol(P, b1, deg1);
	mpz_set_ui(e, 0);
	usign = 1;
	Uspensky_rec(P, e, 0, &deg1, roots, nbroot);
    }

    /* Free memory. */
    for (i = deg - nb_z; i >= 0; i--)
	mpz_clear(P[i]);

    free(P);

    /* Now a bit of cleaning. */
    /* HAHA */
    mpz_clear(e);
    return clean_output(roots, deg, *nbroot);
}

/* 
   Check interval [c/2^k, (c+1)/2^k]. The value of k is returned, this
   is necessary to know from where we come [i.e., what exactly is the 
   current polynomial P] when several recursive calls return in a row. 
   In practice, this is used to update the polynomial at HERE */
long Uspensky_couple_rec(mpz_t * P1, mpz_t *P2, mpz_t c, unsigned long k,
			 unsigned long *Deg1, unsigned long *Deg2,
			 interval * roots, unsigned int *nbroot)
{
    unsigned long i, j, nb1, nb2;
    long oldk;
    long shalf1, shalf2, flag;
    mpz_t tmp;

    if (k > DEPTH_MAX) {
	fprintf(stderr, "Maximal depth reached.\n");
	exit(-1);
    }

    mpz_init(tmp);

    /* Check whether c/2^k is a root of P1 */
    if (mpz_cmp_ui(P1[0], 0) == 0) {
	i = 1;
	while (mpz_cmp_ui(P1[i], 0) == 0) {
	    i++;
	}

	for (j = 0; j < i; j++) {
	    add_root(roots, c, k, 1, *nbroot);
	    (*nbroot)++;
	}

	*Deg1 -= i;		/* Update the polynomial */
	for (j = 0; j <= *Deg1; j++, i++)
	    mpz_set(P1[j], P1[i]);
    }

    /* Check whether c/2^k is a root of P2 */
    if (mpz_cmp_ui(P2[0], 0) == 0) {
	i = 1;
	while (mpz_cmp_ui(P2[i], 0) == 0) {
	    i++;
	}

	for (j = 0; j < i; j++) {
	    add_root(roots, c, k, 1, *nbroot);
	    (*nbroot)++;
	}

	*Deg2 -= i;		/* Update the polynomial */
	for (j = 0; j <= *Deg2; j++, i++)
	    mpz_set(P2[j], P2[i]);
    }

    /* 
       Compute the sign of P(1/2) ; thus if Descartes bound is 2, 
       whereas sign(P(0)) = sign(P(1)) = -sign(P(1/2)) we have
       found two roots. 
     */
    shalf1 = evalhalf(P1, *Deg1);
    shalf2 = evalhalf(P2, *Deg2);

    /* Compute Descartes' bound */
    nb1 = Descartes(P1, *Deg1, shalf1, &flag);
    nb2 = Descartes(P2, *Deg2, shalf2, &flag);

    switch (nb1 + nb2) {
    case 0:			/* no root */
	mpz_clear (tmp);
	return k;

    case 1:			/* exactly one root */
	add_root (roots, c, k, 0, *nbroot);
	(*nbroot)++;
	mpz_clear (tmp);
	return k;

    default:	/* recursive call on each half of the interval */
	mpz_set(tmp, c);

	mpz_mul_2exp(tmp, tmp, 1);
	Homoth(P1, -1, *Deg1);
	Homoth(P2, -1, *Deg2);
	oldk = Uspensky_couple_rec (P1, P2, tmp, k + 1, Deg1, Deg2, roots,
				    nbroot);
	if (oldk == TOT_POS) {
	    mpz_clear(tmp);
	    return TOT_POS;
	}

	mpz_add_ui(tmp, tmp, 1);
	X2XP1(P1, *Deg1);
	X2XP1(P2, *Deg2);

	if (oldk > (long) k + 1) {
	    Homoth(P1, oldk - (k + 1), *Deg1);
	    Homoth(P2, oldk - (k + 1), *Deg2);
	}

	oldk = Uspensky_couple_rec (P1, P2, tmp, k + 1, Deg1, Deg2, roots,
				    nbroot);
	if (oldk == TOT_POS) {
	    mpz_clear(tmp);
	    return TOT_POS;
	}

	mpz_clear(tmp);
	return oldk;
    }
}

interval *Uspensky_couple (mpz_t * Q1, mpz_t *Q2, unsigned long deg1,
			   unsigned long deg2, unsigned int *nbroot)
{
    interval *roots = (interval *) malloc((deg1+deg2) * sizeof(interval));

    mpz_t e;

    mpz_init_set_ui(e, 0);
    *nbroot = 0;

    /* First work on the positive roots. */
    b2 = bound_roots(Q1, deg1);
    /* b2 = bound_roots(P2, deg2); Not needed for Legendre */

    change_pol(Q1, b2, deg1);
    change_pol(Q2, b2, deg2);

    usign = 0;
    Uspensky_couple_rec (Q1, Q2, e, 0, &deg1, &deg2, roots, nbroot);

    mpz_clear(e);
    return clean_output(roots, deg1 + deg2, *nbroot);
}

/* Output things in the right order */
interval *clean_output(interval * roots, unsigned long deg,
		       unsigned int nbroot)
{
    int i;
    if (nbroot == 0)
	return roots;
    else if (mpz_sgn(roots[nbroot - 1].c) >= 0)	/* Nothing to be done */
	return roots;
    else {
	interval *cleaned_roots =
	    (interval *) malloc(deg * sizeof(interval));

	/*check_alloc((!cleaned_roots),"interval *","clean_output"); */

	int j = 0;

	for (i = nbroot - 1; i >= 0; i--) {
	    if (mpz_sgn(roots[i].c) >= 0)
		break;

	    mpz_init_set(cleaned_roots[j].c, roots[i].c);
	    cleaned_roots[j].k = roots[i].k;
	    cleaned_roots[j].isexact = roots[i].isexact;
	    j++;
	}

	for (i = 0; i < nbroot - j; i++) {
	    mpz_init_set(cleaned_roots[i + j].c, roots[i].c);
	    cleaned_roots[i + j].k = roots[i].k;
	    cleaned_roots[i + j].isexact = roots[i].isexact;
	}

	return cleaned_roots;
    }
}

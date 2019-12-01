/* crq-impl.h -- private header file.

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

#ifndef _INTLIB_IMPL_H
#define _INTLIB_IMPL_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "usp.h"

#include "crq_structs.h"

/* Private definitions */

#ifndef DBASE
#define DBASE 10
#endif

#define DICHO(x) ((x) & 1)
#define DICHO_TRUNC(x) ((x) & 8)
#define NEWTON(x) (((x) & 1) == 0)
#define USP_COUPLE(x) ((x) & 2)
#define NEWTON_INT(x) ((x) & 4)

/* #define USP_COUPLE */
/* #define NEWTON_EXPERIMENTAL */
/* #define IMP_DICHO */

/*#define PROGRESS*/

enum { SPA = 0, SPB = 1, SPPA = 2, SPPB = 3};

#define INTERVAL(a,b) do { if (mpfr_cmp ((a), (b)) > 0) \
			    { printf ("Bad interval line %d\n", __LINE__);\
			      mpfr_dump ((a)); mpfr_dump ((b)); \
			      abort (); } } while (0)


void deriv (mpz_t *, mpz_t *, unsigned long);
void horner_eval (mpfr_t res, mpfr_t x, mpz_t *p, unsigned long n, mpfr_t, mpfr_t);
void horner_eval_int2 (mpfr_t, mpfr_t, mpfr_t, mpfr_t,
		       mpfr_t *, mpfr_t *, unsigned long);

int newton_int (mpfr_t, mpz_t *, mpz_t *, unsigned long, interval *,
		mp_prec_t, int *);

extern int newton_good;
extern int newton_bad;
void root_dicho (mpfr_t ret, interval *I, mpz_t *p, unsigned int d,
		 mp_prec_t wp, struct crq_stats *);
int root_dicho_trunc (mpfr_t ret, interval *I, mpz_t *p, mpfr_t *pl, mpfr_t *pr,
		       unsigned int d, mp_prec_t np);
void root_newton (mpfr_t ret, interval *I, mpz_t *p, mpz_t *pp,
			unsigned int d, mp_prec_t wp, struct crq_gl_opts newton_interval,
			struct crq_stats *stats);
void root_secante (mpfr_t root, interval *I, mpz_t *p, unsigned long d,
		   mp_prec_t wp, struct crq_gl_opts options, struct crq_stats *timings);

#endif /* _INTLIB_IMPL_H */

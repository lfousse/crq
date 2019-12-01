/* crq.h -- public header file.

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

#ifndef _INTLIB_H
#define _INTLIB_H

#include <mpfr.h>
#include "crq_structs.h"

/* Integration methods -- Composed */
void compose_gl (mpfr_t, mpfr_t, mpfr_t, mpfr_t, unsigned int,
		 unsigned int, aqfunc_t, mpfr_t, mp_prec_t, struct crq_gl_opts);
void compose_gl_stats (mpfr_t, mpfr_t, mpfr_t, mpfr_t, unsigned int,
		 unsigned int, aqfunc_t, mpfr_t, mp_prec_t, struct crq_gl_opts,
		 struct crq_stats *);
void compose_nc (mpfr_t, mpfr_t, mpfr_t, mpfr_t, unsigned int, unsigned int,
		 aqfunc_t f, mpfr_t, mp_prec_t);

/* Integration methods -- Simple */
void gen_nc_closed (mpfr_t, mpfr_t, mpfr_t, mpfr_t, unsigned int,
	            mpz_t *, mpz_t, aqfunc_t, mpfr_t, mp_prec_t, mpfr_t **);

void gen_gl_closed (mpfr_t res, mpfr_t errb, mpfr_t a, mpfr_t b, unsigned int n,
		    mpfr_t *coeffs, mpfr_t *roots, aqfunc_t f, mpfr_t m,
		    mp_prec_t wp, mpfr_t *fxs[]);

/* Integration methods -- Adaptive */
void adaptive_quad (mpfr_t res, mpfr_t a, mpfr_t b, aqfunc_t f, mp_prec_t prec,
		    void (*quadf)(mpfr_t, mpfr_t, mpfr_t, aqfunc_t, void *),
		    void *params);

/* error bounds */
void gl_error (mpfr_t, mpfr_t, mpfr_t b, unsigned int, unsigned int, mpfr_t M);
void nc_closed_error (mpfr_t, mpfr_t a, mpfr_t b, unsigned int m, unsigned int n, mpfr_t M);

/* Weights and abscissas computations */
void compute_glc (mpfr_t *b, mpfr_t *roots, unsigned long n, mp_prec_t wp,
		  struct crq_gl_opts options, struct crq_stats *timings);

#endif /* _INTLIB_H */

/* interval.h -- interval functions declarations.

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

#ifndef _INTERVAL_H_
#define _INTERVAL_H_
#include <stdio.h>
#include <mpfr.h>
#include "crq-impl.h"

/* Interval arithmetic */

void myint_add (mpfr_t retl, mpfr_t retr, mpfr_t al, mpfr_t ar, mpfr_t bl,
	       mpfr_t br);
void myint_sub (mpfr_t retl, mpfr_t retr, mpfr_t al, mpfr_t ar, mpfr_t bl,
	       mpfr_t br);
void myint_mul (mpfr_t retl, mpfr_t retr, mpfr_t al, mpfr_t ar, mpfr_t bl,
	       mpfr_t br, mpfr_t *tmp);
void myint_div (mpfr_t retl, mpfr_t retr, mpfr_t al, mpfr_t ar, mpfr_t bl,
	       mpfr_t br);
void myint_add_z (mpfr_t resl, mpfr_t resr, mpfr_t xl, mpfr_t xr, mpz_t b);
void horner_eval_int (mpfr_t resl, mpfr_t resr, mpfr_t xl, mpfr_t xr, mpz_t *p,
		      unsigned long n);
void horner_eval_int2 (mpfr_t resl, mpfr_t resr, mpfr_t xl, mpfr_t xr,
		       mpfr_t *pl, mpfr_t *pr, unsigned long n);
void intersect (mpfr_t resl, mpfr_t resr, mpfr_t al, mpfr_t ar, mpfr_t bl,
		mpfr_t br);
int is_in_interval (mpfr_t root, interval *I);

#endif /* _INTERVAL_H_ */

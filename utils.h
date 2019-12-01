/* utils.h -- utility functions declarations.

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
#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdio.h> /* for FILE */
#include "usp.h" /* for interval */
#include <gmp.h>

int cputime (void);

void sign_set (char *signmask, int pos, int value);
int sign_get (const char signmask, int pos);
int sign_eval (mpz_t c, mpz_t *P, unsigned long k, unsigned long n);
void abort_if_nan (mpfr_t, char *, int);
void debug_mpfr_value (mpfr_t foo, char *name);
void mpfr_setulp (mpfr_t, mpfr_t);
void dump_poly (FILE *out, mpz_t *P, unsigned long d);
void load_poly (FILE *in, mpz_t *P, unsigned long d);
void dump_interval (FILE *out, interval I);
void load_interval (FILE *in, interval *I);
#endif /* _UTILS_H_ */

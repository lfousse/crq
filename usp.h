/* usp.h -- Uspensky's root isolation header file.

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

#ifndef _USP_H_
#define _USP_H_

#include <gmp.h>

typedef struct 
{
    mpz_t c; 
    long k; 
    unsigned int isexact; 
    char polysigns;
} interval; 

interval * Uspensky(mpz_t *, unsigned long, unsigned int *);
interval * Uspensky_couple(mpz_t *, mpz_t *, unsigned long, unsigned long, 
			   unsigned int *);

void X2XP1(mpz_t * P, unsigned long deg);
void X2XPC(mpz_t * P, mpz_t c, unsigned long deg);
void X2XM1(mpz_t * P, unsigned long deg);
int Homoth(mpz_t * P, long k, unsigned long deg);

#endif /* _USP_H_ */

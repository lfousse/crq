/* crq_structs.h -- structures definitions.

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
#ifndef _CRQ_STRUCTS_H
#define _CRQ_STRUCTS_H

#include <stdio.h> /* for FILE */

struct crq_stats {
    int pol_time; /* Time spent computing the polynomial */
    int usp_time; /* Time spent in Uspensky */
    int refine_time; /* Time spent refining */
    int coeff_time; /* Time spent computing the coefficients */
    int newton_good; /* Time spent in successfull newton iterations */
    int newton_bad; /* Time spent in useless newton iterations */
    int dich_prep_time; /* Time spent preparing the dichotomy */
    int dicho_time; /* Time spent in the dichotomy */
    int refinement_failures; /* Number of times the refinement iteration failed */
};

typedef enum {
    FULL,   /* full dichotomy with integer operations */
    TRUNC   /* truncated dichotomy */
} dicho_opts;

typedef enum {
    SIMPLE,	/* isolate roots of p only */
    DOUBLE,	/* isolate roots of p and p' */
    FROMFILE	/* load from file */
} isol_opts;

typedef enum {
    NEWTON,	/* newton root refinement */
    SECANTE,	/* secante root refinement */
    NEWTON_INT,	/* interval newton root refinement */
    DICHO	/* use *only* dichotomy for root refinement */
} refine_opts;

struct crq_gl_opts {
    dicho_opts dicho;
    isol_opts iso;
    refine_opts ref;
    FILE *file;
};

struct crq_gl_root {
    mpfr_t *xmin;
    mpfr_t *xmax;   /* roots \in [xmin, xmax] */
    mpfr_t *approx; /* curr approx */
};

/* Integration data types */
typedef void (aqfunc_t)(mpfr_t, mpfr_t);
typedef void (aqderiv_bound_t)(mpfr_t, mpfr_t, mpfr_t, unsigned int);

#endif  /* _CRQ_STRUCTS_H */

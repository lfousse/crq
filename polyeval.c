/* Implements algorithm polyeval and remainder tree using middle product.

  Copyright 2003, 2004, 2005 Laurent Fousse, Alexander Kruppa, Paul Zimmermann.

  This file is part of the ECM Library.

  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
  MA 02111-1307, USA.
*/

#include <stdlib.h>
#include <string.h> /* for strlen */
#include <stdio.h>
#include <gmp.h>
#include "listz.h"

unsigned int
TMulGen (listz_t b, unsigned int n, listz_t a, unsigned int m,
         listz_t c, unsigned int l, listz_t tmp, mpz_t modulus);
unsigned int
TMulGen_space (unsigned int n, unsigned int m, unsigned int l);

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#define ECM_STDOUT stdout

#ifndef NDEBUG
void
print_vect (listz_t t, unsigned int l)
{
    unsigned int i;

    fprintf (ECM_STDOUT, "[");
    for (i = 0; i < l; i++)
    {
        mpz_out_str (ECM_STDOUT, 10, t[i]);
        if (i != l - 1)
            fprintf (ECM_STDOUT, ", ");
        else
            fprintf (ECM_STDOUT, "]");
    }
}
#endif

/* Computes TUpTree as described in ref[1]. k is the degree of the
 * polynomial at the root of the tree. sh is the shift we need to
 * apply to find the actual coefficients of the polynomial at the root
 * of the tree.
 */

void
TUpTree (listz_t b, listz_t *Tree, unsigned int k, listz_t tmp, int dolvl,
         unsigned int sh, mpz_t n)
{
    unsigned int m, l;

    m = k / 2;
    l = k - m;
    
    if (k == 1)
      return;
   
#if 0
    fprintf (ECM_STDOUT, "In TupTree, k = %d.\n", k);

    fprintf (ECM_STDOUT, "b = ");
    print_vect (b, k);
    fprintf (ECM_STDOUT, "\nThe polynomials at that level are: ");
    print_vect (Tree[0] + sh, k);
    fprintf (ECM_STDOUT, "\n");
#endif

    if (dolvl == 0 || dolvl == -1)
      {
#ifdef DEBUG_TREEDATA
            printf ("Got from Tree: ");
            print_vect (Tree[0] + sh, l);
            print_vect (Tree[0] + sh + l, m);
            printf ("\n");
#endif
            TMulGen (tmp + l, m - 1, Tree[0] + sh, l - 1, b, k - 1, tmp + k, n);
            TMulGen (tmp, l - 1, Tree[0] + sh + l, m - 1, b, k - 1, tmp + k, n);

#if defined(DEBUG) || defined (DEBUG_TREEDATA)
        fprintf (ECM_STDOUT, "And the result at that level (before correction) is:");
        print_vect (tmp, k);
        fprintf (ECM_STDOUT, "\n");
#endif

        /* GMP-ECM specific: leading coefficients in the product tree
        * are implicit ones, so we need some extra work here.
        */

        list_add (tmp, tmp, b + m, l);
        list_add (tmp + l, tmp + l, b + l, m);

        list_mod (b, tmp, k, n); /* reduce both parts simultaneously */

#ifdef DEBUG
        fprintf (ECM_STDOUT, "And the result at this level is:");
        print_vect (b, k);
        fprintf (ECM_STDOUT, "\n");
#endif
      }
    
    if (dolvl > 0 || dolvl == -1)
      {
        if (dolvl > 0)
          dolvl--;
        TUpTree (b, Tree + 1, l, tmp, dolvl, sh, n);
        TUpTree (b + l, Tree + 1, m, tmp, dolvl, sh + l, n);
      }
}

static unsigned int
TUpTree_space (unsigned int k)
{

    unsigned int m, l;
    unsigned int r1, r2;

    m = k / 2;
    l = k - m;
    
    if (k == 1)
      return 0;
   
    r1 = TMulGen_space (l - 1, m - 1, k - 1) + l;
    if (m != l)
      {
        r2 = TMulGen_space (m - 1, l - 1, k - 1) + k;
        r1 = MAX (r1, r2);
      }

    r2 = TUpTree_space (l);
    r1 = MAX (r1, r2);
    
    if (m != l)
      {
        r2 = TUpTree_space (m);
        r1 = MAX (r1, r2);
      }

    return r1;
}

/* Same as polyeval. Needs invF as extra argument.
   Return non-zero iff an error occurred.
*/
int
polyeval_tellegen (listz_t b, unsigned int k, listz_t *Tree, listz_t tmp,
                   unsigned int sizeT, listz_t invF, mpz_t n)
{
    unsigned int tupspace;
    unsigned int tkspace;
    int allocated = 0, 
        r = 0; /* return value, 0 = no error */
    listz_t T;

    tupspace = TUpTree_space (k) + k;
#ifndef USE_SHORT_PRODUCT
    tkspace = TMulGen_space (k - 1, k - 1, k - 1) + k;
#else
    tkspace = 2 * k - 1 + list_mul_mem (k);
#endif

    tupspace = MAX (tupspace, tkspace);
    
    if (sizeT >= tupspace)
        T = tmp;
    else
      {
        T = init_list (tupspace);
        allocated = 1;
      }
    
#ifdef TELLEGEN_DEBUG
    fprintf (ECM_STDOUT, "In polyeval_tellegen, k = %d.\n", k);
    fprintf (ECM_STDOUT, "Required memory: %d.\n", 
	     TMulGen_space (k - 1, k - 1, k - 1));
#endif

    /* revert invF for call to TMulGen below */
    list_revert (invF, k - 1);
#ifndef NDEBUG
    printf ("invF = ");
    pari_poly (invF, k - 1, 0);
    printf ("\n");
#endif
    TMulGen (T, k - 1, invF, k - 1, b, k - 1, T + k, n);

    list_revert (T, k - 1);
#ifndef NDEBUG
    printf ("T = ");
    pari_poly (T, k - 1, 0);
    printf ("\n");
#endif
    TUpTree (T, Tree, k, T + k, -1, 0, n);
    list_swap (b, T, k); /* more efficient than list_set, since T is not
                            needed anymore */

    if (allocated)
      clear_list (T, tupspace);

    return r;
}

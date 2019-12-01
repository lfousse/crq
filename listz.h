/* Header file for listz operations.

  Copyright 2006 Laurent Fousse

  This file is part of the CRQ Library.

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

#ifndef _LISTZ_H
typedef mpz_t *listz_t;
void list_set (listz_t p, listz_t q, unsigned int n);
void list_zero (listz_t p, unsigned int n);
void list_add (listz_t p, listz_t q, listz_t r, unsigned int l);
void list_sub (listz_t p, listz_t q, listz_t r, unsigned int l);
void list_mod (listz_t p, listz_t q, unsigned int n, mpz_t mod);
listz_t init_list (unsigned int n);
void list_swap (listz_t p, listz_t q, unsigned int n);
void list_revert (listz_t p, unsigned int n);
void clear_list (listz_t p, unsigned int n);
void PolyInvert (listz_t q, listz_t b, unsigned int K, listz_t t, mpz_t n);
#endif /* _LISTZ_H */

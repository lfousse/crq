
#include <stdio.h>
#include <gmp.h>

/* Given u = x/y mod n
 * where |x| < B and |y| < B
 * and n >= 2B^2, compute a and b
 * such that a/b = x/y.
 */
void frac_rec (mpz_t a, mpz_t b, mpz_t u, mpz_t n, mpz_t B)
{
    mpz_t r[2];
    mpz_t t[2];
    mpz_t s[2];
    mpz_t q, tmp;
    int p;

    mpz_init (r[0]);
    mpz_init (r[1]);
    mpz_init (t[0]);
    mpz_init (t[1]);
    mpz_init (s[0]);
    mpz_init (s[1]);
    mpz_init (q);
    mpz_init (tmp);

    mpz_set (r[0], n);
    mpz_set (r[1], u);
    mpz_set_ui (t[0], 0);
    mpz_set_ui (t[1], 1);
    mpz_set_ui (s[0], 1);
    mpz_set_ui (s[1], 0);
    p = 0;

    for (;;) {
#ifdef DEBUG
	mpz_out_str (stdout, 10, s[p]);
	printf (" * ");
	mpz_out_str (stdout, 10, n);
	printf (" + ");
	mpz_out_str (stdout, 10, t[p]);
	printf (" * ");
	mpz_out_str (stdout, 10, u);
	printf (" = ");
	mpz_out_str (stdout, 10, r[p]);
	printf ("\n");
#endif
	if (mpz_cmp (r[p], B) < 0) {
#if 0
	    /* TODO : calculer l'autre possible solution */
	    {
		mpz_t foo, bar;
		mpz_init (foo);
		mpz_set (foo, r[p]);
		mpz_mul (foo, foo, q);
		mpz_sub (foo, r[1 - p], foo);

		mpz_init (bar);
		mpz_set (bar, t[p]);
		mpz_mul (bar, bar, q);
		mpz_sub (bar, t[1 - p], bar);
		printf ("Autre solution : "); mpz_out_str (stdout, 10, foo);
		printf (" / "); mpz_out_str (stdout, 10, bar); printf ("\n");
	    }
#endif
	    mpz_set (a, r[p]);
	    mpz_set (b, t[p]);
	    break;
	}
	mpz_tdiv_qr (q, r[p], r[p], r[1 - p]);
#ifdef DEBUG
	mpz_out_str (stdout, 10, q); printf ("\n");
#endif
	mpz_mul (tmp, q, t[1 - p]);
	mpz_sub (t[p], t[p], tmp);
	mpz_mul (tmp, q, s[1 - p]);
	mpz_sub (s[p], s[p], tmp);
	p = 1 - p;
    }
    mpz_clear (r[0]);
    mpz_clear (r[1]);
    mpz_clear (t[0]);
    mpz_clear (t[1]);
    mpz_clear (q);
}

#if 0
int main (void)
{
    mpz_t a, b, u, n, B;

    mpz_init (a);
    mpz_init (b);
    mpz_init (u);
    mpz_init (n);
    mpz_init (B);

    mpz_set_ui (u, 3645263);
    mpz_set_ui (B, 2000);
    mpz_set_ui (n, 8000000);
    /* 1001 / 327 */
    frac_rec (a, b, u, n, B);
    mpz_out_str (stdout, 10, a); printf ("\n");
    mpz_out_str (stdout, 10, b); printf ("\n");
    return 0;
}
#endif

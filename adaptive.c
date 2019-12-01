#include "crq.h"
#include "utils.h"

void adaptive_quad_aux (mpfr_t res, mpfr_t global, mpfr_t a, mpfr_t b, aqfunc_t f, mp_prec_t prec,
		        void (*quadf)(mpfr_t, mpfr_t, mpfr_t, aqfunc_t, void *),
		        void *params, mpfr_t eps)
{
    mpfr_t medium, leftv, rightv;

    mpfr_init2 (medium, prec + 10);
    mpfr_init2 (leftv, prec + 10);
    mpfr_init2 (rightv, prec + 10);

    mpfr_add (medium, a, b, GMP_RNDN); /* exact */
    mpfr_div_2ui (medium, medium, 1, GMP_RNDN);

    quadf (leftv, a, medium, f, params);
    quadf (rightv, medium, b, f, params);
    mpfr_add (res, leftv, rightv, GMP_RNDN);

    mpfr_sub (global, res, global, GMP_RNDN);
    mpfr_abs (global, global, GMP_RNDN);
    if (mpfr_cmp (global, eps) > 0) {
#if VERBOSE >= 3
	printf ("diff = ");
	mpfr_dump (global);
	printf ("eps = ");
	mpfr_dump (eps);
#endif
	adaptive_quad_aux (res, leftv, a, medium, f, prec, quadf, params, eps);
	adaptive_quad_aux (global, rightv, medium, b, f, prec, quadf, params, eps);
	mpfr_add (res, res, global, GMP_RNDN);
    }
    mpfr_clears (medium, leftv, rightv, (void *) 0);
}

/* Adaptive quadrature of f on [a, b], using */
void adaptive_quad (mpfr_t res, mpfr_t a, mpfr_t b, aqfunc_t f, mp_prec_t prec,
		    void (*quadf)(mpfr_t, mpfr_t, mpfr_t, aqfunc_t, void *),
		    void *params)
{
    mpfr_t global, eps;

    mpfr_init2 (global, prec);
    mpfr_init (eps);

    /* rough estimate of full integral */
    quadf (global, a, b, f, params);
    
    /* Precision threshold */
    mpfr_set_ui (eps, 1, GMP_RNDN);
    mpfr_mul_2si (eps, eps, -prec - mpfr_get_exp (global), GMP_RNDN);

    adaptive_quad_aux (res, global, a, b, f, prec, quadf, params, eps);
    mpfr_clear (global);
    mpfr_clear (eps);
}

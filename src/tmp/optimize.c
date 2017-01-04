#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define NO_NLS
#include <Defn.h>
#include <float.h>		/* for DBL_MAX */
#include <R_ext/Applic.h>	/* for optif9, fdhess */
#include <R_ext/RS.h>	       	/* for Memcpy */

#include "statsR.h"
#include "stats.h" // R_zeroin2

#undef _
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif


/* Formerly in src/appl/fmim.c */

/* fmin.f -- translated by f2c (version 19990503).
*/

/* R's  optimize() :   function	fmin(ax,bx,f,tol)
   =    ==========		~~~~~~~~~~~~~~~~~

        an approximation  x  to the point where  f  attains a minimum  on
    the interval  (ax,bx)  is determined.

    INPUT..

    ax    left endpoint of initial interval
    bx    right endpoint of initial interval
    f     function which evaluates  f(x, info)  for any  x
          in the interval  (ax,bx)
    tol   desired length of the interval of uncertainty of the final
          result ( >= 0.)

    OUTPUT..

    fmin  abcissa approximating the point where  f  attains a minimum

        The method used is a combination of  golden  section  search  and
    successive parabolic interpolation.  convergence is never much slower
    than  that  for  a  Fibonacci search.  If  f  has a continuous second
    derivative which is positive at the minimum (which is not  at  ax  or
    bx),  then  convergence  is  superlinear, and usually of the order of
    about  1.324....
        The function  f  is never evaluated at two points closer together
    than  eps*abs(fmin)+(tol/3), where eps is  approximately  the  square
    root  of  the  relative  machine  precision.   if   f   is a unimodal
    function and the computed values of   f   are  always  unimodal  when
    separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
    the abcissa of the global minimum of  f  on the interval  ax,bx  with
    an error less than  3*eps*abs(fmin)+tol.  if   f   is  not  unimodal,
    then fmin may approximate a local, but perhaps non-global, minimum to
    the same accuracy.
        This function subprogram is a slightly modified  version  of  the
    Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
    Minimization without Derivatives, Prentice-Hall, Inc. (1973).
*/
#include <math.h>
#include <float.h> /* DBL_EPSILON */

#include <Rmath.h>
#include <R_ext/Applic.h>

static
double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		  void *info, double tol) {
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = (*f)(x, info);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */

	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */

	    d = p / q;
	    u = x + d;

	    /* f must not be evaluated too close to ax or bx */

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = (*f)(u, info);

	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    return x;
}


/* One Dimensional Minimization --- just wrapper for
 * Brent's "fmin" --> ../appl/fmin.c */

struct callinfo {
  SEXP R_fcall;
  SEXP R_env;
} ;

/*static SEXP R_fcall1;
  static SEXP R_env1; */

static double fcn1(double x, struct callinfo *info)
{
  SEXP s, sx;
  PROTECT(sx = ScalarReal(x));
  SETCADR(info->R_fcall, sx);
  s = eval(info->R_fcall, info->R_env);
  UNPROTECT(1);
  switch(TYPEOF(s)) {
    case INTSXP:
    	if (length(s) != 1) goto badvalue;
    	if (INTEGER(s)[0] == NA_INTEGER) {
        warning(_("NA replaced by maximum positive value"));
        return DBL_MAX;
    	} else {
        return INTEGER(s)[0];
      }
    	break;
    case REALSXP:
    	if (length(s) != 1) goto badvalue;
    	if (!R_FINITE(REAL(s)[0])) {
    	    warning(_("NA/Inf replaced by maximum positive value"));
    	    return DBL_MAX;
    	} else {
        return REAL(s)[0];
      }
    	break;
    default:
    	goto badvalue;
  }
 badvalue:
    error(_("invalid function value in 'optimize'"));
    return 0;/* for -Wall */
}

/* fmin(f, xmin, xmax tol) */
SEXP do_fmin(SEXP call, SEXP op, SEXP args, SEXP rho) {
  double xmin, xmax, tol;
  SEXP v, res;
  struct callinfo info;

  args = CDR(args);
  PrintDefaults();

  /* the function to be minimized */

  v = CAR(args);
  if (!isFunction(v)) {
    error(_("attempt to minimize non-function"));
  }
  args = CDR(args);

  /* xmin */

  xmin = asReal(CAR(args));
  if (!R_FINITE(xmin)) {
    error(_("invalid '%s' value"), "xmin");
  }
  args = CDR(args);

  /* xmax */

  xmax = asReal(CAR(args));
  if (!R_FINITE(xmax)) {
  	error(_("invalid '%s' value"), "xmax");
  }
  if (xmin >= xmax) {
  	error(_("'xmin' not less than 'xmax'"));
  }
  args = CDR(args);

  /* tol */

  tol = asReal(CAR(args));
  if (!R_FINITE(tol) || tol <= 0.0) {
  	error(_("invalid '%s' value"), "tol");
  }

  info.R_env = rho;
  PROTECT(info.R_fcall = lang2(v, R_NilValue));
  PROTECT(res = allocVector(REALSXP, 1));
  REAL(res)[0] = Brent_fmin(xmin, xmax,
  	      (double (*)(double, void*)) fcn1, &info, tol);
  UNPROTECT(2);
  return res;
}

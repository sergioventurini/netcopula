#include "netcopula.h"

// [[Rcpp::depends("RcppArmadillo")]]

// TO DO:
//   1. re-implement these functions using a pointer to function instead of
//      Rcpp::Function objects (as in the optim.c code in R stats package)

#define big             1.0e+35   /*a very large number*/

static double fminfn(const Rcpp::NumericVector& p, const Rcpp::List& ex) {
  double fnscale, s, val;
  Rcpp::Function fn = Rcpp::as<Rcpp::Function>(ex["fcall"]);
  Rcpp::List args = Rcpp::as<Rcpp::List>(ex["args"]);

  fnscale = Rcpp::as<double>(ex["fnscale"]);
  s = Rcpp::as<double>(fn(p, args));
  val = s/fnscale;

  return val;
}

static void fmingr(const int& n, const Rcpp::NumericVector& p, Rcpp::NumericVector df, const Rcpp::List& OS) {
  int i;
  double fnscale, s, val1, val2, eps;

  Rcpp::NumericVector x(n);
  Rcpp::NumericVector parscale = Rcpp::as<Rcpp::NumericVector>(OS["parscale"]);
  Rcpp::NumericVector ndeps = Rcpp::as<Rcpp::NumericVector>(OS["ndeps"]);
  Rcpp::Function fn = Rcpp::as<Rcpp::Function>(OS["fcall"]);
  Rcpp::List args = Rcpp::as<Rcpp::List>(OS["args"]);
  fnscale = Rcpp::as<double>(OS["fnscale"]);

  /* numerical derivatives */
  for (i = 0; i < n; i++) {
    if (!R_FINITE(p(i))) {
      error("non-finite value supplied by optim_rcpp");
    }
    x(i) = p(i) * parscale(i);
  }
  for (i = 0; i < n; i++) {
    eps = ndeps(i);
    x(i) = (p(i) + eps) * parscale(i);
    s = Rcpp::as<double>(fn(x, args));
    val1 = s/fnscale;
    x(i) = (p(i) - eps) * parscale(i);
    s = Rcpp::as<double>(fn(x, args));
    val2 = s/fnscale;
    df(i) = (val1 - val2)/(2 * eps);
    // if (!R_FINITE(df(i))) {
    //   error(("non-finite finite-difference value [%d]"), i + 1);
    // }
    x(i) = p(i) * parscale(i);
  }
}

void nmmin_rcpp(int n, Rcpp::NumericVector Bvec, Rcpp::NumericVector X, double *Fmin, int *fail, double abstol, double intol, const Rcpp::List& ex, double alpha, double bet, double gamm, int trace, int *fncount, int maxit) {
  char action[50];
  int C;
  Rboolean calcvert;
  double convtol, f;
  int funcount = 0, H, i, j, L = 0;
  int n1 = 0;
  double oldsize;
  double size, step, temp, trystep;
  char tstr[9]; // allow for 10^8 iters ...
  double VH, VL, VR;

  if (maxit <= 0) {
    *Fmin = fminfn(Bvec, ex);
    *fncount = 0;
    *fail = 0;
    return;
  }
  if (trace) {
    Rprintf("  Nelder-Mead direct search function minimizer\n");
  }
  Rcpp::NumericMatrix P(n + 1, (n + 1) + 1);
  *fail = FALSE;
  f = fminfn(Bvec, ex);
  if (!R_FINITE(f)) {
    // error("function cannot be evaluated at initial parameters");
    *fail = TRUE;
    return;   // added by myself!
  } else {
    if (trace) Rprintf("function value for initial parameters = %f\n", f);
    funcount = 1;
    convtol = intol * (fabs(f) + intol);
    if (trace) Rprintf("  Scaled convergence tolerance is %g\n", convtol);
    n1 = n + 1;
    C = n + 2;
    P(n1 - 1, 0) = f;
    for (i = 0; i < n; i++) {
      P(i, 0) = Bvec(i);
    }

    L = 1;
    size = 0.0;

    step = 0.0;
    for (i = 0; i < n; i++) {
      if (0.1 * fabs(Bvec(i)) > step) {
        step = 0.1 * fabs(Bvec(i));
      }
    }
    if (step == 0.0) step = 0.1;
    if (trace) Rprintf("Stepsize computed as %f\n", step);
    for (j = 2; j <= n1; j++) {
      strcpy(action, "BUILD          ");
      for (i = 0; i < n; i++) {
        P(i, j - 1) = Bvec(i);
      }

      trystep = step;
      while (P(j - 2, j - 1) == Bvec(j - 2)) {
        P(j - 2, j - 1) = Bvec(j - 2) + trystep;
        trystep *= 10;
      }
      size += trystep;
    }
    oldsize = size;
    calcvert = TRUE;
    do {
      if (calcvert) {
        for (j = 0; j < n1; j++) {
          if (j + 1 != L) {
            for (i = 0; i < n; i++) {
              Bvec(i) = P(i, j);
            }
            f = fminfn(Bvec, ex);
            if (!R_FINITE(f)) f = big;
            funcount++;
            P(n1 - 1, j) = f;
          }
        }
        calcvert = FALSE;
      }

      VL = P(n1 - 1, L - 1);
      VH = VL;
      H = L;

      for (j = 1; j <= n1; j++) {
        if (j != L) {
          f = P(n1 - 1, j - 1);
          if (f < VL) {
            L = j;
            VL = f;
          }
          if (f > VH) {
            H = j;
            VH = f;
          }
        }
      }

      if (VH <= VL + convtol || VL <= abstol) break;

      // avoid buffer overflow at 100001 iters. (PR#15240)
      if (trace) {
        snprintf(tstr, 9, "%5d", funcount);
        Rprintf("%s%s %f %f\n", action, tstr, VH, VL);
      }

      for (i = 0; i < n; i++) {
        temp = -P(i, H - 1);
        for (j = 0; j < n1; j++) {
          temp += P(i, j);
        }
        P(i, C - 1) = temp / n;
      }
      for (i = 0; i < n; i++) {
        Bvec(i) = (1.0 + alpha) * P(i, C - 1) - alpha * P(i, H - 1);
      }
      f = fminfn(Bvec, ex);
      if (!R_FINITE(f)) f = big;
      funcount++;
      strcpy(action, "REFLECTION     ");
      VR = f;
      if (VR < VL) {
        P(n1 - 1, C - 1) = f;
        for (i = 0; i < n; i++) {
          f = gamm * Bvec(i) + (1 - gamm) * P(i, C - 1);
          P(i, C - 1) = Bvec(i);
          Bvec(i) = f;
        }
        f = fminfn(Bvec, ex);
        if (!R_FINITE(f)) f = big;
        funcount++;
        if (f < VR) {
          for (i = 0; i < n; i++)
          P(i, H - 1) = Bvec(i);
          P(n1 - 1, H - 1) = f;
          strcpy(action, "EXTENSION      ");
        } else {
          for (i = 0; i < n; i++) {
            P(i, H - 1) = P(i, C - 1);
          }
          P(n1 - 1, H - 1) = VR;
        }
      } else {
        strcpy(action, "HI-REDUCTION   ");
        if (VR < VH) {
          for (i = 0; i < n; i++) {
            P(i, H - 1) = Bvec(i);
          }
          P(n1 - 1, H - 1) = VR;
          strcpy(action, "LO-REDUCTION   ");
        }

        for (i = 0; i < n; i++) {
          Bvec(i) = (1 - bet) * P(i, H - 1) + bet * P(i, C - 1);
        }
        f = fminfn(Bvec, ex);
        if (!R_FINITE(f)) f = big;
        funcount++;

        if (f < P(n1 - 1, H - 1)) {
          for (i = 0; i < n; i++) {
            P(i, H - 1) = Bvec(i);
          }
          P(n1 - 1, H - 1) = f;
        } else {
          if (VR >= VH) {
            strcpy(action, "SHRINK         ");
            calcvert = TRUE;
            size = 0.0;
            for (j = 0; j < n1; j++) {
              if (j + 1 != L) {
                for (i = 0; i < n; i++) {
                  P(i, j) = bet * (P(i, j) - P(i, L - 1)) + P(i, L - 1);
                  size += fabs(P(i, j) - P(i, L - 1));
                }
              }
            }
            if (size < oldsize) {
              oldsize = size;
            } else {
              if (trace)
              Rprintf("Polytope size measure not decreased in shrink\n");
              *fail = 10;
              break;
            }
          }
        }
      }
    } while (funcount <= maxit);
  }

  if (trace) {
    Rprintf("Exiting from Nelder Mead minimizer\n");
    Rprintf("    %d function evaluations used\n", funcount);
  }
  *Fmin = P(n1 - 1, L - 1);
  for (i = 0; i < n; i++) {
    X(i) = P(i, L - 1);
  }
  if (funcount > maxit) {
    *fail = 1;
  }
  *fncount = funcount;
}

//' @export
// [[Rcpp::export]]
Rcpp::List optim_rcpp(const Rcpp::NumericVector& par, const Rcpp::Function& fn, const Rcpp::List& args, const Rcpp::List& options, const bool& hessian = true) {
  int i, npar = 0, trace, maxit, fncount = 0, ifail = 0;
  double fnscale, val = 0.0, abstol, reltol, alpha, beta, gamm;

  Rcpp::List OS = Rcpp::List::create(
      Rcpp::Named("fcall") = fn,
      Rcpp::Named("args") = args,
      Rcpp::Named("ndeps") = Rcpp::as<Rcpp::NumericVector>(options["ndeps"]),
      Rcpp::Named("fnscale") = Rcpp::as<double>(options["fnscale"]),
      Rcpp::Named("parscale") = Rcpp::as<Rcpp::NumericVector>(options["parscale"]));

  npar = par.size();
  Rcpp::NumericVector dpar(npar), opar(npar);
  fnscale = Rcpp::as<double>(options["fnscale"]);
  trace = Rcpp::as<int>(options["trace"]);
  Rcpp::NumericVector parscale = Rcpp::as<Rcpp::NumericVector>(options["parscale"]);
  for (i = 0; i < npar; i++) {
    dpar(i) = par(i) / parscale(i);
  }
  abstol = Rcpp::as<double>(options["abstol"]);
  reltol = Rcpp::as<double>(options["reltol"]);
  maxit = Rcpp::as<double>(options["maxit"]);

  alpha = Rcpp::as<double>(options["alpha"]);
  beta = Rcpp::as<double>(options["beta"]);
  gamm = Rcpp::as<double>(options["gamma"]);
  nmmin_rcpp(npar, dpar, opar, &val, &ifail, abstol, reltol, OS, alpha, beta, gamm, trace, &fncount, maxit);
  for (i = 0; i < npar; i++) {
    opar(i) = opar(i) * parscale(i);
  }

  if (hessian) {
    Rcpp::NumericMatrix hess = optimhess_rcpp(opar, fn, args, options);
    return Rcpp::List::create(Rcpp::Named("par") = opar,
                              Rcpp::Named("value") = val * fnscale,
                              Rcpp::Named("counts") = fncount,
                              Rcpp::Named("convergence") = ifail,
                              Rcpp::Named("hessian") = hess);
  } else {
    return Rcpp::List::create(Rcpp::Named("par") = opar,
                              Rcpp::Named("value") = val * fnscale,
                              Rcpp::Named("counts") = fncount,
                              Rcpp::Named("convergence") = ifail,
                              Rcpp::Named("hessian") = false);
  }
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix optimhess_rcpp(const Rcpp::NumericVector& par, const Rcpp::Function& fn, const Rcpp::List& args, const Rcpp::List& options) {
  int i, j, npar = 0;
  double fnscale, eps;

  Rcpp::List OS = Rcpp::List::create(
      Rcpp::Named("fcall") = fn,
      Rcpp::Named("args") = args,
      Rcpp::Named("ndeps") = Rcpp::as<Rcpp::NumericVector>(options["ndeps"]),
      Rcpp::Named("fnscale") = Rcpp::as<double>(options["fnscale"]),
      Rcpp::Named("parscale") = Rcpp::as<Rcpp::NumericVector>(options["parscale"]));

  npar = par.size();
  fnscale = Rcpp::as<double>(options["fnscale"]);
  Rcpp::NumericVector parscale = Rcpp::as<Rcpp::NumericVector>(options["parscale"]);
  Rcpp::NumericVector ndeps = Rcpp::as<Rcpp::NumericVector>(options["ndeps"]);
  if (ndeps.size() != npar) {
    error("'ndeps' is of the wrong length");
  }
  Rcpp::NumericMatrix ans(npar, npar);
  Rcpp::NumericVector dpar(npar);
  for (i = 0; i < npar; i++) {
    dpar(i) = par(i) / parscale(i);
  }
  Rcpp::NumericVector df1(npar), df2(npar);
  for (i = 0; i < npar; i++) {
    eps = ndeps(i) / parscale(i);
    dpar(i) = dpar(i) + eps;
    fmingr(npar, dpar, df1, OS);
    dpar(i) = dpar(i) - 2 * eps;
    fmingr(npar, dpar, df2, OS);
    for (j = 0; j < npar; j++) {
      ans(i, j) = fnscale * (df1(j) - df2(j))/(2 * eps * parscale(i) * parscale(j));
    }
    dpar(i) = dpar(i) + eps;
  }
  // symmetrize the hessian
  for (i = 0; i < npar; i++) {
    for (j = 0; j < i; j++) {
      double tmp = 0.5 * (ans(i, j) + ans(j, i));
      ans(i, j) = ans(j, i) = tmp;
    }
  }

  return ans;
}

//' @export
// [[Rcpp::export]]
Rcpp::List laplace_rcpp(const Rcpp::Function& logpost, const Rcpp::NumericVector& guess, const Rcpp::List& args, const Rcpp::List& options) {
  double logdet, logdet_sgn, val = 0.0, s = 0.0, eps = 1e-6;
  int p = guess.size();
  Rcpp::List fit = optim_rcpp(guess, logpost, args, options, true);
  Rcpp::NumericVector mode(p);
  arma::mat var(p, p);
  var = eps*arma::eye(p, p);
  bool conv = false;

  if (!Rcpp::as<int>(fit["convergence"])) {
    mode = Rcpp::as<Rcpp::NumericVector>(fit["par"]);
    arma::mat hessian_orig = Rcpp::as<arma::mat>(fit["hessian"]);
    // ***
    // these lines are needed to avoid errors in getting the inverse hessian
    // when the hessian is singular
    arma::mat hessian = hessian_orig;
    while (is_singular(hessian)) {
      Rprintf("*** singular hessian ==> eps = %4.1e ***\n", eps);
      hessian = hessian_orig + eps*arma::eye(p, p);
      eps *= 10;
    }
    // ***
    var = -arma::inv(hessian);
    conv = (Rcpp::as<int>(fit["convergence"]) == 0 ? true : false);
    arma::log_det(logdet, logdet_sgn, var);
    s = Rcpp::as<double>(logpost(mode, args));
    val = p/2.0 * log(2.0 * M_PI) + 0.5 * logdet + s;
  }

  return Rcpp::List::create(Rcpp::Named("mode") = mode,
                            Rcpp::Named("var") = var,
                            Rcpp::Named("int") = val,
                            Rcpp::Named("converge") = conv);
}

//' @export
// [[Rcpp::export]]
Rcpp::List optim_rcpp_example() {
  Rcpp::NumericVector par(2), parscale(2), ndeps(2);
  par(0) = -7.0;
  par(1) = 6.0;
  parscale(0) = 1;
  parscale(1) = 1;
  ndeps(0) = 0.001;
  ndeps(1) = 0.001;
  Rcpp::List args = Rcpp::List::create(
      Rcpp::Named("tau") = 1.0,
      Rcpp::Named("eta") = 2.0,
      Rcpp::Named("y") = 12,
      Rcpp::Named("n") = 20,
      Rcpp::Named("w") = 1.0,
      Rcpp::Named("gamma") = 1.0,
      Rcpp::Named("mu_sigma") = 10.0,
      Rcpp::Named("eps") = 1e-12,
      Rcpp::Named("eps_abs") = 1e-3);
  Rcpp::List options = Rcpp::List::create(
      Rcpp::Named("trace") = 0,
      Rcpp::Named("fnscale") = -1,
      Rcpp::Named("parscale") = Rcpp::as<Rcpp::NumericVector>(parscale),
      Rcpp::Named("ndeps") = Rcpp::as<Rcpp::NumericVector>(ndeps),
      Rcpp::Named("maxit") = 500,
      Rcpp::Named("abstol") = R_NegInf,
      Rcpp::Named("reltol") = 1.490116e-08,
      Rcpp::Named("alpha") = 1.0,
      Rcpp::Named("beta") = 0.5,
      Rcpp::Named("gamma") = 2.0);
  Rcpp::Function mu_delta_logpost("mu_delta_logpost_func");
  Rcpp::List res = optim_rcpp(par, mu_delta_logpost, args, options, true);

  return res;
}

static double Brent_fmax(double ax, double bx, Rcpp::Function fn, Rcpp::List args, double tol) {
  /*  c is the squared inverse of the golden ratio */
  const double c = (3. - sqrt(5.)) * .5;

  /* Local variables */
  double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

  /*  eps is approximately the square root of the relative machine precision. */
  eps = DBL_EPSILON;
  tol1 = eps + 1.; /* the smallest 1.000... > 1 */
  eps = sqrt(eps);

  a = ax;
  b = bx;
  v = a + c * (b - a);
  w = v;
  x = v;

  d = 0.;/* -Wall */
  e = 0.;
  /* -1.* is needed because the original Brent_fmin minimizes a function */
  fx = -1.*Rcpp::as<double>(fn(x, args)); 
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
    } else { /* a parabolic-interpolation step */
      d = p / q;
      u = x + d;

      /* f must not be evaluated too close to ax or bx */

      if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm) d = -d;
      }
    }

    /* f must not be evaluated too close to x */

    if (fabs(d) >= tol1) {
      u = x + d;
    } else if (d > 0.) {
      u = x + tol1;
    } else {
      u = x - tol1;
    }

    /* -1.* is needed because the original Brent_fmin minimizes a function */
    fu = -1.*Rcpp::as<double>(fn(u, args));

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

//' @export
// [[Rcpp::export]]
double optimize_rcpp(const Rcpp::Function& fn, const double& xmin, const double& xmax, const double& tol, const Rcpp::List& args) {
  if (!R_FINITE(xmin)) {
    error("invalid 'xmin' value");
  }
  if (!R_FINITE(xmax)) {
    error("invalid 'xmax' value");
  }
  if (xmin >= xmax) {
    error("'xmin' cannot be less than 'xmax'");
  }
  if (!R_FINITE(tol) || tol <= 0.0) {
    error("invalid 'tol' value");
  }

  return Brent_fmax(xmin, xmax, fn, args, tol);
}

//' @export
// [[Rcpp::export]]
Rcpp::List laplace_u_rcpp(const Rcpp::Function& fn, const double& xmin, const double& xmax, const double& tol, const Rcpp::List& args, const Rcpp::List& options) {
  Rcpp::NumericVector mode(1);
  mode(0) = optimize_rcpp(fn, xmin, xmax, tol, args);
  double var = -1./optimhess_rcpp(mode, fn, args, options)(0, 0);

  return Rcpp::List::create(Rcpp::Named("mode") = mode(0),
                            Rcpp::Named("var") = var);
}

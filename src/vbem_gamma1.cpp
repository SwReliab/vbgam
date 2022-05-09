#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

double dgamma_log(double t0, double t1, double alpha, double rate) {
  return alpha * log(rate) + (alpha-1) * t0 - rate * t1 - lgamma(alpha);
}

// @name vbem
// [[Rcpp::export]]

List vb1fit(double alpha, List params, List prior, List data, List options) {

  const int maxiter = options["maxiter"];
  const double abstol = options["abstol"];
  const double reltol = options["reltol"];

  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  const double mw = as<NumericVector>(prior["omega"])[0];
  const double phiw = as<NumericVector>(prior["omega"])[1];
  const double mb = as<NumericVector>(prior["beta"])[0];
  const double phib = as<NumericVector>(prior["beta"])[1];

  const double omega_shape = as<NumericVector>(params["omega"])[0];
  const double omega_rate = as<NumericVector>(params["omega"])[1];
  const double beta_shape = as<NumericVector>(params["beta"])[0];
  const double beta_rate = as<NumericVector>(params["beta"])[1];

  double psiw = -log(omega_rate) + R::digamma(omega_shape);
  double xiw = omega_shape / omega_rate;
  double psib = -log(beta_rate) + R::digamma(beta_shape);
  double xib = beta_shape / beta_rate;

  const double a0 = lgamma(alpha);
  const double a1 = lgamma(alpha + 1.0);


  double zetaN;
  double zetaT;
  double me;
  double logP;
  double vfe = NAN;
  double aerror;
  double rerror;

  int iter = 1;
  bool conv = false;

  while (1) {
    double prev = vfe;
    double t = 0.0;
    double gam10 = 1.0;
    double gam11 = 1.0;
    double gam20 = 1.0;
    double gam21 = 1.0;
    logP = 0.0;
    me = 0.0;
    zetaN = 0.0;
    zetaT = 0.0;

    for (int i=0; i<dsize; i++) {
      double x = num[i];
      if (time[i] != 0.0) {
        t += time[i];
        gam10 = gam20;
        gam11 = gam21;
        gam20 = R::pgamma(t, alpha, 1.0/xib, false, false);
        gam21 = R::pgamma(t, alpha+1, 1.0/xib, false, false);
        // gam20 = q_gamma(alpha, xib*t, a0);
        // gam21 = q_gamma(alpha+1, xib*t, a1);
      }
      if (x != 0.0) {
        double tmp1 = gam10 - gam20;
        double tmp2 = (alpha / xib) * (gam11 - gam21);
        me += x;
        zetaT += x * tmp2 / tmp1;
        logP += x * log(tmp1) - lgamma(x+1.0);
      }
      if (type[i] == 1) {
        me += 1.0;
        zetaT += t;
        logP += alpha*log(xib) + (alpha-1.0)*log(t) - xib*t - lgamma(alpha);
      }
    }
    zetaN = me + exp(psiw + psib * alpha - alpha*log(xib)) * gam20;
    zetaT += alpha * exp(psiw + psib * alpha - (alpha+1)*log(xib)) * gam21;
    logP += exp(psiw + psib * alpha - alpha*log(xib)) * gam20;

    psiw = -log(phiw + 1) + R::digamma(mw + zetaN);
    xiw = (mw + zetaN) / (phiw + 1);
    psib = -log(phib + zetaT) + R::digamma(mb + zetaN * alpha);
    xib = (mb + zetaN * alpha) / (phib + zetaT);

    vfe = -xiw + me * (psiw + psib * alpha - alpha * log(xib)) + logP
      + dgamma_log(psiw, xiw, mw, phiw)
      + dgamma_log(psib, xib, mb, phib)
      - dgamma_log(psiw, xiw, mw + zetaN, phiw + 1)
      - dgamma_log(psib, xib, mb + zetaN * alpha, phib + zetaT);

    aerror = fabs(prev - vfe);
    rerror = aerror / vfe;

    if (aerror < abstol && rerror < reltol) {
      conv = true;
      break;
    }

    if (iter >= maxiter) {
      break;
    }

    iter += 1;
  }

  List posterior = List::create(
    Named("omega") = NumericVector::create(mw+zetaN, phiw+1),
    Named("beta") = NumericVector::create(mb+zetaN*alpha, phib+zetaT)
  );

  return List::create(
    Named("posterior") = posterior,
    Named("vfe") = vfe,
    Named("zetaN") = zetaN,
    Named("zetaT") = zetaT,
    Named("iter") = iter,
    Named("conv") = conv,
    Named("aerror") = aerror,
    Named("rerror") = rerror
  );
}

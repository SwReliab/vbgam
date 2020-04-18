#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

List vbem_vb2gam(int n, double alpha, double init_zetaT, List prior, List data, List options) {

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
  const double a0 = lgamma(alpha);
  const double a1 = lgamma(alpha + 1.0);

  double zetaT = init_zetaT;
  double xib = (mb + n * alpha) / (phib + zetaT);
  double me;

//  Rcout << "  init zetaT=" << zetaT << " xib=" << xib << std::endl;

  int iter = 1;
  bool conv = false;
  double logP = 0.0;
  double logQ = 0.0;
  double aerror = 0.0;
  double rerror = 0.0;

  while (1) {
    double prev = zetaT;
    double t = 0.0;
    double gam10 = 1.0;
    double gam11 = 1.0;
    double gam20 = 1.0;
    double gam21 = 1.0;
    logP = 0.0;
    me = 0.0;
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
//        Rcout << t << " " << gam10 << " " << gam20 << " " << tmp1 << std::endl;
        logP += x * log(tmp1) - lgamma(x+1.0);
      }
      if (type[i] == 1) {
        me += 1.0;
        zetaT += t;
        logP += alpha*log(xib) + (alpha-1.0)*log(t) - xib*t - lgamma(alpha);
      }
    }
    zetaT += (n - me) * (alpha / xib) * gam21 / gam20;  // gam21 is the last time
    xib = (mb + n * alpha) / (phib + zetaT);
    logP += (n - me) * log(gam20) - lgamma(n-me+1.0);

//    Rcout << "  iter=" << iter << " zetaT=" << zetaT << " xib=" << xib << std::endl;

    aerror = fabs(prev - zetaT);
    rerror = aerror / zetaT;

    if (aerror < abstol && rerror < reltol) {
      conv = true;
      break;
    }

    if (iter >= maxiter) {
      break;
    }

    iter += 1;
  }

  if (mw == 0.0 || mb == 0.0) {
    logQ = lgamma(mw+n) + lgamma(mb+n*alpha) - (mw+n)*log(phiw+1) - (mb+n*alpha)*log(phib+zetaT)
      + logP + xib*zetaT - (n*alpha) * log(xib);
  } else {
    logQ = mw * log(phiw) + mb * log(phib)  - lgamma(mw) - lgamma(mb)
      + lgamma(mw+n) + lgamma(mb+n*alpha) - (mw+n)*log(phiw+1) - (mb+n*alpha)*log(phib+zetaT)
      + logP + xib*zetaT - (n*alpha) * log(xib);
  }

  List posterior = List::create(
    Named("omega") = NumericVector::create(mw+n, phiw+1),
    Named("beta") = NumericVector::create(mb+n*alpha, phib+zetaT)
  );

  return List::create(
    Named("omega_shape") = mw+n,
    Named("omega_rate") = phiw+1,
    Named("beta_shape") = mb+n*alpha,
    Named("beta_rate") = phib+zetaT,
    Named("logQ") = logQ,
    Named("zetaT") = zetaT,
    Named("iter") = iter,
    Named("conv") = conv,
    Named("aerror") = aerror,
    Named("rerror") = rerror
  );
}

// @name vbem
// [[Rcpp::export]]

List vb2fit(double alpha, List prior, List data, List options) {

  const int rmax = options["rmax"];
  const double eps = options["eps"];
  const int me = data["total"];

  const double mw = as<NumericVector>(prior["omega"])[0];
  const double mb = as<NumericVector>(prior["beta"])[0];

  double init_zetaT = as<NumericVector>(data["mean"])[0] * as<NumericVector>(data["len"])[0];
  double sumq = 0.0;
  bool conv = false;

  NumericVector nres = NumericVector(0);
  NumericVector qres = NumericVector(0);
  NumericVector wres1 = NumericVector(0);
  NumericVector wres2 = NumericVector(0);
  NumericVector bres1 = NumericVector(0);
  NumericVector bres2 = NumericVector(0);

  for (int r=0; r <= rmax; r++) {
    List result = vbem_vb2gam(me + r, alpha, init_zetaT, prior, data, options);
    double q = exp(as<double>(result["logQ"]));
    if (as<bool>(result["conv"]) == false) {
      warning("VBEM did not converge at r=%d (alpha=%e q=%e)", r, alpha, q);
    }
    sumq += q;
    nres.push_back(me+r);
    qres.push_back(q);
    wres1.push_back(result["omega_shape"]);
    wres2.push_back(result["omega_rate"]);
    bres1.push_back(result["beta_shape"]);
    bres2.push_back(result["beta_rate"]);

//    Rcout << "n=" << me+r << " iter=" << (int) result["iter"] << " convergence=" << (bool) result["conv"]
//     << " q=" << q << " q/sum(q)=" << q/sumq << std::endl;

    if (q / sumq < eps) {
      conv = true;
      break;
    }

    init_zetaT = result["zetaT"];
  }

  List posterior = List::create(
    Named("n") = nres,
    Named("mix") = qres / sumq,
    Named("omega") = List::create(
      Named("shape") = wres1,
      Named("rate") = wres2
    ),
    Named("beta") = List::create(
      Named("shape") = bres1,
      Named("rate") = bres2
    )
  );

  double vfe;
  if (mw == 0 || mb == 0) {
    vfe = 0.0;
  } else {
    vfe = log(sumq);
  }

  return List::create(
    Named("posterior") = posterior,
    Named("vfe") = log(sumq),
    Named("conv") = conv
  );
}

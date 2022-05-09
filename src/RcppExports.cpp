// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// vb1fit
List vb1fit(double alpha, List params, List prior, List data, List options);
RcppExport SEXP _vbsrm_vb1fit(SEXP alphaSEXP, SEXP paramsSEXP, SEXP priorSEXP, SEXP dataSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(vb1fit(alpha, params, prior, data, options));
    return rcpp_result_gen;
END_RCPP
}
// vb2fit
List vb2fit(double alpha, List prior, List data, List options);
RcppExport SEXP _vbsrm_vb2fit(SEXP alphaSEXP, SEXP priorSEXP, SEXP dataSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< List >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(vb2fit(alpha, prior, data, options));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_vbsrm_vb1fit", (DL_FUNC) &_vbsrm_vb1fit, 5},
    {"_vbsrm_vb2fit", (DL_FUNC) &_vbsrm_vb2fit, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_vbsrm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

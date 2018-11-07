// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// r_maternCov
arma::mat r_maternCov(arma::mat dist, double scale, double range, double smoothness, double nugget);
RcppExport SEXP _smolBayes_r_maternCov(SEXP distSEXP, SEXP scaleSEXP, SEXP rangeSEXP, SEXP smoothnessSEXP, SEXP nuggetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type range(rangeSEXP);
    Rcpp::traits::input_parameter< double >::type smoothness(smoothnessSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    rcpp_result_gen = Rcpp::wrap(r_maternCov(dist, scale, range, smoothness, nugget));
    return rcpp_result_gen;
END_RCPP
}
// r_maternArray
arma::vec r_maternArray(arma::vec dist, double scale, double range, double smoothness, double nugget);
RcppExport SEXP _smolBayes_r_maternArray(SEXP distSEXP, SEXP scaleSEXP, SEXP rangeSEXP, SEXP smoothnessSEXP, SEXP nuggetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type range(rangeSEXP);
    Rcpp::traits::input_parameter< double >::type smoothness(smoothnessSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    rcpp_result_gen = Rcpp::wrap(r_maternArray(dist, scale, range, smoothness, nugget));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP t_sfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP t_spredict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_smolBayes_r_maternCov", (DL_FUNC) &_smolBayes_r_maternCov, 5},
    {"_smolBayes_r_maternArray", (DL_FUNC) &_smolBayes_r_maternArray, 5},
    {"t_sfit",     (DL_FUNC) &t_sfit,     14},
    {"t_spredict", (DL_FUNC) &t_spredict,  7},
    {NULL, NULL, 0}
};

RcppExport void R_init_smolBayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

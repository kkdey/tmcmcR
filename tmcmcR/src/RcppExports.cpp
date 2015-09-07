// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/tmcmcR.h"
#include <Rcpp.h>

using namespace Rcpp;

// tmcmcUpdate
NumericVector tmcmcUpdate(NumericVector x, NumericVector b, double eps, F<double> f);
RcppExport SEXP tmcmcR_tmcmcUpdate(SEXP xSEXP, SEXP bSEXP, SEXP epsSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< F<double> >::type f(fSEXP);
    __result = Rcpp::wrap(tmcmcUpdate(x, b, eps, f));
    return __result;
END_RCPP
}
// rwmhUpdate
NumericVector rwmhUpdate(NumericVector x, NumericVector eps, F<double> f);
RcppExport SEXP tmcmcR_rwmhUpdate(SEXP xSEXP, SEXP epsSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< F<double> >::type f(fSEXP);
    __result = Rcpp::wrap(rwmhUpdate(x, eps, f));
    return __result;
END_RCPP
}

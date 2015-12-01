#include <Rcpp.h>
#include "tmcmcR.h"
using namespace Rcpp;
// [[Rcpp::export]]
List tmcmcUpdate(NumericVector x, NumericVector b, double eps, F<double> f){

//  int n = x.size();
  NumericVector y = x + (b*eps);
  double p1= f(x);
  double p2 = f(y);
  double temp = exp(p2-p1);
  double accept = std::min(1.0,temp);
  double u = R::runif(0,1);
  if (u < accept) {
    return List::create(
      _["chain"] = y,
      _["acc_rate"] = accept
    );
    }
  else{
    return List::create(
      _["chain"] = x,
      _["acc_rate"] = accept
    );
    }
}


// [[Rcpp::export]]
List rwmhUpdate(NumericVector x, NumericVector eps, F<double> f){

//  int n = x.size();
  NumericVector y = x + eps;
  double p1= f(x);
  double p2 = f(y);
  double temp = exp(p2-p1);
  double accept = std::min(1.0,temp);
  double u = R::runif(0,1);
  if (u < accept) {
    return List::create(
      _["chain"] = y,
      _["acc_rate"] = accept
    );
    }else{
      return List::create(
        _["chain"] = x,
        _["acc_rate"] = accept
      );
      }
}

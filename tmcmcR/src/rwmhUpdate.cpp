#include <Rcpp.h>
using namespace Rcpp;

template <typename T>
class F {
public:
  F( SEXP f_) : f(f_){}

  inline T operator()(NumericVector x){
    return as<T>(f(x)) ;
  }

private:
  Function f ;
};

// [[Rcpp::export]]

NumericVector rwmhUpdate(NumericVector x, NumericVector eps, F<double> f){

  int n = x.size();
  NumericVector y = x + eps;
  double p1= f(x);
  double p2 = f(y);
  double temp = exp(p2-p1);
  double accept = std::min(1.0,temp);
  double u = R::runif(0,1);
  if (u < accept) {
    return (y);
  }else{
    return (x);
  }
}


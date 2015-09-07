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

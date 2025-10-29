#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double myFunction(double x) {
  return x * x;
}

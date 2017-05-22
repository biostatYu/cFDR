#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector rects(NumericVector x, NumericVector y) {
  int n=x.length();
  LogicalVector ret(n);
  double y0=y[0];
  double x0=x[0];
  for(int i=0; i<n; i++) {
    if(y[i] < y0) {
      ret[i]=false;
      continue;
    }
    ret[i]=true;
    y0=y[i];
  }
  return ret;
}

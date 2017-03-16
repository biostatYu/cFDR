#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector becdf(NumericVector x, NumericVector y) {
  int n=x.length();
  NumericVector cdf(n);
  for(int i=0; i<n; i++) {
    cdf[i]=0;
    for(int j=0; j<n; j++)
      if(x[j]<=x[i] && y[j]<=y[i])
	cdf[i]++;
    cdf[i]=cdf[i]/n;
  }  
  return cdf;
}


// [[Rcpp::export]]
NumericVector cecdf(NumericVector x, NumericVector y) {
  int n=x.length();
  NumericVector cdf(n);
  for(int i=0; i<n; i++) {
    cdf[i]=0;
    int ni=0;
    for(int j=0; j<n; j++) {
      if(y[j]<=y[i]) {
	ni++;
	if(x[j]<=x[i])
	  cdf[i]++;
      }
    }
    cdf[i]=cdf[i]/ni;
  }  
  return cdf;
}

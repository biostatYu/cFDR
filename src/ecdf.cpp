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
// [[Rcpp::export]]
NumericVector cecdf2(NumericVector x1, NumericVector x2, NumericVector y1, NumericVector y2) {
  int n=x1.length();
  NumericVector cdf(n);
  for(int i=0; i<n; i++) {
    cdf[i]=0;
    int ni=0;
    for(int j=0; j<n; j++) {
      if(y1[j]<=y1[i] & y2[j]<=y2[i]) {
	ni++;
	if(x1[j]<=x1[i] & x2[j]<=x2[i])
	  cdf[i]++;
      }
    }
    cdf[i]=cdf[i]/ni;
  }  
  return cdf;
}

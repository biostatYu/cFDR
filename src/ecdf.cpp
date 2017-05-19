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


//' calculate p(X<x|Y<y) for fixed vectors of points x and y
// [[Rcpp::export]]
NumericVector cecdf(NumericVector x, NumericVector y) {
  int n=x.length();
  NumericVector cdf(n);
  for(int i=0; i<n; i++) {
    int mi=0,ni=0;
    for(int j=0; j<n; j++) {
      if(y[j]<=y[i]) {
	ni++;
	if(x[j]<=x[i])
	  mi++;
      }
    }
    cdf[i]=(float)mi/ni;
  }  
  return cdf;
}

//' find weighted mean of values of x either side of target
// [[Rcpp::export]]
NumericVector bestfit(NumericVector x,NumericVector y,NumericVector target) {
  int n=x.length();  
  int hi=0,lo=0;
  double hival=10,loval=-10;
  NumericVector ret(1);
  for(int i=0; i<n; i++) {
    double yt=y[i]-target[0];
    if(yt>0) { // y above target
      if(yt<hival) {
	hival=yt;
	hi=i;
	Rcpp::Rcout << lo << ' ' << loval << ' ' << hi << ' ' << hival << std::endl;
      }
    } else { // y below target
      if(yt>loval) {
	loval=yt;
	lo=i;
	Rcpp::Rcout << lo << ' ' << loval << ' ' << hi << ' ' << hival << std::endl;
      }
    }
  }
  double d= - loval/(hival-loval);
  Rcpp::Rcout << x[lo] << ' ' << loval << ' ' << ' ' << hi << ' ' << hival << ' ' << d << std::endl;
  //  Rccp::Rcout << x[lo] << ' ' << x[hi] << std::endl;
  ret[0] = x[lo] + d * x[hi];
  return ret;
}

// [[Rcpp::export]]
NumericVector fdrp(NumericVector x) {
  int n=x.length();  
  NumericVector y(n);
  for(int i=0; i<n; i++) {
    int ni=0;
    for(int j=0; j<n; j++) {
	if(x[j]<x[i])
	  ni++;
    }
    y[i] = (double) ni; // exp( log(x[i]) + log( (double) n) - log( (double)ni) );
  }  
  return y;
}

// NumericVector findx(NumericVector x, IntegerVector y, IntegerVector ny, NumericVector target) {
//   int n=x.length();
//   NumericVector targx(n);
//   for(int j=0; j<ny[0]; j++) {
    
//     for(int i=0; i<n; i++) {
      
  
// }

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

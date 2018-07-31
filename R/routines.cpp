
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector CaviarLoop(NumericVector Beta, NumericVector y, double EmpiricalQuantile)
{
  int p = y.size();
  int i;
  NumericVector VaR(p);
  VaR[0] = EmpiricalQuantile;
  for(i = 1; i < p; i++)
  {
    /* CAViaR */
    VaR[i] = Beta[0] + Beta[1] * VaR[i-1] + Beta[2] * (y[i-1]);
  }
  return(VaR);
}


// [[Rcpp::export]]
NumericVector iGarchLoop(NumericVector Beta, NumericVector y, double EmpiricalQuantile)
{
  int p = y.size();
  int i;
  NumericVector VaR(p);
  VaR[0] = EmpiricalQuantile;
  for(i = 1; i < p; i++)
  {
    /* Indirect GARCH */
    VaR[i] =  sqrt(Beta[0] + Beta[1] * pow(VaR[i-1],2) + Beta[2] * pow(y[i-1],2));
  }
  return(VaR);
}

// [[Rcpp::export]]
NumericVector AdaptLoop(NumericVector Beta, NumericVector y, double THETA, double G, double EmpiricalQuantile)
{
  int p = y.size();
  int i;
  NumericVector VaR(p);
  VaR[0] = EmpiricalQuantile;
  for(i = 1; i < p; i++)
  {
    /* Adaptive */
    VaR[i] = VaR[i-1] + Beta[0] * (1/(1 + exp(G*(y[i-1]+ VaR[i-1]))) - THETA);		
  }
  return(VaR);
}

// [[Rcpp::export]]
NumericVector AsymSlopeLoop(NumericVector Beta, NumericVector y,  double EmpiricalQuantile)
{
  int p = y.size();
  int i;
  NumericVector VaR(p);
  VaR[0] = EmpiricalQuantile;
  for(i = 1; i < p; i++)
  {
    /* Asymmetric Slope */
    VaR[i] = Beta[0] + Beta[1] * VaR[i-1] + Beta[2] * y[i-1] * (y[i-1] > 0) - Beta[3] * y[i-1] * (y[i-1] < 0); 
  }
  return(VaR);
}

// [[Rcpp::plugins(cpp11)]]
NumericVector mmult(NumericVector a, NumericMatrix b, int p){
  int j; 
  NumericVector c;
  NumericVector out(p);
  for (j = 0; j < p; ++j) {
    c = b(j,_);
    out(j) = std::inner_product(a.begin(), a.end(),c.begin(),0.);              
  }
  return(out);
}

// [[Rcpp::export]]
List MvMqCavLoop(NumericVector c,NumericMatrix a,NumericMatrix r,NumericMatrix y,NumericMatrix X,NumericVector Theta,NumericVector EmpiricalQuantile) {
  int i,j;
  int p = y.ncol(), n = y.nrow();
  double rq = 0, temp = 0;
  NumericMatrix q(n,p);
  
  q.row(0) = EmpiricalQuantile;
  /*1. X - information set up to t-1, computing regression quantiles */
  for (i = 1; i < n; ++i) {
    q(i,_) = c+mmult(q(i-1,_),a,p)+mmult(X(i-1,_),r,p);
  }
  /*2. computing rq statistic for MVMQCAViaR */
  for (i = 0; i < n; ++i) {
    for (j = 0; j < p; ++j){
      /* summing rq statistic through cross section, allowing different theta value for predictors */
      temp += (y(i,j)-q(i,j))*(Theta(j)-(y(i,j)<q(i,j)));
    }
    /* summing rq through time */
    rq += temp;
    temp = 0;
  }
  /* computing the average over time */
  rq = rq/n;
  return  List::create(Rcpp::Named("q") = q,Rcpp::Named("rq") = rq); /*,Rcpp::Named("s") = s,Rcpp::Named("sn") = sn*/
}



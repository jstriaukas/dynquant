#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
List CAVSAVloop(NumericVector BETA, NumericVector y, double empiricalQuantile)
{
  int p          = y.size();
  NumericVector VaR(p);
  int i;
	
	
	/* Initialize output variables */
	VaR[0] = empiricalQuantile;

	/* Start the loop */
	for(i = 1; i < p; i++)
		{
         /* Symmetric Absolute Value */
         VaR[i] = BETA[0] + BETA[1] * VaR[i-1] + BETA[2] * (y[i-1]*(y[i-1]>0) - y[i-1]*(y[i-1]<0));

         }
	return  List::create(Rcpp::Named("VaR") = VaR);
}

// [[Rcpp::export]]
List CAVGARCHloop(NumericVector BETA, NumericVector y, double empiricalQuantile)
{
  int p          = y.size();
  NumericVector VaR(p);
  int i;
  
  
  /* Initialize output variables */
  VaR[0] = empiricalQuantile;
  
  /* Start the loop */
  for(i = 1; i < p; i++)
  {
    /* Indirect GARCH */
    VaR[i] =  sqrt(BETA[0] + BETA[1] * pow(VaR[i-1],2) + BETA[2] * pow(y[i-1],2));
    
  }
  return  List::create(Rcpp::Named("VaR") = VaR);
}


// [[Rcpp::export]]
List CAVloop(NumericVector BETA, NumericVector y, double empiricalQuantile)
{
  int p          = y.size();
  NumericVector VaR(p);
  int i;
  
  
  /* Initialize output variables */
  VaR[0] = empiricalQuantile;
  
  /* Start the loop */
  for(i = 1; i < p; i++)
  {
     
    VaR[i] = BETA[0] + BETA[1] * VaR[i-1] + BETA[2] * (y[i-1]);
    
  }
  return  List::create(Rcpp::Named("VaR") = VaR);
}

// [[Rcpp::export]]
List ADAPTIVEloop(NumericVector BETA, NumericVector y, double THETA, double K, double empiricalQuantile)
{
  int p          = y.size();
  NumericVector VaR(p);
  int i;
  
  
  /* Initialize output variables */
  VaR[0] = empiricalQuantile;
  
  /* Start the loop */
  for(i = 1; i < p; i++)
  {
    /* Adaptive */
    VaR[i] = VaR[i-1] + BETA[0] * (1/(1 + exp(K*(y[i-1]+ VaR[i-1]))) - THETA);		
    
  }
  return  List::create(Rcpp::Named("VaR") = VaR);
}



// [[Rcpp::export]]
List ASYMloop(NumericVector BETA, NumericVector y,  double empiricalQuantile)
{
  int p          = y.size();
  NumericVector VaR(p);
  int i;
  
  
  /* Initialize output variables */
  VaR[0] = empiricalQuantile;
  
  /* Start the loop */
  for(i = 1; i < p; i++)
  {
    /* Asymmetric Slope */
    VaR[i] = BETA[0] + BETA[1] * VaR[i-1] + BETA[2] * y[i-1] * (y[i-1] > 0) - BETA[3] * y[i-1] * (y[i-1] < 0); 
  }
  return  List::create(Rcpp::Named("VaR") = VaR);
}


// [[Rcpp::export]]
List CAVMIDASloop(NumericVector BETA, NumericVector y, NumericMatrix x, NumericVector seq, double nlag, double empiricalQuantile)
{
  int p          = y.size();
  NumericVector VaR(p);
  NumericVector weights(nlag);
  NumericVector W(nlag);
  NumericVector X(p);
  int i, j, k, d;
  
  for (j = 0; j<nlag; j++) 
  {
     weights[j] = pow(1-seq[j], BETA[3]-1);
  }
  double sumW = sum(weights);
    
  
  for (j = 0; j<nlag; j++) 
  {
    W[j] = weights[j]/sumW;
  }
  double D = 0;
  for (k = 0 ; k<p ; k++){  /* Loop through observations */
    NumericVector M = x(k,_);
    for ( d = 0; d<p ; d++){  /* Loop through HF lags */
       D +=(M[d]*W[d]);
    }
    X[k] = D;
  }
  
  
  /* Initialize output variables */
  VaR[0] = empiricalQuantile;
  for(i = 1; i < p; i++)
  {
    /* Symmetric Absolute Value */ 
    VaR[i] = BETA[0] + BETA[1] * VaR[i-1] + BETA[2] * (X[i-1]*(X[i-1]>0) - X[i-1]*(X[i-1]<0));
    
  }
  return  List::create(Rcpp::Named("VaR") = VaR);
}




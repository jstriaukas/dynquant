
rm(list = ls())
setwd('/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/R')
#require('zoo')
require('parallel')
require('optimx')
require("cmaes")
#require('ggplot2')
source('functions.R')

set.seed(11)
load("/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/Data/dataQuantileMIDAS.RData")
type <- "cav"
quant.type <- "quant"
is.midas <- TRUE
opt.method <- "cma-es"
opt.transform <- "abs"
mc <- TRUE
us <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)
z <- list(y = diff(log(DataQuantile$SNP))*100)
dq.options <- NULL
for(i in 1:length(us)) {
  theta <- us[i]
  cat("quantile = ", format(theta), "..")
  fit <- fit.u.dyn.quant(theta,z,type,is.midas,opt.method,mc=mc,quant.type=quant.type,period=252,nlag=252)
  if (i==1){
    plot(fit$Y, type = 'l', lty=2, col="black")
    quantMat <- array(0,c(length(fit$fitted.values),length(us)))
    ceoffMat <- array(0,c(length(coef(fit)),length(us)))
  }
  lines(fit$fitted.values, lty=2, col="red")
  quantMat[,i] <- fit$fitted.values
  ceoffMat[,i] <- coef(fit)
  cat("\n")
}  
  

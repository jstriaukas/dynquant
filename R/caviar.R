rm(list = ls())
setwd('/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/R')
#require('zoo')
require('parallel')
require('optimx')
require("cmaes")
#require('ggplot2')
source('functions.R')

set.seed(10)

load("/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/Data/dataCAViaR.RData")
GM <- dataCAViaR$V1
IBM <- dataCAViaR$V2
SnP500 <- dataCAViaR$V3


z <- NULL
inSample <- 2892
theta      <- 0.01
# GM
z$y                        <- GM[1:inSample]
temp                       <- sort(z$y[1:300])
empiricalQuantile          <- temp[300*z$theta]
type <- "cav"
quant.type <- "var"
is.midas <- TRUE
opt.method <- "cma-es"
mc <- TRUE
us <- c(0.05,0.25,0.5,0.75,0.95)
z <- list(y = GM[1:inSample])
theta 
cat("quantile = ", format(theta), "..")
fit <- fit.dyn.quant(theta,z,type,is.midas,opt.method,mc=mc,quant.type=quant.type,period=252,nlag=252)
  
est.GM                     <- fit.dyn.quant(z, type, mc = TRUE,  quant.type = "var", empiricalQuantile  = -empiricalQuantile)
# IBM
z$y                        <- IBM[1:inSample]
temp                       <- sort(z$y[1:300])
empiricalQuantile          <- temp[300*z$theta]
est.IBM                    <- fit.dyn.quant(z, type, mc = TRUE,  quant.type = "var",    empiricalQuantile  = -empiricalQuantile)
# SNP500
z$y                        <- SnP500[1:inSample]
temp                       <- sort(z$y[1:300])
empiricalQuantile          <- temp[300*z$theta]
est.SNP500                 <- fit.dyn.quant(z, type, mc = TRUE,  quant.type = "var",    empiricalQuantile  = -empiricalQuantile)




rm(list = ls())
setwd('/Users/striaukas/Documents/GitHub/dynquant')

#require('zoo')
require('parallel')
require('optimx')
require("cmaes")
#change for testing
setwd('/Users/striaukas/Documents/GitHub/dynquant')
#require('ggplot2')
source('R/functions.R')
set.seed(11)
load("Data/dataQuantileMIDAS.RData")
#load("C:/Users/Aq/Documents/R/win-library/3.5/bindrcpp/dataQuantileMIDAS.RData")
type <- "cav"
quant.type <- "quant"
is.midas <- TRUE
opt.method <- "cma-es"
opt.transform <- "abs"
mc <- FALSE
thetas <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)
z <- list(y = diff(log(DataQuantile$SNP))*100, date = DataQuantile$DATE[2:length(DataQuantile$DATE)])
fit <- fit.u.dyn.quants(thetas,z,type,is.midas,opt.method,opt.transform,mc,quant.type,isplot = TRUE, title = "SNP Quantiles")















rm(list = ls())
setwd('/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/R')
#require('zoo')
require('parallel')
require('optimx')
require("cmaes")
#require('ggplot2')
source('functions.R')
load("/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/Data/dataCAViaR.RData")
type <- "cav"
quant.type <- "quant"
is.midas <- FALSE
opt.method <- "cma-es"
opt.transform <- c("abs","sq")
mc <- TRUE
theta <- c(0.01,0.5,0.5)
z <- list(y = as.matrix(dataCAViaR))
fit <- fit.mv.dyn.quant(theta,z,type,is.midas,opt.method,opt.transform,mc,quant.type)

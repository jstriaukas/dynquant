rm(list = ls())
#change for testing
setwd('/Users/striaukas/Documents/GitHub/dynquant')
#require('zoo')
require('parallel')
require('optimx')
require("cmaes")
#require('ggplot2')
source('R/functions.R')
load("Data/dataCAViaR.RData")
type <- "cav"
quant.type <- "quant"
is.midas <- TRUE
opt.method <- "cma-es"
opt.transform <- c("abs","abs")
mc <- TRUE
theta <- c(0.05,0.05)
z <- list(y = as.matrix(dataCAViaR[,2:3]))
dq.options <- NULL
min.evals <- 2
fit <- fit.mv.dyn.quant(theta,z,type,is.midas,opt.method,opt.transform,mc,quant.type,min.evals)

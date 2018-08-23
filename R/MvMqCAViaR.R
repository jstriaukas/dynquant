rm(list = ls())
#change for testing
setwd('/Users/striaukas/Documents/GitHub/dynquant')
#require('zoo')
require('parallel')
require('optimx')
require("cmaes")
#require('ggplot2')
source('R/functions.R')
load("Data/MvMQCAViaRBanks.RData")
y <- rowMeans(data)
r <- log(y[2:length(y)]/y[1:length(y)-1])
x.1 <- data[,116] # wells fargo
x.2 <- data[,65] # jp morgan

r.1 <- log(x.1[2:length(x.1)]/x.1[1:length(x.1)-1])
r.2 <- log(x.2[2:length(x.2)]/x.2[1:length(x.2)-1])
y.date <- date[2:length(date)]
req.lf <- '1-m'
dat.r <- mixed.freq.data.ref.lf(req.lf,r,y.date,22,1,y.date[44],y.date[length(y.date)-22])
dat.x.1 <- mixed.freq.data.ref.lf(req.lf,r.1,y.date,22,1,y.date[44],y.date[length(y.date)-22], FALSE)
dat.x.2 <- mixed.freq.data.ref.lf(req.lf,r.2,y.date,22,1,y.date[44],y.date[length(y.date)-22], FALSE)
X.1 <- X.2 <- NULL
X.1[[1]] <- X.2[[1]] <-  dat.r$est.x
X.1[[2]] <- dat.x.1$est.x
X.2[[2]] <- dat.x.2$est.x

z.1 <- list(y = cbind(dat.r$est.y, dat.x.1$est.y), x = X.1)
z.2 <- list(y = cbind(dat.r$est.y, dat.x.2$est.y), x = X.2)

type <- "cav"
quant.type <- "quant"
is.midas <- TRUE
opt.method <- "cma-es"
opt.transform <- c("abs","abs")
mc <- TRUE
theta <- c(0.05,0.05)
dq.options <- NULL
min.evals <- 5
fit.1 <- fit.mv.dyn.quant(theta,z.1,type,is.midas,opt.method,opt.transform,mc,quant.type,min.evals)
fit.2 <- fit.mv.dyn.quant(theta,z.2,type,is.midas,opt.method,opt.transform,mc,quant.type,min.evals)

plot(fit.1$VaR[,1], type = 'l')
plot(fit.1$VaR[,2], type = 'l')
fit.1$k2

plot(fit.2$VaR[,1], type = 'l')
plot(fit.2$VaR[,2], type = 'l')
fit.2$k2


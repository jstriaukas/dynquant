rm(list = ls())
set.seed(123)
setwd('/Users/striaukas/Documents/GitHub/dynquant')
source('R/data_functions.R')
#data
require("R.matlab")
pathname <- file.path("/Users/striaukas/Dropbox/PhD/Projects/stock-bond-mv-quantile/data/clean/", "data_sw_1961_2017.mat")
data <- readMat(pathname)
daily <- data$data[[1]]
monthly <- data$data[[2]]
x.1 <- daily[,,1]$sp.500.daily.ret[2:length(daily[,,1]$sp.500.daily.ret)]
x.2 <- daily[,,1]$yields[,13]
x.2 <- diff(x.2)
x.date <- daily[,,1]$date
x.date <- as.Date(x.date, origin = "1970-01-01")-719529
x.date <- x.date[2:length(x.date)]
y.1 <- monthly[,,1]$sp.500.monthly.ret
y.2 <- monthly[,,1]$monthly.ret[,5]
y.date <- monthly[,,1]$date
y.date <- as.Date(y.date, origin = "1970-01-01")-719529
### MIDAS data ###
x.lag <- 22
horizon <- 1
est.start <- y.date[2]
est.end <- y.date[length(y.date)]
dat.x.1 <- mixed.freq.data.single(y.date,x.1,x.date,x.lag,horizon,est.start,est.end,disp.flag=TRUE)
dat.x.2 <- mixed.freq.data.single(y.date,x.2,x.date,x.lag,horizon,est.start,est.end,disp.flag=TRUE)

#change for testing
setwd('/Users/striaukas/Documents/GitHub/dynquant')

#require('zoo')
require('parallel')
require('optimx')
require("cmaes")
#require('ggplot2')
source('R/functions.R')

X.1  <- NULL


ind <- c(-length(y.1)+1,-length(y.1))

X.1[[1]] <- dat.x.1$est.x[ind[1],]
X.1[[2]] <- dat.x.2$est.x[ind[1],]
z.1 <- list(y = cbind(y.1[ind], y.2[ind]), x = X.1, date = dat.x.1$est.refdate[2:length(dat.x.1$est.refdate)])

type <- "cav"
quant.type <- "quant"
is.midas <- TRUE
opt.method <- "cma-es"
opt.transform <- c("abs","abs")
mc <- TRUE
thetas <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)
min.evals <- 10
wd <- "/Users/striaukas/Dropbox/PhD/Projects/stock-bond-mv-quantile/results/figures"
fit.1 <- fit.mv.dyn.quants(thetas,z.1,type,is.midas,opt.method,opt.transform,mc,quant.type,min.evals = min.evals,isplot = FALSE,titles = c("SNP Quantiles", "10 Bond Quantiles"),pdftitle = "mv-cav-stock-bond",wd = wd)




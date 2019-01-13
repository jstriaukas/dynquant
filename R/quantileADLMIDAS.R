rm(list = ls())
setwd('/Users/striaukas/Documents/GitHub/dynquant')
### DATA ###
rm(list = ls())
setwd('/Users/striaukas/Documents/GitHub/dynquant/R/')
source('data_functions.R')

est.y <- read.csv("~/Documents/GitHub/dynquant/R/EstY.txt", sep="")
est.ydate.imp <- read.table("~/Documents/GitHub/dynquant/R/EstYdate.csv", quote="\"", comment.char="")
est.x <- read.csv("~/Documents/GitHub/dynquant/R/EstX.txt", sep="")
est.xdate.imp <- read.table("~/Documents/GitHub/dynquant/R/EstXdate.csv", quote="\"", comment.char="")

est.ydate <- est.xdate <- NULL

for (i in 1:dim(est.ydate.imp)[1]) {
  est.ydate[i] <- as.Date(est.ydate.imp[i,], origin = "1899-12-30") 
}
est.ydate <- as.Date(est.ydate, origin = "1970-01-01")
for (j in 1:dim(est.xdate.imp)[1]) {
  est.xdate[j] <- as.Date(est.xdate.imp[j,], origin = "1899-12-30") 
}
est.xdate <- as.Date(est.xdate, origin = "1970-01-01")


data.y <- est.y 
data.ydate <- est.ydate
data.x <- est.x
data.xdate <- est.xdate

est.start <-  as.Date('1964-02-01')
est.end   <- as.Date('2017-09-01')

x.lag <- 66
y.lag <- 1
horizon <- '1d'

data <- mixed.freq.data(data.y,data.ydate,data.x,data.xdate,x.lag,y.lag,horizon,est.start,est.end)





### MODELING ###
require('parallel')
require('optimx')
require("cmaes")
#change for testing
setwd('/Users/striaukas/Documents/GitHub/dynquant')
#require('ggplot2')
source('R/functions.R')
set.seed(11)
type <- "ar.rq"
quant.type <- "quant"
is.midas <- TRUE
opt.method <- "cma-es"
opt.transform <- "lev"
mc <- TRUE
thetas <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)
z <- list(y = data$est.y, x.highfreq = data$est.x, date = data$est.ydate)
fit <- fit.u.dyn.quants(thetas,z,type,is.midas,opt.method,opt.transform,mc,quant.type,isplot = TRUE, title = "SMTH",y.lowfreqlag = data$est.lag.y)
#fit <- fit.u.dyn.quant(thetas[1],z,type,is.midas,opt.method,opt.transform,mc,quant.type,isplot = TRUE, title = "SMTH", y.lowfreqlag = data$est.lag.y)










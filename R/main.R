rm(list = ls(all = T))
# get data #

setwd('C:/Users/Z440/Documents/FINMETRICS/Data')

load('dataqmidas.RData')

setwd('C:/Users/Z440/Documents/FINMETRICS/R')
source('functions.R')

betacoeff <- NULL
betacoeff[1] <- 0.005
betacoeff[2] <- 12
betacoeff[3] <- 15
y       <- dataqmidas$US
.computeQMIDAS(betacoeff, y, nlag = 251, period = 251, theta = 0.05)

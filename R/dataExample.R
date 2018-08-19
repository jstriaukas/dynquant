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

data <- mixed.freq.data.single(data.ydate,data.x,data.xdate,x.lag,horizon,est.start,est.end)

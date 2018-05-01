setwd('/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/R')
require('parallel')
require('neldermead')
require('cmaes')
require('optimx')
source('functions.R')

set.seed(50)
iscaviar <- 1
load("/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/Data/dataCAViaR.RData")
GM  <- dataCAViaR$V1[1:2892]
IBM <- dataCAViaR$V2
SNP <- dataCAViaR$V3

plot(GM,  type = 'l')
plot(IBM, type = 'l')
plot(SNP, type = 'l')

z <- NULL

##### 1 % symmetric absolute value #####
type         <- 'symabs'
z$theta      <- 0.01
# GM
z$y                 <- GM
EST.GM.01.symabs    <- fitdynquant(z, type, iscaviar, mc = TRUE)
# IBM
z$y                 <- IBM
EST.IBM.01.symabs   <- fitdynquant(z, type, iscaviar, mc = TRUE)
# SnP500
z$y                 <- SNP
EST.SNP.01.symabs   <- fitdynquant(z, type, iscaviar, mc = TRUE)
##### 1 % asymmetric slope #####
type         <- 'asymslope'
# GM
z$y                    <- GM
EST.GM.01.asymslope    <- fitdynquant(z, type, iscaviar, mc = TRUE)
# IBM
z$y                    <- IBM
EST.IBM.01.asymslope   <- fitdynquant(z, type, iscaviar, mc = TRUE)
# SnP500
z$y                    <- SNP
EST.SNP.01.asymslope   <- fitdynquant(z, type, iscaviar, mc = TRUE)
##### 1 % integrated garch #####
type         <- 'igarch'
# GM
z$y                    <- GM
EST.GM.01.igarch    <- fitdynquant(z, type, iscaviar, mc = TRUE)
# IBM
z$y                    <- IBM
EST.IBM.01.igarch   <- fitdynquant(z, type, iscaviar, mc = TRUE)
# SnP500
z$y                    <- SNP
EST.SNP.01.igarch   <- fitdynquant(z, type, iscaviar, mc = TRUE)
##### 1 % adaptive #####
type         <- 'adapt'
# GM
z$y                    <- GM
EST.GM.01.adapt    <- fitdynquant(z, type, iscaviar, mc = TRUE)
# IBM
z$y                    <- IBM
EST.IBM.01.adapt   <- fitdynquant(z, type, iscaviar, mc = TRUE)
# SnP500
z$y                    <- SNP
EST.SNP.01.adapt   <- fitdynquant(z, type, iscaviar, mc = TRUE)






##### 5 % symmetric absolute value #####
type         <- 'symabs'
z$theta      <- 0.05
# GM
z$y                 <- GM
EST.GM.05.symabs    <- fitdynquant(z, type, iscaviar, mc = TRUE)
# IBM
z$y                 <- IBM
EST.IBM.05.symabs   <- fitdynquant(z, type, iscaviar, mc = TRUE)
# SnP500
z$y                 <- SNP
EST.SNP.05.symabs   <- fitdynquant(z, type, iscaviar, mc = TRUE)
##### 5 % asymmetric slope #####
type         <- 'asymslope'
# GM
z$y                    <- GM
EST.GM.05.asymslope    <- fitdynquant(z, type, iscaviar, mc = TRUE)
# IBM
z$y                    <- IBM
EST.IBM.05.asymslope   <- fitdynquant(z, type, iscaviar, mc = TRUE)
# SnP500
z$y                    <- SNP
EST.SNP.05.asymslope   <- fitdynquant(z, type, iscaviar, mc = TRUE)
##### 5 % integrated garch #####
type         <- 'igarch'
# GM
z$y                    <- GM
EST.GM.05.igarch    <- fitdynquant(z, type, iscaviar, mc = TRUE)
# IBM
z$y                    <- IBM
EST.IBM.05.igarch   <- fitdynquant(z, type, iscaviar, mc = TRUE)
# SnP500
z$y                    <- SNP
EST.SNP.05.igarch   <- fitdynquant(z, type, iscaviar, mc = TRUE)
##### 5 % adaptive #####
type         <- 'adapt'
# GM
z$y                    <- GM
EST.GM.05.adapt    <- fitdynquant(z, type, iscaviar, mc = TRUE)
# IBM
z$y                    <- IBM
EST.IBM.05.adapt   <- fitdynquant(z, type, iscaviar, mc = TRUE)
# SnP500
z$y                    <- SNP
EST.SNP.05.adapt   <- fitdynquant(z, type, iscaviar, mc = TRUE)

plot(EST.GM.05.symabs$VaR, type = 'l')
plot(EST.GM.05.asymslope$VaR, type = 'l')
plot(EST.GM.05.igarch$VaR, type = 'l')
plot(EST.SNP.05.adapt$VaR, type = 'l')
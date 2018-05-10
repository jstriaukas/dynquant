setwd('/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/R')
require('parallel')
require('optimx')
require('ggplot2')
source('functions.R')

set.seed(10)

load("/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/Data/dataCAViaR.RData")
GM <- dataCAViaR$V1
IBM <- dataCAViaR$V2
SnP500 <- dataCAViaR$V3


z <- NULL
inSample <- 2892
##### 1 % symmetric absolute value #####
type         <- 'symabs'
z$theta      <- 0.01
# GM
z$y                        <- GM[1:inSample]
temp                       <- sort(z$y[1:300])
empiricalQuantile          <- temp[300*z$theta]
fit                     <- fit.dyn.quant(z, type, mc = TRUE,  quant.type = "var", empiricalQuantile  = -empiricalQuantile)
# IBM
z$y                        <- IBM[1:inSample]
temp                       <- sort(z$y[1:300])
empiricalQuantile          <- temp[300*z$theta]
est.IBM                    <- fit.dyn.quant(z, type, mc = TRUE,  quant.type = "var",    empiricalQuantile  = -empiricalQuantile)
# SNP500
z$y                        <- SNP500[1:inSample]
temp                       <- sort(z$y[1:300])
empiricalQuantile          <- temp[300*z$theta]
est.SNP500                 <- fit.dyn.quant(z, type, mc = TRUE,  quant.type = "var",    empiricalQuantile  = -empiricalQuantile)




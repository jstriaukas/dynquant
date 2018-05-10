setwd('/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/R')
require('parallel')
require('optimx')
require('ggplot2')
source('functions.R')

set.seed(10)

load("/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/Data/dataQuantilMIDAS.RData")
snp.ret <- log(DataQuantile$US[2:length(DataQuantile$US)]/DataQuantile$US[1:length(DataQuantile$US)-1])

type         <- 'qmidas'
z            <- NULL
z$theta      <- 0.01
# GM
z$y                     <- snp.ret
temp                    <- sort(z$y[1:300])
empiricalQuantile       <- temp[300*z$theta]
fit                     <- fit.dyn.quant(z, type, mc = TRUE,  quant.type = "var", nlag = 252, period = 252)

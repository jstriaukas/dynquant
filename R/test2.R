load("/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/Data/dataqmidas.RData")
require('parallel')
source('functions.R')

idx <- dataqmidas$US
ret <- log(idx[seq(2,length(idx))]/idx[seq(1,length(idx)-1)])*100
plot(ret, type = 'l')

z <- NULL
z$theta <- 0.05
z$y <- ret

iscaviar <- 0
type <- 'symabs'
EST <- fitdynquant(z, type, iscaviar, mc = TRUE)

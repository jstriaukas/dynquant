

betacoeff <- NULL
betacoeff[1] <- 0.02
betacoeff[2] <- 0.6
betacoeff[3] <- 0.9
type <- 'simple'
EST <-.computeCAViaR(betacoeff, z, type)
plot(EST$VaR, type ='l')

betacoeff <- NULL
betacoeff[1] <- 0.02
betacoeff[2] <- 0.6
betacoeff[3] <- 0.9
type <- 'symabs'
EST <-.computeCAViaR(betacoeff, z, type)
plot(EST$VaR, type ='l')

betacoeff <- NULL
betacoeff[1] <- 0.02
betacoeff[2] <- 0.6
betacoeff[3] <- 0.9
type <- 'igarch'
EST <-.computeCAViaR(betacoeff, z, type)
plot(EST$VaR, type ='l')


betacoeff <- NULL
betacoeff[1] <- 0.02
betacoeff[2] <- 0.6
betacoeff[3] <- 0.9
betacoeff[4] <- 0.9
type <- 'asymslope'
EST <-.computeCAViaR(betacoeff, z, type)
plot(EST$VaR, type ='l')



betacoeff <- NULL
betacoeff[1] <- 0.02
betacoeff[2] <- 0.6
betacoeff[3] <- 0.9
type <- 'hybrid'
EST <-.computeCAViaR(betacoeff, z, type)
plot(EST$VaR, type ='l')


betacoeff <- NULL
betacoeff[1] <- 0.02
type <- 'adapt'
EST <-.computeCAViaR(betacoeff, z, type, K = 200)
plot(EST$VaR, type ='l')


midasdata <- getMIDASstucture(ret, nlag = 22, period = 22, theta = z$theta)

betacoeff <- NULL
betacoeff[1] <- 0.02
betacoeff[2] <- 0.6
betacoeff[3] <- 0.9
EST <- .computeQMIDAS(betacoeff, midasdata)



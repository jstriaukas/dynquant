#setwd('C:/Users/Z440/Documents/FINMETRICS/R')
setwd('/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/R')
Rcpp::sourceCpp('routines.cpp')

load("/Users/striaukas/Dropbox/PhD/Projects/FINMETRICS/Data/dataqmidas.RData")

fitdynquant <- function(z, type, iscaviar, mc = FALSE, ...) {
  
  Z <- list(...)
  
  if(length(Z$num.bestevals)==0){Z$num.bestevals = 40}
  
  CONST         <- getconstraints(type, iscaviar)
  starting.vals <- getinitial(type, iscaviar)
  RQ            <- NULL
  
  if (iscaviar == 0){
    z <- getMIDASstucture(z$y, period = Z$period, nlag = Z$nlag)
  }
  
  if(!mc)RQ <-lapply(1:dim(starting.vals)[1], getbestintial, 
                     type, z, iscaviar, starting.vals) 
  
  if(mc){ if(length(Z$ncores) == 0) { ncores = detectCores()} }
  if(mc)RQ <- mclapply(1:dim(starting.vals)[1], getbestintial, 
                       type, z, iscaviar, starting.vals, mc.cores = ncores)
  
  RQ       <- unlist(RQ)
  all.vals <- cbind(RQ, starting.vals)
  all.vals <- all.vals[order(RQ),]
  
  best.evals <- all.vals[1:Z$num.bestevals,2:dim(all.vals)[2]]

  if(!mc) est  <-  lapply(1:dim(best.evals)[1], getoptimal, 
                          best.evals, CONST, z, type, iscaviar)
  if(mc)  est  <-  mclapply(1:dim(best.evals)[1], getoptimal, 
                            best.evals, CONST, z, type, iscaviar, mc.cores = ncores)
  
  val <- NULL
  for(j in 1:length(est)){val[j] <- est[[j]]$value }
  
  temp         <- est[min(val)==val]
  z$est.params <- as.numeric(temp[[1]][1:dim(best.evals)[2]])
  
  
  if(iscaviar == 1){
    z <- computeCAViaR(z$est.params, z, type, out.val = 1)
  } else {
    z <- computeQMIDAS(z$est.params, z, out.val = 1)
  }
 
  return(z = z)
}

computeQMIDAS <- function(betacoeff, z, out.val) {
  
  
  Y         <- z$yLowFreq
  X         <- z$xHighFreq
  THETA     <- z$theta
  
  intercept <- betacoeff[1]
  slope     <- betacoeff[2]
  k2        <- betacoeff[3]
  nlag      <- ncol(z$xHighFreq)
  k1        <- 1
  
  weights   <- betacoeffWeights(nlag, k1, k2)
  W         <- weights$weights
  
  CondQuant <- intercept + slope * (z$xHighFreq%*%W)
  
  VaR <- -CondQuant
  Hit <- (Y < -VaR) - THETA
  RQ  <-  -t(Hit)%*%(Y + VaR)
  
  if (out.val == 1){
    z$VaR <- VaR
    z$Hit <- Hit
    z$RQ  <- RQ
    z$CondQuant <- CondQuant
    class(z) <- 'dymquant'
  } else {
    z <- RQ
  }
  
  return(z)
}

computeCAViaR <- function(betacoeff, z, type=c('symabs','igarch','simple','adapt','asymslope', 'hybrid'), out.val = 0) {
  
  # Preliminaries
  Y                    <- z$y
  Tobs                 <- length(Y)
  ysort                <- sort(Y) 
  THETA                <- z$theta
  empiricalQuantile    <- ysort[round(Tobs*THETA)]
  type                 <- match.arg(type)
  
  if(type == 'symabs') {
    if(length(betacoeff) != 3){stop("betacoeff is of the size 3x1")}
    VaR   <- CAVSAVloop(betacoeff, Y, -empiricalQuantile)
  } else if(type == 'igarch'){
    if(length(betacoeff) != 3){stop("betacoeff is of the size 3x1")}
    VaR   <- CAVGARCHloop(betacoeff, Y, -empiricalQuantile)
  } else if(type == 'simple'){
    if(length(betacoeff) != 3){stop("betacoeff is of the size 3x1")}
    VaR   <- CAVloop(betacoeff, Y, -empiricalQuantile)
  } else if(type == 'adapt'){
    if(length(betacoeff) != 2){stop("betacoeff is of the size 2x1")}
    G     <- betacoeff[2]
    VaR   <- ADAPTIVEloop(betacoeff[1], Y, THETA, G, -empiricalQuantile)
  } else if(type == 'asymslope'){ 
    if(length(betacoeff) != 4){stop("betacoeff is of the size 4x1")}
    VaR   <- ASYMloop(betacoeff, Y, -empiricalQuantile)
  } else {
    stop('CAViaR type either not implemented or does not exist')
  }
  VaR <- VaR$VaR
  Hit <- (Y < -VaR) - THETA
  RQ  <-  -t(Hit)%*%(Y + VaR)
  
  if (out.val == 1){
    z$VaR <- VaR
    z$Hit <- Hit
    z$RQ <- RQ
    z$CondQuant <- -VaR
    class(z) <- 'dymquant'
  } else if(out.val == 0){
    z <- as.numeric(RQ)
  }
  
  
  return(z)
  
}

getMIDASstucture <- function(y, ...){ 
  
  Z         <- list(...)
  NLAG      <- Z$nlag
  PERIOD    <- Z$period
  if (length(PERIOD) == 0) { PERIOD = 22 }
  if (length(NLAG) == 0) { NLAG = 22 }
  
  Tobs      <- length(y)
  
  nobsShort <- Tobs-NLAG-PERIOD+1
  data <- getmidasdata(y, nobsShort, NLAG, PERIOD)
  z$yLowFreq  <- data$yLowFreq
  z$xHighFreq <- data$xHighFreq
  
  return(z)
  
}

getmidasdata     <- function(y, nobsShort, nlag = 22, period = 22){
  
  Tobs      <- length(y)
  yLowFreq  <- array(0, c(nobsShort,1))
  xHighFreq <- array(0, c(nobsShort,nlag))
  Regressor <- abs(y)
  
  for (t in (nlag+1):(Tobs-period+1)) {
    yLowFreq[t-nlag,1]  = sum(y[t:t+period-1])  
    xHighFreq[t-nlag, ] = t(Regressor[seq(t-1, t-nlag, by = -1)])
  }
  
  return(list(yLowFreq = yLowFreq, xHighFreq = xHighFreq))
}

betacoeffWeights <- function(nLag, coeff1, coeff2) {
  temp <- seq(0,1,length=nLag)
  if (coeff1==1){
    weights <- (1-temp)^(coeff2-1)
  } else {
    weights <- (1-temp)^(coeff2-1) * (temp^(coeff1-1))
  }
  weights <- weights/sum(weights)
  
  return(list(weights = weights))
}

getconstraints   <- function(type, iscaviar){
  
  const  <- NULL
  
  if (iscaviar == 1) {
    
    if(type == 'symabs') {

      const$LB[1] <- -100  
      const$LB[2] <- -1.5  
      const$LB[3] <- -100
      const$UB[1] <-  100  
      const$UB[2] <-  1.5  
      const$UB[3] <-  100
      
    } else if(type == 'igarch'){
     
      const$LB[1] <- -100  
      const$LB[2] <-  0  
      const$LB[3] <-  0
      const$UB[1] <-  100  
      const$UB[2] <- 100  
      const$UB[3] <- 100
      
    } else if(type == 'simple'){

      const$LB[1] <- -100  
      const$LB[2] <- -1.5  
      const$LB[3] <- -100
      const$UB[1] <-  100  
      const$UB[2] <-  1.5  
      const$UB[3] <-  100
      
    } else if(type == 'adapt'){

      const$LB[1] <- -100  
      const$LB[2] <- 0  
      const$UB[1] <-  100  
      const$UB[2] <- 100000        
            
    } else if(type == 'asymslope'){ 

      const$LB[1] <- -100  
      const$LB[2] <- -1.5  
      const$LB[3] <- -100  
      const$LB[4] <- -100
      const$UB[1] <-  100  
      const$UB[2] <-  1.5  
      const$UB[3] <-  100  
      const$UB[4] <-  100     
      
    } else {
      stop('CAViaR type either not implemented or does not exist')
    }
    
    
    
  } else {
    
    const$LB[1] <- -100  
    const$LB[2] <- -100  
    const$LB[3] <-    1
    const$UB[1] <-  100  
    const$UB[2] <-  100  
    const$UB[3] <-  100
    
  }
  
  return(const)
  
}

getinitial       <- function(type, iscaviar, num.test = 1e5){
  
  CONST         <- getconstraints(type, iscaviar)
  
  
  if (iscaviar == 1) {
    if(type == 'symabs') {
      starting.vals     <- array(NA, c(num.test,3))
      starting.vals[,1] <- runif(num.test, min =  0, 1) 
      starting.vals[,2] <- runif(num.test, min =  0, 1) 
      starting.vals[,3] <- runif(num.test, min =  0, 1)
      
    } else if(type == 'igarch'){
      starting.vals     <- array(NA, c(num.test,3))
      starting.vals[,1] <- runif(num.test, min =  0, 1) 
      starting.vals[,2] <- runif(num.test, min =  0, 1) 
      starting.vals[,3] <- runif(num.test, min =  0, 1)
      
    } else if(type == 'simple'){
      starting.vals     <- array(NA, c(num.test,3))
      starting.vals[,1] <- runif(num.test, min =  0, 1) 
      starting.vals[,2] <- runif(num.test, min =  0, 1) 
      starting.vals[,3] <- runif(num.test, min =  0, 1)

      
    } else if(type == 'adapt'){
      starting.vals     <- array(NA, c(num.test,2))
      starting.vals[,1] <- runif(num.test, min =  0,     1) 
      starting.vals[,2] <- runif(num.test, min =  0,  1000) 
      
    } else if(type == 'asymslope'){ 
      starting.vals     <- array(NA, c(num.test,4))
      starting.vals[,1] <- runif(num.test, min =  0, 1) 
      starting.vals[,2] <- runif(num.test, min =  0, 1) 
      starting.vals[,3] <- runif(num.test, min =  0, 1)   
      starting.vals[,4] <- runif(num.test, min =  0, 1)  
    } else {
      stop('CAViaR type either not implemented or does not exist')
    }
    
  } else {
    
    starting.vals     <- array(NA, c(num.test,3))
    starting.vals[,1] <- runif(num.test, min =  -1,   1)
    starting.vals[,2] <- runif(num.test, min =  -1,   1)
    starting.vals[,3] <- runif(num.test, min =    1, 30)   
    
    
  }
  
  return(starting.vals)
  
}


getbestintial    <- function(dt, type, z, iscaviar, starting.vals) {

  if (iscaviar ==  1) {
      temp  <- computeCAViaR(starting.vals[dt,], z, type, out.val = 1)
      RQ <- temp$RQ
  } else {
      temp  <- computeQMIDAS(starting.vals[dt,], z, out.val = 1)
      RQ <- temp$RQ
  }
  return(RQ = RQ)
}

getoptimal       <- function(dt, param,  CONST, z, type, iscaviar) {
  
REP <- 10
if (iscaviar == 1) {
        est <-  optimx(param[dt,], computeCAViaR,  z=z, type=type, out.val = 0, method = c("Nelder-Mead")) 
      for (i in 1:REP) {
        est <-  optimx(as.numeric(est[1:length(param[dt,])]), computeCAViaR,  z=z, type=type, out.val = 0, method = c("Nelder-Mead"))
        est <-  optimx(as.numeric(est[1:length(param[dt,])]), computeCAViaR,  z=z, type=type, out.val = 0, method = c("BFGS"))
      }
        est <-  optimx(param[dt,], computeCAViaR,  z=z, type=type, out.val = 0, method = c("Nelder-Mead"))   
  } else {
      est <- optimx(param[dt, ], computeQMIDAS, z=z, out.val = 0, out.val = 0, method = c("Nelder-Mead")) 
      for (i in 1:REP) {
        est <-  cma_es(est$pars, computeQMIDAS, z=z, out.val = 0, out.val = 0, method = c("Nelder-Mead")) 
      }      
  }
      return(est = est)
}

VarianceCovariance <- function(z, type, ...) {
  Z   <- list(...)
  
  lags <- Z$lags
  if (length(lags) == 0) { lags = 4 }
  
  BETA  <- z$est.params
  THETA <- z$theta
  y     <- z$y
  Tobs  <- length(y)
  VaR   <- z$VaR
  
  # Compute the quantile residuals.
  residuals <- y + VaR
  
  # Set up the bandwidth for the k-Nearest neighbor estimator.
  SortedRes <- sort(abs(residuals), decreasing =FALSE)
  
  if (length(residuals) < 40) {
      k <- 30
  }  else if (THETA == 0.01) {
      k <- 40
  }  else if (THETA == 0.025) {
      k <- 50
  }  else {
      k <- 60
  }

  BANDWIDTH <- SortedRes[k]
  
  # Initialize matrices.
  derivative1 <- array(0, c(Tobs,1))
  derivative2 <- array(0, c(Tobs,1))
  derivative3 <- array(0, c(Tobs,1))
  derivative4 <- array(0, c(Tobs,1))
  
  D <- array(0, c(length(BETA),length(BETA)))
  A <- D
  t <- 0 # counter 
  ##################################
  # Model: Simple lagged value.
  if (type == 'simple') {
    gradient <- array(0, c(Tobs, 3))
    
    for (i in 2:Tobs){
      derivative1[i] <- 1 + BETA[2] * derivative1[i-1]
      derivative2[i] <- VaR[i-1] + BETA[2] * derivative2[i-1]
      derivative3[i] <- BETA[2] * derivative3[i-1] + y[i-1]
      
      gradient[i,]  <- c(derivative1[i], derivative2[i], derivative3[i])
      
      A <- A + gradient[i,]%*%t(gradient[i,])
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gradient[i,]%*%t(gradient[i,])
      }
      
    }
  }
  # Model: Symmetric Absolute Value.
  if (type == 'symabs') {
    gradient <- array(0, c(Tobs, 3))
    
    for (i in 2:Tobs){
      derivative1[i] <- 1 + BETA[2] * derivative1[i-1]
      derivative2[i] <- VaR[i-1] + BETA[2] * derivative2[i-1]
      derivative3[i] <- BETA[2] * derivative3[i-1] + abs(y[i-1])
      
      gradient[i,]  <- c(derivative1[i], derivative2[i], derivative3[i])
      
      A <- A + gradient[i,]%*%t(gradient[i,])
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gradient[i,]%*%t(gradient[i,])
      }
      
    }
  }
  
  # Model: Asymmetric Slope.
  if (type == 'asymslope') {
    gradient <- array(0, c(Tobs, 4))
    
    for (i in 2:Tobs){
      derivative1[i] <- 1 + BETA[2] * derivative1[i-1]
      derivative2[i] <- VaR[i-1] + BETA[2] * derivative2[i-1]
      derivative3[i] <- BETA[2]*derivative3[i-1] + y[i-1]*(y[i-1]>0)
      derivative4[i] <- BETA[2]*derivative4[i-1] - y[i-1]*(y[i-1]<0)
      
      gradient[i,]  <- c(derivative1[i], derivative2[i], derivative3[i], derivative4[i])
      
      A <- A + gradient[i,]%*%t(gradient[i,])
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gradient[i,]%*%t(gradient[i,])
      }
      
    }
  }
  
  # Model: iGARCH.
  if (type == 'igarch') {
    gradient <- array(0, c(Tobs, 3))
    
    for (i in 2:Tobs){
      derivative1[i] <- (1 + 2*BETA[2]*VaR[i-1]*derivative1[i-1]) / (2*VaR[i])
      derivative2[i] = (VaR[i-1]^2 + 2*BETA[2]*VaR[i-1]*derivative2[i-1]) / (2*VaR[i]);
      derivative3[i] = (2*BETA[2]*VaR[i-1]*derivative3[i] + y[i-1]^2) / (2*VaR[i])
      
      gradient[i,]  <- c(derivative1[i], derivative2[i], derivative3[i])
      
      A <- A + gradient[i,]%*%t(gradient[i,])
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gradient[i,]%*%t(gradient[i,])
      }
      
    }
  }
  
  # Model: Adaptive.
  if (type == 'adapt') {
    gradient <- array(0, c(Tobs, 2))
    
    for (i in 2:Tobs){
      indicator <-  exp(BETA[2] * (y[i-1] + VaR[i-1]))
      if (indicator == Inf){
        derivative1[i] <- derivative1[i-1]
      } else {
        derivative1[i] <- derivative1[i-1] + 1/(1 + indicator) - THETA + BETA[1] * BETA[2] * indicator * derivative1[i-1] / (1 + indicator)^2
        derivative2[i] <- BETA[1] * indicator * (y[i-1] + VaR[i-1]) / (1 + indicator)^2
      }
      
      gradient[i,]  <- c(derivative1[i], derivative2[i])
      
      A <- A + gradient[i,]%*%t(gradient[i,])
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gradient[i,]%*%t(gradient[i,])
      }
      
    }
  }
  
  z$tStdErr  <- t # Check the k-NN bandwidth
  A          <- A/Tobs
  z$D        <- D  / (2*BANDWIDTH*Tobs)
  z$gradient <- gradient
  z$VCmatrix <- THETA * (1-THETA) * solve(z$D) * A * solve(z$D) / Tobs
  
  return(z = z)
}


output       <- DQtest(OUT, MODEL, T, z, THETA, VaR, Hit, D, gradient)



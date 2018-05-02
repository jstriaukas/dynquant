Rcpp::sourceCpp('routines.cpp')

fit.dyn.quant <- function(z, type, iscaviar, mc = FALSE, quant.type = c("var", "quant"), ...) {
  call <- match.call()
  Z <- list(...)
  
  if(mc){ if(length(Z$ncores) == 0) { ncores = detectCores()} }
  
  if(length(Z$num.bestevals)==0){Z$num.bestevals = 10}
  
  if(length(Z$empiricalQuantile)==0){
    if(iscaviar==1){
      temp <- order(z$y)
      Z$empiricalQuantile <- temp[100*z$theta]
    }
      Z$empiricalQuantile <- NULL
  }
  empiricalQuantile <- Z$empiricalQuantile

  quant.type           <- match.arg(quant.type)
  if (quant.type != "var" && quant.type != "quant"){
    quant.type = "var"
    warning("Either quant.type not set or this option is not available. In this case, VaR is defined as a positive")
  }
  
  constraints         <- get.constraints(type, iscaviar)
  starting.vals <- get.initial(type, iscaviar)
  RQ            <- NULL
  
  if (iscaviar == 0){
    if (length(Z$yLowFreq) ==0 || length(Z$xHighFreq)){
    z <- get.midas.structure(z$y, period = Z$period, nlag = Z$nlag)
    } else {
    z$xHighFreq <- Z$xHighFreq  
    z$yLowFreq  <- Z$yLowFreq
  }
  } 
  
  if(!mc)RQ <-lapply(1:dim(starting.vals)[1], get.best.intial, 
                     type, z, iscaviar, starting.vals, quant.type, empiricalQuantile) 
  

  if(mc)RQ <- mclapply(1:dim(starting.vals)[1], get.best.intial, 
                       type, z, iscaviar, starting.vals, quant.type, empiricalQuantile, mc.cores = ncores)
  
  
  
  RQ       <- unlist(RQ)
  all.vals <- cbind(RQ, starting.vals)
  all.vals <- all.vals[order(RQ),]
  
  best.evals <- all.vals[1:Z$num.bestevals,2:dim(all.vals)[2]]
  
  if(!mc) est  <-  lapply(1:dim(best.evals)[1], get.optimal, 
                          best.evals, CONST, z, type, iscaviar, quant.type, empiricalQuantile)
  if(mc)  est  <-  mclapply(1:dim(best.evals)[1], get.optimal, 
                            best.evals, CONST, z, type, iscaviar, quant.type, empiricalQuantile, mc.cores = ncores)
  
  val <- NULL
  for(j in 1:length(est)){val[j] <- est[[j]]$value }
  
  temp         <- est[min(val)==val]
  fit$est.params <- as.numeric(temp[[1]][1:dim(best.evals)[2]])
  
  
  if(iscaviar == 1){
    fit <- compute.caviar(fit$est.params, z, type, out.val = 1,  quant.type, empiricalQuantile)
  } else {
    fit <- compute.qmidas(fit$est.params, z, out.val = 1, quant.type)
  }
  
  fit <- dq.stat(fit, type, quant.type, lags = 4)
  fit$std.error <- sqrt(diag(fit$VCmatrix))
  fit$pval      <- pnorm(-abs(fit$est.params)/fit$std.error)
  fit$call      <- call
  return(fit = fit)
}

compute.qmidas <- function(betacoeff, z, out.val, quant.type) {
  
  
  Y         <- z$yLowFreq
  X         <- z$xHighFreq
  THETA     <- z$theta
  
  quant.type           <- match.arg(quant.type)
  
  intercept <- betacoeff[1]
  slope     <- betacoeff[2]
  k2        <- betacoeff[3]
  nlag      <- ncol(z$xHighFreq)
  k1        <- 1
  
  weights   <- beta.coeff.w(nlag, k1, k2)
  W         <- weights$weights
  
  CondQuant <- intercept + slope * (z$xHighFreq%*%W)
  
  VaR <- -CondQuant
  VaR <- CondQuant
  
  Hit <- (Y < -VaR) - THETA
  RQ  <-  -t(Hit)%*%(Y + VaR)
  
  if (out.val == 1){
    if(quant.type == "var"){
      z$VaR <- VaR
    } else {
      z$VaR <- -VaR
    }
    z$Hit <- Hit
    z$RQ <- RQ
    z$CondQuant <- CondQuant
    class(z) <- "dymquant"
  } else if(out.val == 0){
    z <- as.numeric(RQ)
  }
  
  return(z)
}

compute.caviar <- function(betacoeff, z, type=c('symabs','igarch','simple','adapt','asymslope', 'hybrid'), out.val = 0, quant.type, empiricalQuantile) {
  
  Y                    <- z$y
  THETA                <- z$theta
  type                 <- match.arg(type)
  
  if(type == 'symabs') {
    if(length(betacoeff) != 3){stop("betacoeff is of the size 3x1")}
    VaR   <- CAVSAVloop(betacoeff, Y, empiricalQuantile)
  } else if(type == 'igarch'){
    if(length(betacoeff) != 3){stop("betacoeff is of the size 3x1")}
    VaR   <- CAVGARCHloop(betacoeff, Y, empiricalQuantile)
  } else if(type == 'simple'){
    if(length(betacoeff) != 3){stop("betacoeff is of the size 3x1")}
    VaR   <- CAVloop(betacoeff, Y, empiricalQuantile)
  } else if(type == 'adapt'){
    if(length(betacoeff) != 2){stop("betacoeff is of the size 2x1")}
    G     <- betacoeff[2]
    VaR   <- ADAPTIVEloop(betacoeff[1], Y, THETA, G, empiricalQuantile)
  } else if(type == 'asymslope'){ 
    if(length(betacoeff) != 4){stop("betacoeff is of the size 4x1")}
    VaR   <- ASYMloop(betacoeff, Y, empiricalQuantile)
  } else {
    stop('Please choose a different type of CAViaR specification (see documentation for more details)')
  }
  VaR <-  VaR$VaR
  Hit <- (Y < -VaR) - THETA
  RQ  <-  -t(Hit)%*%(Y + VaR)
  
  if (out.val == 1){
    if(quant.type == "var"){
      z$VaR <- VaR
    } else {
      z$VaR <- -VaR
    }
    z$Hit <- Hit
    z$RQ <- RQ
    z$CondQuant <- -VaR
    z$quant.type <- quant.type
    class(z) <- "dymquant"
    
  } else if(out.val == 0){
    z <- as.numeric(RQ)
  }
  
  
  return(z)
  
}

get.midas.structure <- function(y, ...){ 
  
  Z         <- list(...)
  nlag      <- Z$nlag
  PERIOD    <- Z$period
  if (length(PERIOD) == 0) { PERIOD = 22 }
  if (length(nlag) == 0) { nlag = 22 }
  
  Tobs      <- length(y)
  
  nobsShort <- Tobs-nlag-PERIOD+1
  data <- get.midas.data(y, nobsShort, nlag, PERIOD)
  z$yLowFreq  <- data$yLowFreq
  z$xHighFreq <- data$xHighFreq
  
  return(z)
  
}

get.midas.data     <- function(y, nobsShort, nlag = 22, period = 22){
  
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

beta.coeff.w <- function(nlag, coeff1, coeff2) {
  temp <- seq(0,1,length=nlag)
  if (coeff1==1){
    weights <- (1-temp)^(coeff2-1)
  } else {
    weights <- (1-temp)^(coeff2-1) * (temp^(coeff1-1))
  }
  weights <- weights/sum(weights)
  
  return(list(weights = weights))
}

get.constraints   <- function(type, iscaviar){
  
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

get.initial       <- function(type, iscaviar, num.test = 1e4){
  
  CONST         <- get.constraints(type, iscaviar)
  
  
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
    starting.vals[,3] <- runif(num.test, min =   1,  30)   
    
    
  }
  
  return(starting.vals)
  
}

get.best.intial    <- function(dt, type, z, iscaviar, starting.vals, quant.type, empiricalQuantile) {
  
  if (iscaviar ==  1) {
      temp  <- compute.caviar(starting.vals[dt,], z, type, out.val = 1, quant.type, empiricalQuantile)
      RQ <- temp$RQ
  } else {
      temp  <- compute.qmidas(starting.vals[dt,], z, out.val = 1, quant.type)
      RQ <- temp$RQ
  }
  return(RQ = RQ)
}

get.optimal       <- function(dt, param,  CONST, z, type, iscaviar, quant.type, empiricalQuantile) {
REP <- 10
if (iscaviar == 1) {
        est <-  optimx(param[dt,], compute.caviar,  z=z, type=type, out.val = 0, quant.type=quant.type, empiricalQuantile=empiricalQuantile, method = c("Nelder-Mead")) 
      for (i in 1:REP) {
        est <-  optimx(as.numeric(est[1:length(param[dt,])]), compute.caviar,  z=z, type=type, out.val = 0, quant.type=quant.type, empiricalQuantile=empiricalQuantile, method = c("Nelder-Mead"))
        est <-  optimx(as.numeric(est[1:length(param[dt,])]), compute.caviar,  z=z, type=type, out.val = 0, quant.type=quant.type, empiricalQuantile=empiricalQuantile, method = c("BFGS"))
      }
        est <-  optimx(as.numeric(est[1:length(param[dt,])]), compute.caviar,  z=z, type=type, out.val = 0, quant.type=quant.type, empiricalQuantile=empiricalQuantile, method = c("Nelder-Mead"))
  } else {
        est <- optimx(param[dt, ], compute.qmidas, z=z, out.val = 0, out.val = 0, quant.type=quant.type, method = c("Nelder-Mead")) 
      for (i in 1:REP) {
        est <-  optimx(as.numeric(est[1:length(param[dt,])]), compute.qmidas,  z=z, type=type, out.val = 0, quant.type, method = c("Nelder-Mead"))
        est <-  optimx(as.numeric(est[1:length(param[dt,])]), compute.qmidas,  z=z, type=type, out.val = 0, quant.type, method = c("BFGS"))
      } 
        est <-  optimx(as.numeric(est[1:length(param[dt,])]), compute.qmidas,  z=z, type=type, out.val = 0, quant.type, method = c("Nelder-Mead"))
  }
      return(est = est)
}

var.cov.mat.est <- function(z, type) {
  
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
    gradient <- array(0, c(Tobs,3))
    
    for (i in 2:Tobs){
      derivative1[i] <- 1 + BETA[2] * derivative1[i-1]
      derivative2[i] <- VaR[i-1] + BETA[2] * derivative2[i-1]
      derivative3[i] <- BETA[2] * derivative3[i-1] + y[i-1]
      
      gradient[i,]  <- c(derivative1[i], derivative2[i], derivative3[i])
      gtg           <- gradient[i,]%*%t(gradient[i,])
      A             <- A + gtg     
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gtg
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
      gtg           <- gradient[i,]%*%t(gradient[i,])
      A             <- A + gtg
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gtg
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
      gtg           <- gradient[i,]%*%t(gradient[i,])
      A             <- A + gtg     
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gtg
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
      gtg           <- gradient[i,]%*%t(gradient[i,])
      A             <- A + gtg
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gtg
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
      gtg           <- gradient[i,]%*%t(gradient[i,])
      A             <- A + gtg
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gtg
      }
      
    }
  }
  z$k         <- k
  z$BANDWIDTH <- BANDWIDTH
  z$tStdErr   <- t # Check the k-NN bandwidth
  z$A         <- A/Tobs
  z$D         <- D/(2*BANDWIDTH*Tobs)
  z$gradient  <- gradient
  z$VCmatrix  <- THETA*(1-THETA)*solve(z$D)%*%z$A%*%solve(z$D)/Tobs
  
  return(z = z)
}

dq.stat <- function(z, type, quant.type, ...){
  opt <- list(...)
  if(length(opt$lags) == 0){lags = 4} else {lags = opt$lags}
  
  z         <- var.cov.mat.est(z, type)
  y         <- z$y
  VaR       <- z$VaR
  gradient  <- z$gradient
  Hit       <- z$Hit
  BANDWIDTH <- z$BANDWIDTH
  D         <- z$D
  THETA     <- z$theta
  Tobs      <- length(y)
  if(quant.type == "var"){
    residuals <- y + VaR
  } else {
    residuals <- y - VaR
  }
  constant    <- matrix(1,Tobs-lags,1)
  HIT         <- Hit[seq(lags+1,Tobs, by = 1)]
  VaRforecast <- VaR[seq(lags+1,Tobs, by = 1)]
  yLag        <- y[seq(lags,Tobs-1, by = 1)]
  VaRforecastlag <- VaR[seq(lags,Tobs-1, by = 1)]
  Z           <- array(0, c(Tobs-lags, lags))

  for (s in 1:lags){
    Z[,s]              <- Hit[seq(s,Tobs-(lags+1-s), by = 1)]
  }
  
  Xout <- cbind(constant, VaRforecast, Z) # Instruments for the out of sample test.
  Xin  <- Z                        # Instruments for the in sample test.
  
  XHNABLA   <- array(0,c(dim(Xin)[2],length(z$est.params)))
  NABLA     <- gradient[seq(lags+1,Tobs, by=1),]
  
  #-- Estimate the matrices that enter the In Sample DQ test --#
  
  for(j in 2:(Tobs-lags)){
    if(abs(residuals[j]) <= BANDWIDTH){
      XHNABLA <- XHNABLA + Xin[j,]%*%t(gradient[j,])
    }
  }
  
  XHNABLA <- XHNABLA/(2*BANDWIDTH*Tobs)
  M         <- t(Xin) - XHNABLA%*%solve(D)%*%t(NABLA)
  
  #-- Compute the DQ tests --#
  z$dq.statIn  <- (t(HIT)%*%Xin%*%solve(M%*%t(M))%*%t(Xin)%*%HIT) / (THETA*(1-THETA))
  z$dq.statOut <- (t(HIT)%*%Xout%*%solve(t(Xout)%*%Xout)%*%t(Xout)%*%HIT) / (THETA*(1-THETA))
  
  z$DQin      <- 1 - pchisq(z$dq.statIn, dim(Xin)[2]) # Compute the P-value of the in sample DQ test.
  z$DQout     <- 1 - pchisq(z$dq.statOut, dim(Xout)[2]) # Compute the P-value of the out of sample DQ test.

  return(z)
}



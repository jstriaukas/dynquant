Rcpp::sourceCpp('R/routines.cpp')
source('R/data_functions.R')
source('R/plot_functions.R')

get.mv.asymmetry <- function(thetas,z,type=c("cav"),is.midas,opt.method=c("nelder-mead","cma-es"),opt.transform,mc=FALSE,quant.type=c("var","quant"),...){
  call <- match.call()
  type <- match.arg(type)
  opt.method <- match.arg(opt.method)
  quant.type <- match.arg(quant.type)
  dq.options <- list(...)
  if(is.null(dq.options$isplot)){dq.options$isplot = TRUE}
  fit <- NULL
  p <- ncol(z$y)
  for(i in 1:length(thetas)){
    theta <- thetas[i]
    cat("quantile = ", format(theta), "..")
    fit.mv <- fit.mv.dyn.quant(theta,z,type,is.midas,opt.method,opt.transform,mc,quant.type,dq.options)
    if (i==1){
      fit$quant.mat <- array(0,c(dim(fit.mv$fitted.values)[1],dim(fit.mv$fitted.values)[2],length(thetas)))
      fit$ceoff.mat <- array(0,c(length(coef(fit.mv)),dim(fit.mv$fitted.values)[2],length(thetas)))
    }
    fit$quant.mat[,,i] <- fit.mv$fitted.values
    fit$all[[i]] <- fit.mv
    cat("\n")
  }
  if(dq.options$isplot) {
    for (i in 1:p){
      plot.asymmetry(fit,z,
                     title = dq.options$titles[i], 
                     xlabelaxis = dq.options$xlabelaxis,
                     ylabelaxis = dq.options$ylabelaxis,
                     pdftitle = paste(dq.options$pdftitle,i,".pdf", sep="-"),
                     wd = dq.options$wd)
    }
  }
  return(fit)
  
}


fit.mv.dyn.quants <- function(thetas,z,type=c("cav"),is.midas,opt.method=c("nelder-mead","cma-es"),opt.transform,mc=FALSE,quant.type=c("var","quant"),...){
  call <- match.call()
  type <- match.arg(type)
  opt.method <- match.arg(opt.method)
  quant.type <- match.arg(quant.type)
  dq.options <- list(...)
  if(is.null(dq.options$isplot)){dq.options$isplot = TRUE}
  fit <- NULL
  p <- ncol(z$y)
  for(i in 1:length(thetas)){
    theta <- thetas[i]
    cat("quantile = ", format(theta), "..")
    fit.mv <- fit.mv.dyn.quant(theta,z,type,is.midas,opt.method,opt.transform,mc,quant.type,dq.options)
    if (i==1){
      fit$quant.mat <- array(0,c(dim(fit.mv$fitted.values)[1],dim(fit.mv$fitted.values)[2],length(thetas)))
    }
    fit$quant.mat[,,i] <- fit.mv$fitted.values
    fit$all[[i]] <- fit.mv
    cat("\n")
  }
  if(dq.options$isplot) {
    for (i in 1:p){
    fit.single <- NULL 
    fit.single$quant.mat <- fit$quant.mat[,i,]
    z.single <- NULL
    z.single$y <- z$y[,i]
    z.single$date <- z$date
    plot.quantiles(fit.single,z.single,
                   title = dq.options$titles[i], 
                   xlabelaxis = dq.options$xlabelaxis,
                   ylabelaxis = dq.options$ylabelaxis,
                   pdftitle = paste(dq.options$pdftitle,i,".pdf", sep="-"),
                   wd = dq.options$wd)
    }
  }
  return(fit)
  
}


fit.mv.dyn.quant <- function(theta,z,type=c("cav"),is.midas,opt.method=c("nelder-mead","cma-es"),opt.transform,mc=FALSE,quant.type=c("var","quant"),...) {
  call <- match.call()
  type <- match.arg(type)
  opt.method <- match.arg(opt.method)
  quant.type <- match.arg(quant.type)
  dq.options <- list(...)
  if(mc){ 
    if(is.null(dq.options$ncores)) {ncores = detectCores()
    } else {ncores = dq.options$ncores}
  }
  p <- ncol(z$y)
  if(is.null(dq.options$min.evals)){dq.options$min.evals = 1}
  if(is.null(dq.options$num.test)){dq.options$num.test <- 1e3}
  if(length(theta)==1){theta=rep(theta,p)}
  if(length(theta)!=p){stop("number of quantile levels must equal to the number of covariates")}
  if(length(is.midas)==1){is.midas=rep(is.midas,p)}
  if(length(is.midas)>p){stop("number of MIDAS weighted covariates must not exceed total number of covariates")}
  if(length(opt.transform)<p){opt.transform=c(opt.transform,rep("lev",p-length(opt.transform)))}
  if(length(opt.transform)>p){stop("number of transformations exceed total number of covariates")}
  min.evals <- dq.options$min.evals
  Y<-X<-a<-r<-c<-fit.u<-dat<-e.q<-kappa<-NULL
  for(i in 1:p){
    if (is.midas[i]==TRUE){
      dat$y.lowfreq <- z$y[,i]
      dat$x.highfreq <- z$x[[i]]
    } else{
      dat$y <- z$y[,i]
    }
    fit.u[[i]] <- fit.u.dyn.quant(theta[i],dat,type,is.midas[i],opt.method,opt.transform[i],mc,quant.type,dq.options)
    coeffs <- coef(fit.u[[i]])
    c[i] <- coeffs[1]
    a[i] <- coeffs[2] # autoregressive coefficient
    r[i] <- coeffs[3] # realized coefficient 
    kappa[i] <- coeffs[4]
    e.q[i] <- fit.u[[i]]$empirical.quantile
    Y[[i]] <- fit.u[[i]]$Y
    X[[i]] <- fit.u[[i]]$X
  }
  if(is.null(dq.options$empirical.quantile)){dq.options$empirical.quantile <- e.q}
  empirical.quantile <- dq.options$empirical.quantile
  
  pars0 <- c(c,c(diag(a)),c(diag(r)),c(kappa))
  Y <- matrix(unlist(Y),ncol=p,byrow=FALSE)
  z$Y <- Y
  z$X <- X
  starting.vals <- get.mv.initial(pars0,type,is.midas,p,num.test = 1e3)
  if(!mc) rq.stat <-lapply(1:dim(starting.vals)[1],get.mv.best.init, 
                           theta,z,is.midas,starting.vals,
                           quant.type,empirical.quantile) 
  if(mc) rq.stat <- mclapply(1:dim(starting.vals)[1],get.mv.best.init, 
                             theta,z,is.midas,starting.vals,
                             quant.type,empirical.quantile,
                             mc.cores = ncores)
  rq.stat <- unlist(rq.stat)
  all.evals <- cbind(rq.stat,starting.vals)
  all.evals <- all.evals[order(rq.stat),]
  min.rq.stat.evals <- array(all.evals[1:min.evals,2:dim(all.evals)[2]],c(min.evals,dim(all.evals)[2]-1))
  min.rq.stat.evals <- min.rq.stat.evals[,!is.na(pars0)]
  min.rq.stat.evals <- array(matrix(min.rq.stat.evals),c(min.evals,dim(all.evals)[2]-1))
  constraints <- get.mv.constraints("cav",is.midas,p)
  constraints$LB <- constraints$LB[!is.na(pars0)]
  constraints$UB <- constraints$UB[!is.na(pars0)]
  if(!mc) est  <-  lapply(1:min.evals,get.mv.optimal, 
                          min.rq.stat.evals,theta,p,z,is.midas,quant.type,empirical.quantile,opt.method,constraints)
  if(mc)  est  <-  mclapply(1:min.evals,get.mv.optimal, 
                            min.rq.stat.evals,theta,p,z,is.midas,quant.type,empirical.quantile,opt.method,constraints,mc.cores=ncores)
  val <- NULL
  for(j in 1:length(est)){val[j] <- est[[j]]$value }
  temp <- est[min(val)==val]
  if(opt.method=="cma-es"){coeff <- temp[[1]]$par}
  if(opt.method=="nelder-mead"){coeff <- as.numeric(temp[[1]][1:dim(min.rq.stat.evals)[2]])}
  fit <- compute.mv.quantile(coeff,theta,z,is.midas,out.val=1,quant.type,empirical.quantile)
  return(fit)
}

compute.mv.quantile <- function(pars0,theta,z,is.midas,out.val=0,quant.type="var",empirical.quantile) {
  Y <- z$Y
  X <- z$X
  p <- ncol(Y)
  n <- nrow(Y)
  Xt <- array(0,c(n,p))
  k2 <- NULL
  for(i in 1:p) { 
    if(is.midas[i]){
      k1 <- 1
      k2[i] <- pars0[length(pars0)]
      pars0 <- pars0[1:length(pars0)-1]
      nlag <- ncol(X[[i]])
      weights <- beta.poli.w(nlag, k1, k2)
      W <- weights$weights
      Xt[,i] <- (X[[i]]%*%W)
    } else{
      Xt[,i] <- X[[i]]
    }
  }
  C <- pars0[1:p]
  A <- matrix(pars0[(p+1):(p^2+p)],nrow=p,ncol=p)
  R <- matrix(pars0[(p^2+p+1):(2*p^2+p)],nrow=p,ncol=p)
  res <- MvMqCavLoop(C,A,R,Y,Xt,theta,empirical.quantile)
  hit.stat <- theta - (Y < res$q)
  rq.stat <- res$rq
  if (out.val == 1){
    if(quant.type == "var"){
      z$VaR <- -res$q
    } else {
      z$VaR <- res$q
    }
    z$hit.stat <- hit.stat
    z$rq.stat <- rq.stat
    z$residuals <- Y - res$q
    z$fitted.values <- res$q
    z$c <- C
    z$a <- A
    z$r <- R
    z$k2 <- k2
    class(z) <- "dynquant"
  } else if(out.val == 0){
    z <- as.numeric(rq.stat)
  }
  return(z)
  
}

fit.u.dyn.quants <- function(thetas,z,type=c("cav","igarch","adapt","asymslope","rq"),is.midas,opt.method=c("nelder-mead","cma-es"),opt.transform=c("lev","abs","sq"),mc=FALSE,quant.type=c("var","quant"),...){
call <- match.call()
type <- match.arg(type)
opt.method <- match.arg(opt.method)
opt.transform <- match.arg(opt.transform)
quant.type <- match.arg(quant.type)
dq.options <- list(...)
if(is.null(dq.options$isplot)){dq.options$isplot = TRUE}
fit <- NULL
for(i in 1:length(thetas)){
  theta <- thetas[i]
  cat("quantile = ", format(theta), "..")
  fit.u <- fit.u.dyn.quant(theta,z,type,is.midas,opt.method,opt.transform,mc,quant.type,dq.options)
  if (i==1){
    fit$quant.mat <- array(0,c(length(fit.u$fitted.values),length(thetas)))
  }
  fit$quant.mat[,i] <- fit.u$fitted.values
  fit$all[[i]] <- fit.u
  cat("\n")
}
  if(dq.options$isplot) {
    plot.quantiles(fit,z,
                   title = dq.options$title, 
                   xlabelaxis = dq.options$xlabelaxis,
                   ylabelaxis = dq.options$ylabelaxis,
                   pdftitle = dq.options$pdftitle,
                   wd = dq.options$wd)
  }
  return(fit)

}



fit.u.dyn.quant <- function(theta,z,type=c("cav","igarch","adapt","asymslope","rq"),is.midas,opt.method=c("nelder-mead","cma-es"),opt.transform=c("lev","abs","sq"),mc=FALSE,quant.type=c("var","quant","es"),...) {
  call <- match.call()
  type <- match.arg(type)
  opt.method <- match.arg(opt.method)
  opt.transform <- match.arg(opt.transform)
  quant.type <- match.arg(quant.type)
  dq.options <- list(...)
  if(mc){ 
    if(is.null(dq.options$ncores)) {
      ncores = detectCores()
    }
    else {
      ncores = dq.options$ncores
    }
  }
  if(is.null(dq.options$min.evals)){dq.options$min.evals = 5}
  min.evals <- dq.options$min.evals
  if(is.null(dq.options$empirical.quantile)){
      temp <- sort(z$y[1:100],decreasing=FALSE)
    dq.options$empirical.quantile <- temp[100*theta]
  }
  empirical.quantile <- dq.options$empirical.quantile

  if(is.midas){  
    if(opt.transform=="lev"){z$X <- z$x.highfreq}
    if(opt.transform=="abs"){z$X <- abs(z$x.highfreq)}
    if(opt.transform=="sq"){z$X <- (z$x.highfreq)^2}
    z$Y <- z$y.lowfreq
  } else {
    if(opt.transform=="lev"){z$X <- z$y}
    if(opt.transform=="abs"){z$X <- abs(z$y)}
    if(opt.transform=="sq"){z$X <- (z$y)^2}
    z$Y <- z$y
  }
  starting.vals <- get.u.initial(type,is.midas)
  if(!mc) rq.stat <-lapply(1:dim(starting.vals)[1], get.u.best.init, 
                           theta,z,type,is.midas,starting.vals,
                           quant.type,empirical.quantile) 
  if(mc) rq.stat <- mclapply(1:dim(starting.vals)[1], get.u.best.init, 
                             theta,z,type,is.midas,starting.vals,
                             quant.type,empirical.quantile,
                             mc.cores = ncores)
  rq.stat <- unlist(rq.stat)
  all.evals <- cbind(rq.stat,starting.vals)
  all.evals <- all.evals[order(rq.stat),]
  min.rq.stat.evals <- all.evals[1:dq.options$min.evals,2:dim(all.evals)[2]]
  min.rq.stat.evals <- matrix(as.numeric(min.rq.stat.evals), dq.options$min.evals,dim(all.evals)[2]-1)
  constraints <- get.u.constraints(type,is.midas)
  if(!mc) est  <-  lapply(1:min.evals,get.u.optimal, 
                          min.rq.stat.evals,theta,z,type,is.midas,quant.type,empirical.quantile,opt.method,constraints)
  if(mc)  est  <-  mclapply(1:min.evals,get.u.optimal, 
                            min.rq.stat.evals,theta,z,type,is.midas,quant.type,empirical.quantile,opt.method,constraints,mc.cores=ncores)
  val <- NULL
  for(j in 1:length(est)){val[j] <- est[[j]]$value }
  temp <- est[min(val)==val]
  if(opt.method=="cma-es"){coeff <- temp[[1]]$par}
  if(opt.method=="nelder-mead"){coeff <- as.numeric(temp[[1]][1:dim(min.rq.stat.evals)[2]])}
  fit <- compute.u.quantile(coeff,theta,z,type,is.midas,out.val=1,quant.type,empirical.quantile)
  fit$coefficients <- coeff
  fit$empirical.quantile <- empirical.quantile
  #fit$constraints <- constraints
  #fit <- dq.stat(fit, type, quant.type, lags = 4)
  #fit$std.error <- sqrt(diag(fit$VCmatrix))
  #fit$pval <- pnorm(-abs(fit$est.params)/fit$std.error)
  fit$call <- call
  return(fit = fit)
}

compute.u.quantile <- function(pars0,theta,z,type=c("symabs","igarch","cav","adapt","asymslope","rq"),is.midas,out.val=0,quant.type="var",empirical.quantile) {
  if(is.midas){
    k1 <- 1
    k2 <- pars0[length(pars0)]
    pars0 <- pars0[1:length(pars0)-1]
    nlag <- ncol(z$x.highfreq)
    weights <- beta.poli.w(nlag, k1, k2)
    W <- weights$weights
    X <- (z$X%*%W)
    Y <- z$Y
  } else{
    X <- z$X
    Y <- z$Y
  }
  type <- match.arg(type)
  if(type=="cav"){
    if(length(pars0)!=3){stop("beta is of the size 3x1")}
    q <- CaviarLoop(pars0,X,empirical.quantile)
  } else if(type=="igarch"){
    if(length(pars0)!=3){stop("beta is of the size 3x1")}
    q <- iGarchLoop(pars0,X,empirical.quantile)
  } else if(type=="adapt"){
    if(length(pars0)!=2){stop("beta is of the size 2x1")}
    G <- pars0[2]
    q <- AdaptLoop(pars0[1],X,theta,G,empirical.quantile)
  } else if(type=="asymslope"){ 
    if(length(pars0)!=4){stop("beta is of the size 4x1")}
    q <- AsymSlopeLoop(pars0,X,empirical.quantile)
  } else if(type=="rq"){ 
    intercept <- pars0[1]
    slope <- pars0[2]
    q <- intercept+slope*X 
  } else {
    stop('Please choose a different type of CAViaR specification (see documentation for more details)')
  }
  hit.stat <- theta-(Y<q)
  rq.stat <- t(hit.stat)%*%(Y-q)
  if (out.val==1){
    if(quant.type=="var"){
      z$VaR <- -q
    } else {
      z$VaR <- q
    }
    z$hit.stat <- hit.stat
    z$rq.stat <- rq.stat
    z$residuals <- Y - q
    z$fitted.values <- q
    class(z) <- "dynquant"
  } else if(out.val == 0){
    z <- as.numeric(rq.stat)
  }
  return(z)
  
}


beta.poli.w <- function(nlag,coeff.1,coeff.2) {
  temp <- seq(0,1,length=nlag)
  if (coeff.1==1){
    weights <- (1-temp)^(coeff.2-1)
  } else {
    weights <- (1-temp)^(coeff.2-1) * (temp^(coeff.1-1))
  }
  weights <- weights/sum(weights)
  return(list(weights = weights))
}

get.u.best.init <- function(dt,theta,z,type,is.midas,starting.vals,quant.type,empirical.quantile) {
  rq.stat <- compute.u.quantile(starting.vals[dt,],theta,z,type,is.midas,out.val=0,quant.type,empirical.quantile)
  return(rq.stat = rq.stat)
}

get.mv.best.init <- function(dt,theta,z,is.midas,starting.vals,quant.type,empirical.quantile) {
  rq.stat <- compute.mv.quantile(starting.vals[dt,],theta,z,is.midas,out.val=0,quant.type,empirical.quantile)
  return(rq.stat = rq.stat)
}

get.u.optimal <- function(dt,param,theta,z,type,is.midas,quant.type,empirical.quantile,opt.method=c("nelder-mead","cma-es"),...) {
  opt.method <- match.arg(opt.method)   
  opt.options <- list(...)
  if(is.null(opt.options$rep)) {rep = 10} else {rep <- opt.options$rep}
  if(is.null(opt.options$constraints)) {constraints <- get.u.constraints(type,is.midas)} else {constraints <- opt.options$constraints}
  if(opt.method=="nelder-mead"){
    est <-  optimx::optimx(param[dt,],compute.u.quantile,theta=theta,z=z,type=type,is.midas=is.midas,out.val=0,quant.type=quant.type,empirical.quantile=empirical.quantile,method = c("Nelder-Mead")) 
    for (i in 1:rep) {
      est <-  optimx::optimx(as.numeric(est[1:length(param[dt,])]),compute.u.quantile,theta=theta,z=z,type=type,is.midas=is.midas,out.val=0,quant.type=quant.type,empirical.quantile=empirical.quantile,method = c("Nelder-Mead"))
    }
    est <-  optimx::optimx(as.numeric(est[1:length(param[dt,])]),compute.u.quantile,theta=theta,z=z,type=type,is.midas=is.midas,out.val=0,quant.type=quant.type,empirical.quantile=empirical.quantile,method = c("Nelder-Mead"))
  } else if(opt.method=="cma-es"){
    est <- cmaes::cma_es(param[dt,],compute.u.quantile,theta=theta,z=z,type=type,is.midas=is.midas,out.val=0,quant.type=quant.type,empirical.quantile=empirical.quantile,lower=constraints$LB,upper=constraints$UB)
  }
  
  return(est = est)
}

get.mv.optimal <- function(dt,param,theta,p,z,is.midas,quant.type,empirical.quantile,opt.method=c("nelder-mead","cma-es"),...) {
  opt.method <- match.arg(opt.method)   
  opt.options <- list(...)
  
  if(is.null(opt.options$rep)) {rep = 10} else {rep <- opt.options$rep}
  if(is.null(opt.options$constraints)) {constraints <- get.mv.constraints("cav",is.midas,p)} else {constraints <- opt.options$constraints}
  if(opt.method=="nelder-mead"){ #
    est <-  optimx::optimx(param[dt,],compute.mv.quantile,theta=theta,z=z,is.midas=is.midas,out.val=0,quant.type=quant.type,empirical.quantile=empirical.quantile,method = c("Nelder-Mead")) 
    for (i in 1:rep) {
      est <-  optimx::optimx(as.numeric(est[1:length(param[dt,])]),compute.mv.quantile,theta=theta,z=z,is.midas=is.midas,out.val=0,quant.type=quant.type,empirical.quantile=empirical.quantile,method = c("Nelder-Mead"))
    }
    est <-  optimx::optimx(as.numeric(est[1:length(param[dt,])]),compute.mv.quantile,theta=theta,z=z,is.midas=is.midas,out.val=0,quant.type=quant.type,empirical.quantile=empirical.quantile,method = c("Nelder-Mead"))
  } else if(opt.method=="cma-es"){
    est <- cmaes::cma_es(param[dt,],compute.mv.quantile,theta=theta,z=z,is.midas=is.midas,out.val=0,quant.type=quant.type,empirical.quantile=empirical.quantile,lower=constraints$LB,upper=constraints$UB)
  }
  
  return(est = est)
}

var.cov.mat.est <- function(z, type) {
  beta  <- z$est.params
  theta <- z$theta
  y     <- z$y
  tobs  <- length(y)
  VaR   <- z$VaR
  # Compute the quantile residuals.
  residuals <- y + VaR
  
  # Set up the bandwidth for the k-Nearest neighbor estimator.
  SortedRes <- sort(abs(residuals), decreasing =FALSE)
  
  if (length(residuals) < 40) {
    k <- 30
  }  else if (theta == 0.01) {
    k <- 40
  }  else if (theta == 0.025) {
    k <- 50
  }  else {
    k <- 60
  }
  
  BANDWIDTH <- SortedRes[k]
  
  # Initialize matrices.
  derivative1 <- array(0, c(tobs,1))
  derivative2 <- array(0, c(tobs,1))
  derivative3 <- array(0, c(tobs,1))
  derivative4 <- array(0, c(tobs,1))
  
  D <- array(0, c(length(beta),length(beta)))
  A <- D
  t <- 0 # counter 
  ##################################
  
  # model: caviar
  if (type == "cav") {
    gradient <- array(0, c(tobs,3))
    
    for (i in 2:tobs){
      derivative1[i] <- 1 + beta[2] * derivative1[i-1]
      derivative2[i] <- VaR[i-1] + beta[2] * derivative2[i-1]
      derivative3[i] <- beta[2] * derivative3[i-1] + y[i-1]
      
      gradient[i,]  <- c(derivative1[i], derivative2[i], derivative3[i])
      gtg           <- gradient[i,]%*%t(gradient[i,])
      A             <- A + gtg     
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gtg
      }
      
    }
  }
  
  # model: asymmetric slope
  if (type == "asymslope") {
    gradient <- array(0, c(tobs, 4))
    
    for (i in 2:tobs){
      derivative1[i] <- 1 + beta[2] * derivative1[i-1]
      derivative2[i] <- VaR[i-1] + beta[2] * derivative2[i-1]
      derivative3[i] <- beta[2]*derivative3[i-1] + y[i-1]*(y[i-1]>0)
      derivative4[i] <- beta[2]*derivative4[i-1] - y[i-1]*(y[i-1]<0)
      
      gradient[i,]  <- c(derivative1[i], derivative2[i], derivative3[i], derivative4[i])
      gtg           <- gradient[i,]%*%t(gradient[i,])
      A             <- A + gtg     
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gtg
      }
      
    }
  }
  
  # model: igarch
  if (type == "igarch") {
    gradient <- array(0, c(tobs, 3))
    
    for (i in 2:tobs){
      derivative1[i] <- (1 + 2*beta[2]*VaR[i-1]*derivative1[i-1]) / (2*VaR[i])
      derivative2[i] = (VaR[i-1]^2 + 2*beta[2]*VaR[i-1]*derivative2[i-1]) / (2*VaR[i]);
      derivative3[i] = (2*beta[2]*VaR[i-1]*derivative3[i] + y[i-1]^2) / (2*VaR[i])
      
      gradient[i,]  <- c(derivative1[i], derivative2[i], derivative3[i])
      gtg           <- gradient[i,]%*%t(gradient[i,])
      A             <- A + gtg
      
      if (abs(residuals[i]) <= BANDWIDTH){
        t <- t+1
        D <- D + gtg
      }
      
    }
  }
  
  # model: adaptive
  if (type == "adapt") {
    gradient <- array(0, c(tobs, 2))
    
    for (i in 2:tobs){
      indicator <-  exp(beta[2] * (y[i-1] + VaR[i-1]))
      if (indicator == Inf){
        derivative1[i] <- derivative1[i-1]
      } else {
        derivative1[i] <- derivative1[i-1] + 1/(1 + indicator) - theta + beta[1] * beta[2] * indicator * derivative1[i-1] / (1 + indicator)^2
        derivative2[i] <- beta[1] * indicator * (y[i-1] + VaR[i-1]) / (1 + indicator)^2
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
  z$A         <- A/tobs
  z$D         <- D/(2*BANDWIDTH*tobs)
  z$gradient  <- gradient
  z$VCmatrix  <- theta*(1-theta)*solve(z$D)%*%z$A%*%solve(z$D)/tobs
  
  return(z = z)
}

dq.stat <- function(z, type, quant.type, ...){
  opt <- list(...)
  if(length(opt$lags) == 0){lags = 4} else {lags = opt$lags}
  
  z         <- var.cov.mat.est(z, type)
  y         <- z$y
  VaR       <- z$VaR
  gradient  <- z$gradient
  hit.stat       <- z$hit.stat
  BANDWIDTH <- z$BANDWIDTH
  D         <- z$D
  theta     <- z$theta
  tobs      <- length(y)
  if(quant.type == "var"){
    residuals <- y + VaR
  } else {
    residuals <- y - VaR
  }
  constant    <- matrix(1,tobs-lags,1)
  hit.stat         <- hit.stat[seq(lags+1,tobs, by = 1)]
  VaRforecast <- VaR[seq(lags+1,tobs, by = 1)]
  yLag        <- y[seq(lags,tobs-1, by = 1)]
  VaRforecastlag <- VaR[seq(lags,tobs-1, by = 1)]
  Z           <- array(0, c(tobs-lags, lags))
  
  for (s in 1:lags){
    Z[,s]              <- hit.stat[seq(s,tobs-(lags+1-s), by = 1)]
  }
  
  Xout <- cbind(constant, VaRforecast, Z) # Instruments for the out of sample test.
  Xin  <- Z                               # Instruments for the in sample test.
  
  XHNABLA   <- array(0,c(dim(Xin)[2],length(z$est.params)))
  NABLA     <- gradient[seq(lags+1,tobs, by=1),]
  
  #-- Estimate the matrices that enter the In Sample DQ test --#
  
  for(j in 2:(tobs-lags)){
    if(abs(residuals[j]) <= BANDWIDTH){
      XHNABLA <- XHNABLA + Xin[j,]%*%t(gradient[j,])
    }
  }
  
  XHNABLA <- XHNABLA/(2*BANDWIDTH*tobs)
  M         <- t(Xin) - XHNABLA%*%solve(D)%*%t(NABLA)
  
  #-- Compute the DQ tests --#
  z$dq.statIn  <- (t(hit.stat)%*%Xin%*%solve(M%*%t(M))%*%t(Xin)%*%hit.stat) / (theta*(1-theta))
  z$dq.statOut <- (t(hit.stat)%*%Xout%*%solve(t(Xout)%*%Xout)%*%t(Xout)%*%hit.stat) / (theta*(1-theta))
  
  z$DQin      <- 1 - pchisq(z$dq.statIn, dim(Xin)[2])   # Compute the P-value of the in     sample DQ test.
  z$DQout     <- 1 - pchisq(z$dq.statOut, dim(Xout)[2]) # Compute the P-value of the out of sample DQ test.
  
  return(z)
}

get.u.constraints <- function(type,is.midas){
  const  <- NULL
  if(type=="cav"){
    const$LB[1] <- -100  
    const$LB[2] <- -100#-1.5  
    const$LB[3] <- -100
    const$UB[1] <-  100  
    const$UB[2] <-  100#1.5  
    const$UB[3] <-  100
  } else if(type=="igarch"){
    const$LB[1] <- -100  
    const$LB[2] <-    0  
    const$LB[3] <-    0
    const$UB[1] <-  100  
    const$UB[2] <-  100  
    const$UB[3] <-  100
  } else if(type=="adapt"){
    const$LB[1] <- -100  
    const$LB[2] <- 0  
    const$UB[1] <  100  
    const$UB[2] <- 100000        
  } else if(type=="asymslope"){ 
    const$LB[1] <- -100  
    const$LB[2] <- -100#-1.5  
    const$LB[3] <- -100  
    const$LB[4] <- -100
    const$UB[1] <-  100  
    const$UB[2] <-  100#1.5  
    const$UB[3] <-  100  
    const$UB[4] <-  100     
  } else if(type=="rq"){ 
    const$LB[1] <- -100  
    const$LB[2] <- -100  
    const$UB[1] <-  100  
    const$UB[2] <-  100  
  } else {
    stop('CAViaR type either not implemented or does not exist')
  }
  if(is.midas) {
    const$LB <- c(const$LB,1)
    const$UB <- c(const$UB,100)
  }
  return(const)
  
}

get.mv.constraints <- function(type,is.midas,p){
  const  <- NULL
  if(type=="cav"){
    const$LB[1:p] <- -100           # constant
    const$LB[(p+1):(p^2+p)] <- -100#-1.5     # autoregressive part
    const$LB[(p^2+p+1):(2*p^2+p)] <- -100 # realized part
    const$UB[1:p] <-  100  
    const$UB[(p+1):(p^2+p)] <- 100#1.5  
    const$UB[(p^2+p+1):(2*p^2+p)] <-  100
  } else {
    stop('MVMQCAViaR type either not implemented or does not exist')
  }
  for(i in 1:p) {
    if(is.midas[i]) {
      const$LB <- c(const$LB,1)
      const$UB <- c(const$UB,100)
    }
  }
  return(const)
  
}

get.u.initial <- function(type,is.midas,num.test=1e3){
  if(type=="cav"){
    starting.vals     <- array(NA, c(num.test,3))
    starting.vals[,1] <- runif(num.test, min =  0, 1) 
    starting.vals[,2] <- runif(num.test, min =  0, 1) 
    starting.vals[,3] <- runif(num.test, min =  0, 1)
  } else if(type=="igarch"){
    starting.vals     <- array(NA, c(num.test,3))
    starting.vals[,1] <- runif(num.test, min =  0, 1) 
    starting.vals[,2] <- runif(num.test, min =  0, 1) 
    starting.vals[,3] <- runif(num.test, min =  0, 1)
  } else if(type=="adapt"){
    starting.vals     <- array(NA, c(num.test,2))
    starting.vals[,1] <- runif(num.test, min =  0,     1) 
    starting.vals[,2] <- runif(num.test, min =  0,  1000) 
  } else if(type=="asymslope"){ 
    starting.vals     <- array(NA, c(num.test,4))
    starting.vals[,1] <- runif(num.test, min =  0, 1) 
    starting.vals[,2] <- runif(num.test, min =  0, 1) 
    starting.vals[,3] <- runif(num.test, min =  0, 1)   
    starting.vals[,4] <- runif(num.test, min =  0, 1)  
  } else if(type=="rq"){
    starting.vals     <- array(NA, c(num.test,2))
    starting.vals[,1] <- runif(num.test, min =  0, 1) 
    starting.vals[,2] <- runif(num.test, min =  0, 1)
  }
  else {
    stop('CAViaR type either not implemented or does not exist')
  }
  if (is.midas) {
    starting.vals     <- cbind(starting.vals,NA)
    starting.vals[,ncol(starting.vals)] <- runif(num.test, min =  1,   30)
  }
  return(starting.vals)
}

get.mv.initial <- function(pars0,type,is.midas,p,num.test=1e3){
  pars <- pars0[1:(2*p^2+p)]
  pars.midas <- pars0[(2*p^2+p+1):length(pars0)]
  if(type=="cav"){
    draw <- matrix(runif(num.test*(2*p^2+p),min=-0.1,0.1),nrow=num.test,ncol=(2*p^2+p))
    starting.vals <- sweep(draw,2,pars,"+")
  } else {
    stop('MVMQCAViaR type either not implemented or does not exist')
  }
  draw <- matrix(runif(num.test*p,min=1,20),nrow=num.test,ncol=p)
  starting.vals <- cbind(starting.vals, sweep(draw,2,pars.midas,"+"))
  return(starting.vals)
}


if (.Platform$OS.type=="windows"){
  foldername <- "C:/Users/Utilisateur/Documents/GitHub/inflation-risk-nowcasting/caviar"
  foldername <- "C:/Users/User/Documents/GitHub/inflation-risk-nowcasting/caviar"
} else {
  foldername <- "~/Documents/GitHub/inflation-risk-nowcasting/caviar"
}
setwd(foldername)
rm(foldername)
Rcpp::sourceCpp("caviar_routines.cpp")
setwd("../")
require("quantreg")
require('foreach')
require("optimx")

compute_caviar <- function(pars0, tau, y, emp_quant, which_caviar, return_quantile = FALSE, is_midas = FALSE){
  if(is_midas){
    theta <- pars0[length(pars0)]
    if(theta<1)
      theta <- 1
    if(theta>100)
      theta <- 100
    y <- as.numeric(y%*%rbeta_w(theta,dayLag = dim(y)[2]))
    pars0 <- pars0[-length(pars0)]
  }
  if(which_caviar == "cav")
    q <- CaviarLoop(pars0, y, emp_quant) # 3-pars
  if(which_caviar == "asym")
    q <- AsymSlopeLoop(pars0, y, emp_quant) # 4-par
  
  hit_stat <- tau-(y<q)
  rq_stat <- t(hit_stat)%*%(y-q)/length(y)
  if(return_quantile) {
    ret <- q
  } else {
    ret <- as.numeric(rq_stat)
  }
  ret
}

generate_init <- function(which_caviar, how_many = 10, sim_hm = 1e4, tau, y, emp_quant, is_midas, seed = 100){
  set.seed(seed)
  if (which_caviar == "cav"){
    num_par <- 3
    try_parameters <- matrix(runif(sim_hm*num_par, min = -1, max = 1), ncol = num_par, nrow = sim_hm)
  }
  if (which_caviar == "asym"){
    num_par <- 4
    try_parameters <- matrix(runif(sim_hm*num_par, min = -1, max = 1), ncol = num_par, nrow = sim_hm)
  }
  if (is_midas)
    try_parameters <- cbind(try_parameters, runif(sim_hm, min = 1, max = 70))
  rq_stats <- apply(try_parameters, 1, compute_caviar, tau, y, emp_quant, which_caviar, is_midas = is_midas)
  res <- cbind(rq_stats, try_parameters)
  res <- res[order(rq_stats),]
  res <- res[1:how_many, ]
  init_params <- res[,-1]
  init_params
}

fit_caviar <- function(y, tau, emp_quant = NULL, which_caviar = "cav", is_midas = FALSE, how_many = 50, sim_hm = 1e4, verbose = FALSE, outx = NULL, seed = 100) {
  # check if midas data was inputed:
  if(!is.null(dim(y)))
    if(!is_midas){
      is_midas <- TRUE
      message("MIDAS data was inputed for y but the 'is_midas' choice was set to FALSE. 'is_midas' is reset to TRUE and model is implemented with MIDAS")
    }
  if (is.null(emp_quant)){
    # fit QAR(1) and use the first fitted value:
    if(is_midas){
      y_tmp <- colMeans(y)
      fit0 <- quantreg::rq(y_tmp[-1]~y_tmp[-length(y_tmp)], tau = tau)
      emp_quant <- as.numeric(fit0$fitted.values[1])
    } else {
      fit0 <- quantreg::rq(y[-1]~y[-length(y)], tau = tau)
      emp_quant <- as.numeric(fit0$fitted.values[1])
    }
  }
  # get initial values:
  beta0 <- generate_init(which_caviar = which_caviar, how_many = how_many, sim_hm = sim_hm, tau = tau, y = y, emp_quant = emp_quant, is_midas = is_midas, seed = seed)
  if(!is_midas){
    lower <- c(-Inf,-1+.Machine$double.eps,-Inf)
    upper <- c(Inf,1-.Machine$double.eps,Inf)
    output <- foreach(k = 1:how_many, .packages = c("optimx")) %do% { 
      fit <- suppressWarnings(optimx(beta0[k, ], compute_caviar, lower = lower, upper = upper, tau = tau, y = y, emp_quant = emp_quant, which_caviar = which_caviar, is_midas = FALSE, method = c("L-BFGS-B")))
      fit
    }
  } else {
    lower <- c(-Inf,-1+.Machine$double.eps,-Inf,1)
    upper <- c(Inf,1-.Machine$double.eps,Inf,70)
    output <- foreach(k = 1:how_many, .packages = c("optimx")) %do% { 
      fit <- suppressWarnings(optimx(beta0[k, ], compute_caviar, lower = lower, upper = upper, tau = tau, y = y, emp_quant = emp_quant, which_caviar = which_caviar, is_midas = TRUE, method = c("L-BFGS-B")))
      fit
    }
  }
  vals <- numeric(how_many)
  for (i in seq(how_many))
    vals[i] <- output[[i]]$value
  fit <- output[[which(vals==min(vals))]]
  if(which_caviar=="cav")
    num_par <- 3
  if(which_caviar == "asym")
    num_par <- 4
  fit <- as.numeric(fit)
  coef <- fit[1:num_par]
  theta <- NULL
  if(is_midas)
    theta <- fit[num_par+1]
  if(is_midas)
    fitted.values <- compute_caviar(c(coef,theta), tau = tau, y = y, emp_quant = emp_quant, which_caviar = which_caviar, return_quantile = TRUE, is_midas = TRUE)
  if(!is_midas)
    fitted.values <- compute_caviar(coef, tau = tau, y = y, emp_quant = emp_quant, which_caviar = which_caviar, return_quantile = TRUE, is_midas = FALSE)
  residuals <- y - fitted.values 
  rq <- compute_caviar(coef, tau = tau, y = y, emp_quant = emp_quant, which_caviar = which_caviar, return_quantile = FALSE)
  if(!is.null(outx)){
    if(is_midas)
      forecast <- compute_caviar(c(coef,theta), tau = tau, y = c(0,outx$y), emp_quant = fitted.values[length(fitted.values)], which_caviar = which_caviar, return_quantile = TRUE, is_midas = TRUE)[2]
    if(!is_midas)
      forecast <- compute_caviar(coef, tau = tau, y = c(0,outx$y), emp_quant = fitted.values[length(fitted.values)], which_caviar = which_caviar, return_quantile = TRUE, is_midas = FALSE)[2]
  } else {
    forecast <- NULL
  }
  fit <- list(coef = coef, theta = theta, fitted.values = fitted.values, residuals = residuals, rq = rq, emp_quant = emp_quant, forecast = forecast)
  return(fit)
}

fit_mvmqcaviar <- function(y1, y2, tau1, tau2, emp_quant1 = NULL, emp_quant2 = NULL, is_midas = FALSE, how_many = 10, sim_hm = 1e4, outx = NULL, seed = 100){
  
  # check if midas data was inputed:
  if(!is.null(dim(y2)))
    if(!is_midas){
      is_midas <- TRUE
      message("MIDAS data was inputed for y2 but the 'is_midas' choice was set to FALSE. 'is_midas' is reset to TRUE and model is implemented with MIDAS")
    }
  
  
  # fit univariate CAViaR model (note: MIDAS variable should be ordered second, if used)
  fit1 <- fit_caviar(y = y1, tau = tau1, emp_quant = emp_quant1, which_caviar = "cav", is_midas = FALSE)
  fit2 <- fit_caviar(y = y2, tau = tau2, emp_quant = emp_quant2, which_caviar = "cav", is_midas = is_midas)
  
  c0 <- c(fit1$coef[1], fit2$coef[1])
  a0 <- cbind(c(fit1$coef[2], 0), c(0, fit2$coef[2]))
  r0 <- cbind(c(fit1$coef[3], 0), c(0, fit2$coef[3]))
  
  emp_quant1 <- fit1$emp_quant
  emp_quant2 <- fit2$emp_quant
  params0 = c(c0,as.vector(a0),as.vector(r0))
  lower <- c(-Inf,-Inf,-1+.Machine$double.eps,-Inf,-1+.Machine$double.eps,-Inf,-Inf,-Inf,-Inf,-Inf)
  upper <- c( Inf, Inf, 1-.Machine$double.eps, Inf, 1-.Machine$double.eps, Inf, Inf, Inf, Inf, Inf)
  if(is_midas){
    params0 <- c(params0, fit2$theta)
    lower <- c(lower,1)
    upper <- c(upper,70)
  }
  
  mvmqobj <- function(params, y1, y2, tau1, tau2, emp_quant1, emp_quant2, return_quantile = FALSE){
    c <- params[1:2]
    params <- params[-c(1:2)]
    a <- params[1:4]
    params <- params[-c(1:4)]
    r <- params[1:4]
    if (length(params)>4)
      y2 <- as.numeric(y2%*%rbeta_w(params[5], dayLag = dim(y2)[2]))
    a <- cbind(c(a[1],a[2]),c(a[3],a[4]))
    r <- cbind(c(r[1],r[2]),c(r[3],r[4]))
    rq <- MvMqCavLoop(c, a, r, y = cbind(y1, y2), c(tau1, tau2), c(emp_quant1, emp_quant2))$rq
    quant <- MvMqCavLoop(c, a, r, y = cbind(y1, y2), c(tau1, tau2), c(emp_quant1, emp_quant2))$q
    if (return_quantile){
      ret <- quant
    } else {
      ret <- rq
    }
    ret
  }
  
  set.seed(100)
  try_parameters <- t(replicate(sim_hm, params0))+(matrix(rnorm(sim_hm*length(params0)),nrow=sim_hm,ncol=length(params0))/50)
  # check if MIDAS parameter satisfies constraint:
  idx <- try_parameters[,dim(try_parameters)[2]]<1
  if(any(idx))
    try_parameters[idx,dim(try_parameters)[2]] <- 1
  
  rq_stats <- apply(try_parameters, 1, mvmqobj, y1 = y1, y2 = y2, tau1 = tau1, tau2 = tau2, emp_quant1 = emp_quant1, emp_quant2 = emp_quant2)
  res <- cbind(rq_stats, try_parameters)
  res <- res[order(rq_stats),]
  res <- res[1:how_many, ]
  init_params <- res[,-1]
  
  fit <- NULL
  for (i in seq(how_many))
    fit[[i]] <- suppressWarnings(optimx(init_params[i, ], mvmqobj, lower = lower, upper = upper, y1 = y1, y2 = y2, tau1 = tau1, tau2 = tau2, emp_quant1 = emp_quant1, emp_quant2 = emp_quant2))
  
  vals <- numeric(how_many)
  for (i in seq(how_many))
    vals[i] <- fit[[i]]$value
  fit <- fit[[which(vals==min(vals))]]
  if(is_midas){
    coef <- as.numeric(fit[1:11])
    rq <- as.numeric(fit[12])
  } else {
    coef <- as.numeric(fit[1:10])
    rq <- as.numeric(fit[11])
  }
  fitted.values <- mvmqobj(coef, y1 = y1, y2 = y2, tau1 = tau1, tau2 = tau2, emp_quant1 = emp_quant1, emp_quant2 = emp_quant2, return_quantile = TRUE)
  colnames(fitted.values) <- c("y1","y2")
  fitted.values <- data.frame(fitted.values)
  if(is_midas){
    residuals <- fitted.values[,1] - y1
  } else {
    residuals <- fitted.values-cbind(y1,y2)
    colnames(residuals) <- c("y1","y2")
    residuals <- data.frame(residuals)
  }
  c <- as.numeric(coef[1:2])
  a <- as.numeric(coef[3:6])
  r <- as.numeric(coef[7:10])
  theta <- coef[11]
  a <- cbind(c(a[1],a[2]),c(a[3],a[4]))
  r <- cbind(c(r[1],r[2]),c(r[3],r[4]))
  colnames(a) <- c("a11", "a12")
  rownames(a) <- c("a21", "a22")
  colnames(r) <- c("r11", "r12")
  rownames(r) <- c("r21", "r22")
  
  if (!is.null(outx)){
    if(is_midas)
      forecast <- mvmqobj(coef, y1 = c(0,outx$y1), y2 = rbind(numeric(dim(outx$y2)[2]), outx$y2), tau1 = tau1, tau2 = tau2, emp_quant1 = fitted.values[dim(fitted.values)[1],1], emp_quant2 = fitted.values[dim(fitted.values)[1],2], return_quantile = TRUE)[2, ]
    if(!is_midas) 
      forecast <- mvmqobj(coef, y1 = c(0,outx$y1), y2 = c(0,outx$y2), tau1 = tau1, tau2 = tau2, emp_quant1 = fitted.values[dim(fitted.values)[1],1], emp_quant2 = fitted.values[dim(fitted.values)[1],2], return_quantile = TRUE)[2, ]
    forecast <- matrix(forecast,ncol=2,nrow=1)
    colnames(forecast) <- c("y1", "y2")
    forecast <- data.frame(forecast)
  } else {
    forecast <- NULL
  }
  return(list(coef = coef, rq = rq, c = c,a = a,r = r, theta = theta, fitted.values = fitted.values, forecast = forecast))
}

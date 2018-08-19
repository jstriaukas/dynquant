mixed.freq.data.single <- function(data.refdate,data.x,data.xdate,x.lag,horizon,est.start,est.end,disp.flag=TRUE) {
  # complete data
  mask.na <- !is.na(data.refdate)
  data.refdate <- data.refdate[mask.na]
  mask.na <- !is.na(data.x)
  data.x <- data.x[mask.na]
  data.xdate <- data.xdate[mask.na]
  data.refdate <- as.Date(data.refdate)
  data.x <- as.vector(data.x)
  data.xdate <- as.Date(data.xdate)
  
  data.refdate.vec <- as.Date(data.refdate)
  data.xdate.vec <- as.Date(data.xdate)
  
  est.start <- as.Date(est.start)
  est.end <- as.Date(est.end)
  
  data.refdate.vec <- date.vec(data.refdate.vec)
  data.xdate.vec <- date.vec(data.xdate.vec)
  data.refdate.vec <- matrix(unlist(data.refdate.vec),nrow=length(data.refdate))
  data.xdate.vec <- matrix(unlist(data.xdate.vec),nrow=length(data.xdate))
  data.refdate.num <- data.refdate
  data.xdate.num <- data.xdate
  
  date.format = c('year(s)','month(s)','day(s)','hour(s)','minute(s)','second(s)')
  period.ref <- data.freq(data.refdate.vec)$period 
  unit.ref <- data.freq(data.refdate.vec)$unit 
  period.x <- data.freq(data.xdate.vec)$period  
  unit.x <- data.freq(data.xdate.vec)$unit
  
  ref.lag <- lag.num(1,period.ref,unit.ref)
  x.lag <- lag.num(x.lag,period.x,unit.x)
  horizon <- lag.num(horizon,period.x,unit.x)
  if (ref.lag < 0){
    stop('ref.lag cannot be negative.')
  }
  if (x.lag < 0) {
    stop('x.lag cannot be negative')
  }
  # Minimum and maximum dates that data support
  min.date.ref <- data.refdate.num[ref.lag+1]
  min.date.x <- data.xdate.num[max(1,x.lag+horizon)]
  if (min.date.ref > min.date.x){
    min.date <- min.date.ref
  } else {
    min.date <- min.date.x
  }
  max.date.ref <- data.refdate.num[length(data.refdate.num)]
  max.date.x = data.xdate.num[length(data.xdate.num)]
  if (horizon < 0){
    max.date.x <- data.xdate.vec[dim(data.xdate.vec)[1],]
    max.date.x[unit.x] <- max.date.x[unit.x] + period.x * horizon
    max.date.x <- ISOdate(max.date.x[1],max.date.x[2],max.date.x[3],max.date.x[4],max.date.x[5],max.date.x[6])
    max.date.x <- as.Date(max.date.x)
  }
  if (max.date.ref > max.date.x){
    max.date <- max.date.x
  } else {
    max.date <- max.date.ref
  }
  # Check and set default sample period
  if (is.null(est.start)){
    est.start <- min.date
  } else { if(est.start < min.date) {warning('Start date cannot be earlier than possible due to lagged regressors. Reset start date to most recent possible: ', paste(min.date))
    est.start <- min.date}
  }
  if (is.null(est.end)){
    est.end <- max.date
  } else { if(est.end > max.date) {warning('Terminal date cannot be later than largest date account for lags. Reset to largest date possible: ', paste(max.date))
    est.end <- max.date}
  }
  # Construct reference date data
  tol <- 1e-10
  loc.start <- min(which((data.refdate.num >= est.start-tol) == TRUE))
  loc.end <- min(which((data.refdate.num >= est.end-tol) == TRUE)) 
  est.refdate <- data.refdate.num[loc.start:loc.end]
  
  loc.forecast.end <- min(which((data.refdate.num >= max.date-tol) == TRUE))
  if(loc.end+1<=loc.forecast.end){
    out.refdate <- data.refdate.num[seq(loc.end+1,loc.forecast.end,by=1)]
    n.forecast <- length(out.refdate)
  } else {
    out.refdate <- NULL
    n.forecast <- length(out.refdate)
  }
  nobs <- loc.end - loc.start + 1
  
  est.x <- est.xdate <- matrix(NaN,nrow=nobs,ncol=x.lag) 
  for (t in 1:nobs){
    loc <- min(which((data.xdate.num >= est.refdate[t]-tol) == TRUE)) 
    if (is.null(loc)) {
      loc <- length(data.xdate.num)
    }
    
    if(loc-horizon > length(data.x)){    
      nobs <- t - 1
      est.ref = est.ref[seq(1,nobs,1)]
      est.refdate = est.refdate[seq(1,nobs,1)]
      est.lag.ref = est.lag.ref[seq(1,nobs,1)]
      est.lag.refdate = est.lag.refdate[seq(1,nobs,1)]
      est.x = est.x[seq(1,nobs,1)]
      est.xdate = est.xdate[seq(1,nobs,1)]
      max.date = est.refdate[length(est.refdate)]
      warning('Horizon is a large negative number. Observations are further truncated to max date possible: ', paste(max.date))
      break
    } else  {      
      est.x[t,] <- data.x[seq(loc-horizon,loc-horizon-x.lag+1,-1)]
      est.xdate[t,] <- data.xdate.num[seq(loc-horizon,loc-horizon-x.lag+1,-1)]
    }
  }
  if(loc.end+1<=loc.forecast.end){
    out.x <- out.xdate <- matrix(NaN,nrow=n.forecast,ncol=x.lag) 
    for(t in 1:n.forecast){
      loc <- min(which((data.xdate.num >= out.refdate[t]-tol) == TRUE))  
      if (is.null(loc)) {
        loc <- length(data.xdate.num)
      }
      
      if(loc-horizon > length(data.x)){      
        n.forecast <- t - 1
        out.ref = out.ref[seq(1,n.forecast,1)]
        out.refdate = out.refdate[seq(1,n.forecast,1)]
        out.lag.ref = out.lag.ref[seq(1,n.forecast,1)]
        out.lag.refdate = out.lag.refdate[seq(1,n.forecast,1)]
        out.x = out.x[seq(1,n.forecast,1)]
        out.xdate = out.xdate[seq(1,n.forecast,1)]
        break
      } else {
        out.x[t,] <- data.x[seq(loc-horizon,loc-horizon-x.lag+1,-1)] 
        out.xdate[t,] <- data.xdate.num[seq(loc-horizon,loc-horizon-x.lag+1,-1)] 
      } 
    }
  } else {
    out.x <- out.xdate <- NULL
  }
  
  if (disp.flag==T){
    # Display mixed frequency data
    cat('Frequency of Reference Date:',period.ref,date.format[unit.ref], "\n") 
    cat('Frequency of Data X:',period.x,date.format[unit.x], "\n") 
    cat('Start Date: ', paste(est.start), "\n") 
    cat('Terminal Date: ', paste(est.end), "\n") 
    
    # Display timeframe of mixed frequency regression
    cat('Mixed frequency data structure time frame:', "\n") 
    for(m in c(1,2,nobs)){
      cat("\n")
      cat(paste('Ref date(',as.Date(est.refdate[m],origin="1970-01-01"),')`s on: ', sep="")) 
      if (x.lag == 1) {
        cat(paste(' X(',as.Date(est.xdate[m],origin="1970-01-01"),')`s', sep="")) 
      }
      if (x.lag == 2){
        cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
      }
      if (x.lag == 3) {
        cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,2],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
      }
      if (x.lag > 3){
        cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,2],origin="1970-01-01"),')`s ... X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
      }
    }
  }
  output = list(est.refdate = est.refdate, est.x = est.x, est.xdate = est.xdate,
                out.refdate = out.refdate, out.x = out.x, out.xdate = out.xdate,
                x.lag = x.lag,min.date = min.date, max.date = max.date)
  return(output)
}



mixed.freq.data <- function(data.y,data.ydate,data.x,data.xdate,x.lag,y.lag,horizon,est.start,est.end,disp.flag=TRUE) {
  # complete data
  mask.na <- !is.na(data.y)
  data.y <- data.y[mask.na]
  data.ydate <- data.ydate[mask.na]
  mask.na <- !is.na(data.x)
  data.x <- data.x[mask.na]
  data.xdate <- data.xdate[mask.na]
  data.y <- as.vector(data.y)
  data.ydate <- as.Date(data.ydate)
  data.x <- as.vector(data.x)
  data.xdate <- as.Date(data.xdate)
  
  
  data.ydate.vec <- as.Date(data.ydate)
  data.xdate.vec <- as.Date(data.xdate)
  
  
  est.start <- as.Date(est.start)
  est.end <- as.Date(est.end)
  
  
  data.ydate.vec <- date.vec(data.ydate.vec)
  data.xdate.vec <- date.vec(data.xdate.vec)
  data.ydate.vec <- matrix(unlist(data.ydate.vec),nrow=length(data.ydate))
  data.xdate.vec <- matrix(unlist(data.xdate.vec),nrow=length(data.xdate))
  data.ydate.num <- data.ydate
  data.xdate.num <- data.xdate
  
  date.format = c('year(s)','month(s)','day(s)','hour(s)','minute(s)','second(s)')
  period.y <- data.freq(data.ydate.vec)$period 
  unit.y <- data.freq(data.ydate.vec)$unit 
  period.x <- data.freq(data.xdate.vec)$period  
  unit.x <- data.freq(data.xdate.vec)$unit
  
  
  y.lag <- lag.num(y.lag,period.y,unit.y)
  x.lag <- lag.num(x.lag,period.x,unit.x)
  horizon <- lag.num(horizon,period.x,unit.x)
  if (y.lag < 0){
  stop('y.lag cannot be negative.')
  }
  if (x.lag < 0) {
  stop('x.lag cannot be negative')
  }
  
  # Minimum and maximum dates that data support
  min.date.y <- data.ydate.num[y.lag+1]
  min.date.x <- data.xdate.num[max(1,x.lag+horizon)]
  if (min.date.y > min.date.x){
    min.date <- min.date.y
  } else {
    min.date <- min.date.x
  }
  max.date.y <- data.ydate.num[length(data.ydate.num)]
  max.date.x = data.xdate.num[length(data.xdate.num)]
  if (horizon < 0){
  max.date.x <- data.xdate.vec[dim(data.xdate.vec)[1],]
  max.date.x[unit.x] <- max.date.x[unit.x] + period.x * horizon
  max.date.x <- ISOdate(max.date.x[1],max.date.x[2],max.date.x[3],max.date.x[4],max.date.x[5],max.date.x[6])
  max.date.x <- as.Date(max.date.x)
  }
  if (max.date.y > max.date.x){
  max.date <- max.date.x
  } else {
    max.date <- max.date.y
  }
  # Check and set default sample period
  if (is.null(est.start)){
  est.start <- min.date
  } else { if(est.start < min.date) {warning('Start date cannot be earlier than possible due to lagged regressors. Reset start date to most recent possible.')
    est.start <- min.date}
  }
  if (is.null(est.end)){
  est.end <- max.date
  } else { if(est.end > max.date) {warning('Terminal date cannot be later than largest date account for lags. Reset to largest date.')
    est.end <- max.date}
  }
  # Construct Y data
  tol <- 1e-10
  loc.start <- min(which((data.ydate.num >= est.start-tol) == TRUE))
  loc.end <- min(which((data.ydate.num >= est.end-tol) == TRUE)) 
  est.y <- data.y[loc.start:loc.end]
  est.ydate <- data.ydate.num[loc.start:loc.end]
  
  loc.forecast.end <- min(which((data.ydate.num >= max.date-tol) == TRUE))
  if(loc.end+1<=loc.forecast.end){
    out.y <- data.y[seq(loc.end+1,loc.forecast.end,by=1)]
    out.ydate <- data.ydate.num[seq(loc.end+1,loc.forecast.end,by=1)]
    n.forecast <- length(out.y)
  } else {
    out.y <- out.ydate <- NULL
    n.forecast <- length(out.y)
  }
  nobs <- loc.end - loc.start + 1
  # Construct lagged Y data
  est.lag.y <- est.lag.ydate <- matrix(NaN,nrow=nobs,ncol=y.lag)
  for (m in 1:y.lag){
  est.lag.y[,m] <- data.y[seq(loc.start-m,loc.end-m,1)]
  est.lag.ydate[,m] <- data.ydate.num[seq(loc.start-m,loc.end-m,1)]
  }
  
  if(loc.end+1<=loc.forecast.end){
  out.lag.y <- out.lag.ydate <- matrix(NaN,nrow=n.forecast,ncol=y.lag) 
  for (m in 1:y.lag){
  out.lag.y[,m] <- data.y[seq(loc.end-m,loc.forecast.end-m,1)]
  out.lag.ydate[,m] <- data.ydate.num[seq(loc.end-m,loc.forecast.end-m,1)]
  }
  } else {
    out.lag.y <- out.lag.ydate <- NULL
  }
  
  est.x <- est.xdate <- matrix(NaN,nrow=nobs,ncol=x.lag) 
  for (t in 1:nobs){
  loc <- min(which((data.xdate.num >= est.ydate[t]-tol) == TRUE)) 
  if (is.null(loc)) {
  loc <- length(data.xdate.num)
  }
  
  if(loc-horizon > length(data.x)){    
  nobs <- t - 1
  est.y = est.y[seq(1,nobs,1)]
  est.ydate = est.ydate[seq(1,nobs,1)]
  est.lag.y = est.lag.y[seq(1,nobs,1)]
  est.lag.ydate = est.lag.ydate[seq(1,nobs,1)]
  est.x = est.x[seq(1,nobs,1)]
  est.xdate = est.xdate[seq(1,nobs,1)]
  max.date = est.ydate[length(est.ydate)]
  warning('Horizon is a large negative number. Observations are further truncated to max date possible')
  break
  } else  {      
    est.x[t,] <- data.x[seq(loc-horizon,loc-horizon-x.lag+1,-1)]
    est.xdate[t,] <- data.xdate.num[seq(loc-horizon,loc-horizon-x.lag+1,-1)]
    }
  }
  if(loc.end+1<=loc.forecast.end){
  out.x <- out.xdate <- matrix(NaN,nrow=n.forecast,ncol=x.lag) 
  for(t in 1:n.forecast){
  loc <- min(which((data.xdate.num >= out.ydate[t]-tol) == TRUE))  
  if (is.null(loc)) {
    loc <- length(data.xdate.num)
  }
  
  if(loc-horizon > length(data.x)){      
  n.forecast <- t - 1
  out.y = out.y[seq(1,n.forecast,1)]
  out.ydate = out.ydate[seq(1,n.forecast,1)]
  out.lag.y = out.lag.y[seq(1,n.forecast,1)]
  out.lag.ydate = out.lag.ydate[seq(1,n.forecast,1)]
  out.x = out.x[seq(1,n.forecast,1)]
  out.xdate = out.xdate[seq(1,n.forecast,1)]
  break
  } else {
    out.x[t,] <- data.x[seq(loc-horizon,loc-horizon-x.lag+1,-1)] 
    out.xdate[t,] <- data.xdate.num[seq(loc-horizon,loc-horizon-x.lag+1,-1)] 
  } 
  }
  } else {
    out.x <- out.xdate <- NULL
  }
  
  if (disp.flag==T){
  # Display mixed frequency data
  cat('Frequency of Data Y:',period.y,date.format[unit.y], "\n") 
  cat('Frequency of Data X:',period.x,date.format[unit.x], "\n") 
  cat('Start Date: ', paste(est.start), "\n") 
  cat('Terminal Date: ', paste(est.end), "\n") 
  
  # Display timeframe of mixed frequency regression
  cat('Mixed frequency regression time frame:', "\n") 
  for(m in c(1,2,nobs)){
    cat("\n")
    cat(paste('Reg Y(',as.Date(est.ydate[m],origin="1970-01-01"),')`s on: ', sep="")) 
  if (y.lag == 1){
    cat(paste('Y(',as.Date(est.lag.ydate[m,],origin="1970-01-01"),')`s', sep="")) 
  }
  if (y.lag == 2){
    cat(paste('Y(',as.Date(est.lag.ydate[m,],origin="1970-01-01"),')`s Y(',as.Date(est.lag.ydate[m,dim(est.lag.ydate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  if (y.lag > 2){
    cat(paste('Y(',as.Date(est.lag.ydate[m,],origin="1970-01-01"),')`s ... Y(',as.Date(est.lag.ydate[m,dim(est.lag.ydate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  if (x.lag == 1) {
    cat(paste(' X(',as.Date(est.xdate[m],origin="1970-01-01"),')`s', sep="")) 
  }
  if (x.lag == 2){
    cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  if (x.lag == 3) {
    cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,2],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  if (x.lag > 3){
    cat(paste(' X(',as.Date(est.xdate[m,1],origin="1970-01-01"),')`s X(',as.Date(est.xdate[m,2],origin="1970-01-01"),')`s ... X(',as.Date(est.xdate[m,dim(est.xdate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  }
  }
  output = list(est.y = est.y,est.ydate = est.ydate, est.x = est.x, est.xdate = est.xdate,
                  est.lag.y = est.lag.y, est.lag.ydate = est.lag.ydate,
                  out.y = out.y, out.ydate = out.ydate, out.x = out.x, out.xdate = out.xdate,
                  out.lag.y = out.lag.y, out.lag.ydate = out.lag.ydate, x.lag = x.lag, y.lag = y.lag,
                  min.date = min.date, max.date = max.date)
  return(output)
}

data.freq <- function(DateVec) {
  # data.freq: Identify data frequency
  #
  # Input Arguments:
  #  %
  # DateVec: T-by-6 R vector format data: [year,month,day,hour,min,sec]
  #
  # Output Arguments:
  #  %
  #% period: length of two consecutive dates
  #%
  #% unit: unit of length measure
#%       o 1 = year
#%       o 2 = month
#%       o 3 = day
#%       o 4 = hour
#%       o 5 = minutes
#%       o 6 = seconds
#%
#% Notes:
#  %
#% Frequency   period   unit
#% yearly         1      1  
#% semiannual     6      2 
#% quarterly      3      2 
#% monthly        1      2 
#% biweekly       14     3
#% weekly         7      3
#% daily          1      3
#% hourly         1      4
#% minutely       1      5
#% secondly       1      6

DateDiff <- diff(DateVec)

# Check annual or lower frequency
modeUse = mode(DateDiff[,1])$dataMode
if(modeUse >= 1) {
period <- modeUse
unit <- 1
return(list(period = period,unit = unit))
}

# Check monthly frequency, quarter = 3 months, semiannual = 6 months
modeUse <- mode(DateDiff[,2])$dataMode
mask <- isTRUE(modeUse < 0)
modeUse[mask] <- modeUse[mask] + 12
if(modeUse >= 1){
period <- modeUse
unit <- 2
return(list(period = period,unit = unit))
}

# Check daily frequency, week = 7 days, biweekly = 14 days
modeUse <- mode(DateDiff[,3])$dataMode
mask <- isTRUE(modeUse < 0)
modeUse[mask] <- modeUse[mask] +30
if(modeUse >= 1){
  period <- modeUse
  unit <- 3
  return(list(period = period,unit = unit))
}

# Check hourly frequency

modeUse <- mode(DateDiff[,4])$dataMode
mask <- isTRUE(modeUse < 0)
modeUse[mask] <- modeUse[mask] + 24
if(modeUse >= 1){
  period <- modeUse
  unit <- 4
  return(list(period = period,unit = unit))
}

# Check minutely frequency
modeUse <- mode(DateDiff[,5])$dataMode
mask <- isTRUE(modeUse < 0)
modeUse[mask] <- modeUse[mask] + 60
if(modeUse >= 1){
  period <- modeUse
  unit <- 5
  return(list(period = period,unit = unit))
}


# Check secondly frequency
elapse <- diff.time.mf(DateVec[2:dim(DateVec)[1],6],DateVec[1:dim(DateVec)[1]-1,6],origin = "1970-01-01",units = "secs")
period <- mean(elapse)
unit   <- 6



}

mode <- function(data) {

# mode: mode (most frequent value) of a vector
#
# Input Arguments:
#
# data: a vector of data

# Output Arguments:

# dataMode: mode of the vector
#
# countMax: frequency count at the mode


nobs <- length(data)
data <- sort(data)
count <- 1
countMax <- 1
dataMode = data[1]

for (t in 2:nobs){
  if(data[t]==data[t-1]){
    count <- count + 1
  }  
  if(data[t]!=data[t-1]){ 
      if(count > countMax){
        countMax <- count 
        dataMode <- data[t-1]
        } 
    count <- 1
  }
} # end of for
if (count > countMax) {
countMax <- count
dataMode  <- data[nobs]
}
return(list(dataMode = dataMode, countMax = countMax))
} # end of mode

date.vec <- function(s) {
  mat <- matrix(0, nrow = length(s), ncol = 6 )
  a <- as.POSIXlt(as.Date(s, "1970-01-01"))
  mat[,1] <- a$year + 1900
  mat[,2] <- a$mon + 1
  mat[,3] <- a$mday 
  mat[,4] <- a$hour
  mat[,5] <- a$min
  mat[,6] <- a$sec
  
  list(mat)
}

lag.num <- function(x.lag,period,unit){

if(is.numeric(x.lag)==T && is.atomic(x.lag)==T && length(x.lag) == 1L) {
    return(x.lag)
}
  multiplier <- as.double(substr(x.lag, start = 1, stop = nchar(x.lag)-1))
  if (is.na(multiplier)) {
    stop('The description of lags cannot be recognized. The format should be 3m, 1q, etc')
  }
  #Convert multiplier to daily frequency (business days)
  ndaysPerYear <- 260
  ndaysPerQuarter <- 65
  ndaysPerMonth <- 22
  nhoursPerDay <- 8
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'y'){
    multiplier <- multiplier * ndaysPerYear
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'q'){
    multiplier <- multiplier * ndaysPerQuarter
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'q'){
    multiplier <- multiplier * ndaysPerQuarter
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'q'){
    multiplier <- multiplier * ndaysPerQuarter
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'm'){
    multiplier <- multiplier * ndaysPerMonth
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'd'){
    multiplier <- multiplier * 1
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 'h'){
    multiplier <- multiplier / nhoursPerDay
  }
  if(substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 's'){
    multiplier <- multiplier /  (nhoursPerDay*60*60)
  }
  if(unit == 1) {
    x.lag <- round(multiplier / (ndaysPerYear * period))
  }
  if(unit == 2) {
    x.lag <- round(multiplier / (ndaysPerMonth * period))
  }
  if(unit == 3) {
    x.lag <- round(multiplier / period)
  }
  if(unit == 4) {
    x.lag <- round(multiplier / (period / nhoursPerDay))
  }
  if(unit == 5) {
    x.lag <- round(multiplier / (period / nhoursPerDay / 60))
  }
  if(unit == 6) {
    x.lag <- round(multiplier / (period / nhoursPerDay / 60 / 60))
  }
  return(x.lag)
}

diff.time.mf <- function(time1, time2, origin, units = c("auto", "secs", "mins", "hours", "days", "weeks")) {
  if (missing(origin)) {
    time1 <- as.POSIXct(time1)
    time2 <- as.POSIXct(time2)
  }
  else {
    time1 <- as.POSIXct(time1, origin = origin)
    time2 <- as.POSIXct(time2, origin = origin)
  }
  z <- unclass(time1) - unclass(time2)
  attr(z, "tzone") <- NULL
  units <- match.arg(units)
  if (units == "auto") 
    units <- if (all(is.na(z))) 
      "secs"
  else {
    zz <- min(abs(z), na.rm = TRUE)
    if (!is.finite(zz) || zz < 60) 
      "secs"
    else if (zz < 3600) 
      "mins"
    else if (zz < 86400) 
      "hours"
    else "days"
  }
  switch(units, secs = .difftime(z, units = "secs"), mins = .difftime(z/60, units = "mins"), 
         hours = .difftime(z/3600, units = "hours"), 
         days = .difftime(z/86400, units = "days"), 
         weeks = .difftime(z/(7 * 86400), units = "weeks"))
}






rm(list = ls())

EstY <- read.csv("~/Documents/GitHub/dynquant/R/EstY.txt", sep="")
EstYdate.imp <- read.table("~/Documents/GitHub/dynquant/R/EstYdate.csv", quote="\"", comment.char="")
EstX <- read.csv("~/Documents/GitHub/dynquant/R/EstX.txt", sep="")
EstXdate.imp <- read.table("~/Documents/GitHub/dynquant/R/EstXdate.csv", quote="\"", comment.char="")

EstYdate <- EstXdate <- NULL

for (i in 1:dim(EstYdate.imp)[1]) {
  EstYdate[i] <- as.Date(EstYdate.imp[i,], origin = "1899-12-30") 
}
EstYdate <- as.Date(EstYdate, origin = "1970-01-01")
for (j in 1:dim(EstXdate.imp)[1]) {
  EstXdate[j] <- as.Date(EstXdate.imp[j,], origin = "1899-12-30") 
}
EstXdate <- as.Date(EstXdate, origin = "1970-01-01")


DataY <- EstY 
DataYdate <- EstYdate
DataX <- EstX
DataXdate <- EstXdate

estStart <-  as.Date('1964-02-01')
estEnd   <- as.Date('2017-09-01')

xlag <- 66
ylag <- 1
horizon <- '1d'

data <- MixFreqData(DataY,DataYdate,DataX,DataXdate,xlag,ylag,horizon,estStart,estEnd)

MixFreqData <- function(DataY,DataYdate,DataX,DataXdate,xlag,ylag,horizon,estStart,estEnd,dispFlag=TRUE) {
  # complete data
  mask.na <- !is.na(DataY)
  DataY <- DataY[mask.na]
  DataYdate <- DataYdate[mask.na]
  mask.na <- !is.na(DataX)
  DataX <- DataX[mask.na]
  DataXdate <- DataXdate[mask.na]
  DataY <- as.vector(DataY)
  DataYdate <- as.Date(DataYdate)
  DataX <- as.vector(DataX)
  DataXdate <- as.Date(DataXdate)
  
  
  DataYdateVec <- as.Date(DataYdate)
  DataXdateVec <- as.Date(DataXdate)
  
  
  estStart <- as.Date(estStart)
  estEnd <- as.Date(estEnd)
  
  
  DataYdateVec <- DatevecV(DataYdateVec)
  DataXdateVec <- DatevecV(DataXdateVec)
  DataYdateVec <- matrix(unlist(DataYdateVec),nrow=length(DataYdate))
  DataXdateVec <- matrix(unlist(DataXdateVec),nrow=length(DataXdate))
  DataYdateNum <- DataYdate
  DataXdateNum <- DataXdate
  Yinfo <- DateFreq(DataYdateVec)
  periodY <-Yinfo$period  
  unitY <- Yinfo$unit
  
  dateFormat = c('year(s)','month(s)','day(s)','hour(s)','minute(s)','second(s)')
  Xinfo <- DateFreq(DataXdateVec)
  periodX <- Xinfo$period  
  unitX <- Xinfo$unit
  ylag <- lagNum(ylag,periodY,unitY)
  xlag <- lagNum(xlag,periodX,unitX)
  horizon <- lagNum(horizon,periodX,unitX)
  if (ylag < 0){
  stop('ylag cannot be negative.')
  }
  if (xlag < 0) {
  stop('xlag cannot be negative')
  }
  
  # Minimum and maximum dates that data support
  minDateY <- DataYdateNum[ylag+1]
  minDateX <- DataXdateNum[max(1,xlag+horizon)]
  if (minDateY > minDateX){
    minDate <- minDateY
  } else {
    minDate <- minDateX
  }
  maxDateY <- DataYdateNum[length(DataYdateNum)]
  maxDateX = DataXdateNum[length(DataXdateNum)]
  if (horizon < 0){
  maxDateX <- DataXdateVec[dim(DataXdateVec)[1],]
  maxDateX[unitX] <- maxDateX[unitX] + periodX * horizon
  maxDateX <- ISOdate(maxDateX[1],maxDateX[2],maxDateX[3],maxDateX[4],maxDateX[5],maxDateX[6])
  maxDateX <- as.Date(maxDateX)
  }
  if (maxDateY > maxDateX){
  maxDate <- maxDateX
  } else {
    maxDate <- maxDateY
  }
  # Check and set default sample period
  if (is.null(estStart)){
  estStart <- minDate
  } else { if(estStart < minDate) {warning('Start date cannot be earlier than possible due to lagged regressors. Reset start date to most recent possible.')
    estStart <- minDate}
  }
  if (is.null(estEnd)){
  estEnd <- maxDate
  } else { if(estEnd > maxDate) {warning('Terminal date cannot be later than largest date account for lags. Reset to largest date.')
    estEnd <- maxDate}
  }
  # Construct Y data
  tol <- 1e-10
  locStart <- min(which((DataYdateNum >= estStart-tol) == TRUE))
  locEnd <- min(which((DataYdateNum >= estEnd-tol) == TRUE)) 
  EstY <- DataY[locStart:locEnd]
  EstYdate <- DataYdateNum[locStart:locEnd]
  
  locForecastEnd <- min(which((DataYdateNum >= maxDate-tol) == TRUE))
  if(locEnd+1<=locForecastEnd){
    OutY <- DataY[seq(locEnd+1,locForecastEnd,by=1)]
    OutYdate <- DataYdateNum[seq(locEnd+1,locForecastEnd,by=1)]
    nforecast <- length(OutY)
  } else {
    OutX <- OutXdate <- NULL
    nforecast <- length(OutY)
  }
  nobs <- locEnd - locStart + 1
  # Construct lagged Y data
  EstLagY <- EstLagYdate <- matrix(NaN,nrow=nobs,ncol=ylag)
  for (m in 1:ylag){
  EstLagY[,m] <- DataY[seq(locStart-m,locEnd-m,1)]
  EstLagYdate[,m] <- DataYdateNum[seq(locStart-m,locEnd-m,1)]
  }
  OutLagY <- OutLagYdate <- matrix(NaN,nrow=nforecast,ncol=ylag) 
  for (m in 1:ylag){
  OutLagY[,m] <- DataY[seq(locEnd-m,locForecastEnd-m,1)]
  OutLagYdate[,m] <- DataYdateNum[seq(locEnd-m,locForecastEnd-m,1)]
  }
  
  EstX <- EstXdate <- matrix(NaN,nrow=nobs,ncol=xlag) 
  for (t in 1:nobs){
  loc <- min(which((DataXdateNum >= EstYdate[t]-tol) == TRUE)) 
  if (is.null(loc)) {
  loc <- length(DataXdateNum)
  }
  
  if(loc-horizon > length(DataX)){    
  nobs <- t - 1
  EstY = EstY[seq(1,nobs,1)]
  EstYdate = EstYdate[seq(1,nobs,1)]
  EstLagY = EstLagY[seq(1,nobs,1)]
  EstLagYdate = EstLagYdate[seq(1,nobs,1)]
  EstX = EstX[seq(1,nobs,1)]
  EstXdate = EstXdate[seq(1,nobs,1)]
  maxDate = EstYdate[length(EstYdate)]
  warning('Horizon is a large negative number. Observations are further truncated to max date possible')
  break
  } else  {      
    EstX[t,] <- DataX[seq(loc-horizon,loc-horizon-xlag+1,-1)]
    EstXdate[t,] <- DataXdateNum[seq(loc-horizon,loc-horizon-xlag+1,-1)]
    }
  }
  if(locEnd+1<=locForecastEnd){
  OutX <- OutXdate <- matrix(NaN,nrow=nforecast,ncol=xlag) 
  for(t in 1:nforecast){
  loc <- min(which((DataXdateNum >= OutYdate[t]-tol) == TRUE))  
  if (is.null(loc)) {
    loc <- length(DataXdateNum)
  }
  
  if(loc-horizon > length(DataX)){      
  nforecast <- t - 1
  OutY = OutY[seq(1,nforecast,1)]
  OutYdate = OutYdate[seq(1,nforecast,1)]
  OutLagY = OutLagY[seq(1,nforecast,1)]
  OutLagYdate = OutLagYdate[seq(1,nforecast,1)]
  OutX = OutX[seq(1,nforecast,1)]
  OutXdate = OutXdate[seq(1,nforecast,1)]
  break
  } else {
    OutX[t,] <- DataX[seq(loc-horizon,loc-horizon-xlag+1,-1)] 
    OutXdate[t,] <- DataXdateNum[seq(loc-horizon,loc-horizon-xlag+1,-1)] 
  } 
  }
  } else {
    OutX <- OutXdate <- NULL
  }
  
  if (dispFlag==T){
  # Display mixed frequency data
  cat('Frequency of Data Y:',periodY,dateFormat[unitY], "\n") 
  cat('Frequency of Data X:',periodX,dateFormat[unitX], "\n") 
  cat('Start Date: ', paste(estStart), "\n") 
  cat('Terminal Date: ', paste(estEnd), "\n") 
  
  # Display timeframe of mixed frequency regression
  cat('Mixed frequency regression time frame:', "\n") 
  for(m in c(1,2,nobs)){
    cat("\n")
    cat(paste('Reg Y(',as.Date(EstYdate[m],origin="1970-01-01"),')`s on: ', sep="")) 
  if (ylag == 1){
    cat(paste('Y(',as.Date(EstLagYdate[m,],origin="1970-01-01"),')`s', sep="")) 
  }
  if (ylag == 2){
    cat(paste('Y(',as.Date(EstLagYdate[m,],origin="1970-01-01"),')`s Y(',as.Date(EstLagYdate[m,dim(EstLagYdate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  if (ylag > 2){
    cat(paste('Y(',as.Date(EstLagYdate[m,],origin="1970-01-01"),')`s ... Y(',as.Date(EstLagYdate[m,dim(EstLagYdate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  if (xlag == 1) {
    cat(paste(' X(',as.Date(EstXdate[m],origin="1970-01-01"),')`s', sep="")) 
  }
  if (xlag == 2){
    cat(paste(' X(',as.Date(EstXdate[m,1],origin="1970-01-01"),')`s X(',as.Date(EstXdate[m,dim(EstXdate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  if (xlag == 3) {
    cat(paste(' X(',as.Date(EstXdate[m,1],origin="1970-01-01"),')`s X(',as.Date(EstXdate[m,2],origin="1970-01-01"),')`s X(',as.Date(EstXdate[m,dim(EstXdate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  if (xlag > 3){
    cat(paste(' X(',as.Date(EstXdate[m,1],origin="1970-01-01"),')`s X(',as.Date(EstXdate[m,2],origin="1970-01-01"),')`s ... X(',as.Date(EstXdate[m,dim(EstXdate)[2]],origin="1970-01-01"),')`s', sep=""))  
  }
  }
  }
  Output = list(Est = EstY,EstYdate = EstYdate, EstX = EstX, EstXdate = EstXdate,
                  EstLagY = EstLagY, EstLagYdate = EstLagYdate,
                  OutY = OutY, OutYdate = OutYdate, OutX = OutX, OutXdate = OutXdate,
                  OutLagY = OutLagY, OutLagYdate = OutLagYdate, Xlag = xlag, Ylag = ylag,
                  MinDate = minDate, MaxDate = maxDate)
  return(Output)
}

DateFreq <- function(DateVec) {
  # DateFreq: Identify data frequency
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

DatevecV <- function(s) {
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

lagNum <- function(xlag,period,unit){

if(is.numeric(xlag)==T && is.atomic(xlag)==T && length(xlag) == 1L) {
    return(xlag)
}
  multiplier <- as.double(substr(xlag, start = 1, stop = nchar(xlag)-1))
  if (is.na(multiplier)) {
    stop('The description of lags cannot be recognized. The format should be 3m, 1q, etc')
  }
  #Convert multiplier to daily frequency (business days)
  ndaysPerYear <- 260
  ndaysPerQuarter <- 65
  ndaysPerMonth <- 22
  nhoursPerDay <- 8
  if(substr(xlag, start = nchar(xlag), stop = nchar(xlag)) == 'y'){
    multiplier <- multiplier * ndaysPerYear
  }
  if(substr(xlag, start = nchar(xlag), stop = nchar(xlag)) == 'q'){
    multiplier <- multiplier * ndaysPerQuarter
  }
  if(substr(xlag, start = nchar(xlag), stop = nchar(xlag)) == 'q'){
    multiplier <- multiplier * ndaysPerQuarter
  }
  if(substr(xlag, start = nchar(xlag), stop = nchar(xlag)) == 'q'){
    multiplier <- multiplier * ndaysPerQuarter
  }
  if(substr(xlag, start = nchar(xlag), stop = nchar(xlag)) == 'm'){
    multiplier <- multiplier * ndaysPerMonth
  }
  if(substr(xlag, start = nchar(xlag), stop = nchar(xlag)) == 'd'){
    multiplier <- multiplier * 1
  }
  if(substr(xlag, start = nchar(xlag), stop = nchar(xlag)) == 'h'){
    multiplier <- multiplier / nhoursPerDay
  }
  if(substr(xlag, start = nchar(xlag), stop = nchar(xlag)) == 's'){
    multiplier <- multiplier /  (nhoursPerDay*60*60)
  }
  if(unit == 1) {
    xlag <- round(multiplier / (ndaysPerYear * period))
  }
  if(unit == 2) {
    xlag <- round(multiplier / (ndaysPerMonth * period))
  }
  if(unit == 3) {
    xlag <- round(multiplier / period)
  }
  if(unit == 4) {
    xlag <- round(multiplier / (period / nhoursPerDay))
  }
  if(unit == 5) {
    xlag <- round(multiplier / (period / nhoursPerDay / 60))
  }
  if(unit == 6) {
    xlag <- round(multiplier / (period / nhoursPerDay / 60 / 60))
  }
  return(xlag)
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





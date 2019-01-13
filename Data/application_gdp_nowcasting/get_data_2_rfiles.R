rm(list = ls(all = TRUE)) # Clear library
require('alfred')
require('R.matlab')


## sort monthly data ##
wd <- '/Users/striaukas/Documents/GitHub/dynquant/Data/application_gdp_nowcasting'
setwd(wd)
foldername <- paste(wd, "/matfiles_montly/", sep = "")
list.files <- list.files(path = foldername)
N <- length(list.files)
for (j in 1:N) {
  filename <- paste(wd, "/matfiles/", list.files[j], sep = "")
  dat <- readMat(filename)
  filename.save <- paste(wd, "/rfiles_monthly/", substr(list.files[j], 1, nchar(list.files[j])-4), ".Rdata", sep = "")
  save(dat, file = filename.save)
}

filename <- paste(wd, "/matfiles_daily/financials.mat", sep = "")
dat <- readMat(filename)
filename.save <- paste(wd, "/rfiles_daily/financials.Rdata", sep = "")
save(dat, file = filename.save)

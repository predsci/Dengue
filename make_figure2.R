##
## Wavelet analysis and plot for Dengue data from Sri Lanka 
## The code plots it for all the regions but only the first page is used in the manuscript
## Output plot filename: LK.wavelet.by.state.new.pdf
##
rm(list=ls())

require(DICE)

library("WaveletComp")
library("Hmisc")
library("zoo")
library("RColorBrewer")

mod_name = c(NAME_2="LK")

##
## Fitting/Forecasting Method either compartmental models ('mech') or statistical ('stat')
##
method = 'mech' 


## 
## Define the SARIMA model if method = 'stat' Default is NULL and code uses auto.arima
##
arima_model = list(p = 1, d = 0, q = 0, P = 3, D = 1, Q = 0) # NULL

##
## Optional co-variate for SARIMA. Used ONLY if 'arima_model' is defined by the User and method = 'stat'
##

covar = 'precip' # or 'temp' or 'precip' default is FALSE

##
## Lag time (in cadence of data units) for co-variate
## 
covar_lag = 1 

##
## Disease to model (flu or dengue)
##

disease ='dengue' #'dengue' #'flu'


##
## An integer or string-abbreviation selecting which data source should be used.  A list of data sources appears in DICE_data_sources.  If value is NULL, ## DICE will attempt to auto-choose a data source.
##

data_source = NULL
## 
## Use the MySQl database (dengue) or not (flu)
##

db_opts=list(DICE_db="predsci", CDC_server=FALSE)

##
## Select RegState
##


##
## Start year of the disease season
##

year  = 2010

##
## Spatial level of model/fit data - 2 entire usa, 3 HHS regions, 4 States
## For CDC we only support level 2  
## For GFT we support levels 2 (entire usa) and  3 (a single HHS region)
##

mod_level = 2

##
## Spatial level of data used to fit the model data >= mod_level
##

fit_level = 3

##
## Number of data points that are fitted - default is to fit all of the available data
##
nfit = 12

##
## Coupled (isingle = 0) or Uncoupled (isingle = 1) fit for the fit_level spatial regions
isingle     = 0

## 
## SIR (1) SEIR (2) can also accept string: sir, seir, vsir , vseir or: SIR, SEIR, VSIR, VSEIR
## 
epi_model = 2

##
## Average Infectious Period  (in days)
##
Tg = 10.0

## 
## Model Number (1-5), 1- SH Only, 2 - School Vacation only, 3 - Both SH and SV, 4 - Constant, 5 -  A two value model
## Models 2 and 3 should not used for dengue
## Model 1 also makes little sense as does model 5..
## The code will reset the model to 4 if the user chose 2 or 3 
##
model  = 4

## Select Method for Data Augmentation
## 0 No data augmentation
## 1 Use historic monthly average, NULL model
## 2 Use the most similar Season
da = 0
##
## Use a prior (prior = 1) or no (prior = 0). Currently relevant only for flu
##
prior = 0
##
## Temperature for LLK in MCMC procedure
##
Temp = 100

##
## Number of steps/trials in MCMC chain
##
nMCMC = 2e6

##
## Number of MCMC chains
##
nreal = 1

plot = TRUE

##
## Optional- Name of sub-directory where all the output files/plots will be saved. If set to NULL code will build the directory
## name based on the parameters selection
##
subDir = NULL # 'test'

complete_mydata <- get.DICE.data(data_source = NULL, mod_level = mod_level, fit_level = fit_level,
                                 year = year, model = model, nperiodsFit = nfit, mod_name=mod_name,
                                 RegState = NULL, fit_names='all', isingle = isingle,
                                 db_opts=db_opts, disease = disease, epi_model = epi_model,
                                 method = method, all_years_flag=T, all_cad_clim=T, raw_col=NULL)
                                 

mydata=complete_mydata$mydata
all_years_epi=complete_mydata$all_years_epi
all_season_dates=complete_mydata$all_season_dates
all_years_clim=complete_mydata$all_years_clim

dat = all_years_epi$fit$raw
sh  = all_years_epi$fit$sh
precip = all_years_epi$fit$precip
dates  = all_years_epi$dates
years  = all_years_epi$years
months = all_years_epi$months

nts = dim(dat)[1]
NA_ind = NULL
for (i in 1:nts) {
  if (all(is.na(dat[i,]))) NA_ind = rbind(NA_ind, i)
}
if (!is.null(NA_ind)) {
  dat = dat[-NA_ind,]
  sh  = sh[-NA_ind,]
  precip = precip[-NA_ind,]
  dates  = dates[-NA_ind]
  years  = years[-NA_ind]
  months = months[-NA_ind]

}

## remove the last two rows since they are all NA

# Population weighted centroid
latlon = data.frame(lat=mydata$fit$attr$sedac_lat, lon = mydata$fit$attr$sedac_lon)
pop=mydata$fit$pop
state_names = mydata$fit$name

nstates = mydata$fit$nregions

cases   = rep(0, nstates)
icount  = 0
names(cases) = state_names

for (id in 1:nstates) {
	icount = icount + 1
	cases[icount] = sum(dat[,id], na.rm = TRUE)

}
#

iorder = order(cases, decreasing = TRUE)
cases = cases[iorder]


## For now we take the entire data set - but this can be an input from the User even in the form of month/year for start/end
## If we actually set iend to an earlier date we can score the prediction..
##


noTPts = length(dates)


## Find the index of the start of each year
nyears = length(unique(years))
unqYear = unique(years)
startYear = rep(0,nyears)
for(i in 1:nyears) {
	startYear[i] = which(years == unqYear[i])[1]
}
## Make a 2D plot of the state  level data - just to see how synchronous it is.

stateData = dat
colnames(stateData) = state_names
## convert to a matrix
stateData = as.matrix(stateData)

state.df= data.frame(name=rep(0,nstates),lon=rep(0,nstates),lat=rep(0,nstates))

for (i in 1:nstates) {
	state.df$name[i] = state_names[i]
	state.df$lon[i] = latlon[i,'lon']
	state.df$lat[i] = latlon[i,'lat']
}

## Reorder the state.df and state Data by the number of cases/

#iorder = order(latlon$lat)

latlon = latlon[iorder,]
state.df = state.df[iorder,]
stateData = stateData[,iorder]


range = range(stateData,na.rm=TRUE)


## For LK 2^6 for BR 2^7

npadd  = noTPts - 2^6

if ((npadd%%2) == 0) {
	npadd2p = npadd2m = npadd/2
} else {
	npadd2p = floor(npadd/2)
	npadd2m = floor(npadd/2) + 1
}


## Maybe plot the log
logState = log(stateData)

for (i in 1:nstates) {
	stateData[is.na(stateData[,i]),i] = 0
	logState[,i] = log(stateData[,i]+1)
	logState[,i] = logState[,i] - mean(logState[,i])
	logState[,i] = logState[,i]/sd(logState[,i])
}

lrange = range(logState,na.rm=TRUE)

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))   # make colors
r <- rf(64)

## For LK

my.date = dates # seq(from = as.Date("2010-01-01"), to = as.Date("2017-07-01"), by = "1 month")
#date    = seq(from = as.Date("2008-11-01"), to = as.Date("2018-08-01"), by = "1 month")

dateP = seq(from=my.date[length(my.date)], length = (npadd2p+1), by = '1 month')
dateP = dateP[-1]

dateM = seq(from=my.date[1], length = (npadd2m+1), by = '-1 month')
dateM = rev(dateM)
dateM = dateM[-length(dateM)]

date = c(dateM, my.date, dateP)

date = seq(from=date[1], to=date[length(date)], by = '1 month')

index = which(state.df$name == names(cases)[1])

ctrLog = as.numeric(logState[,index])

y = c(rep(0,npadd2m),ctrLog,rep(0,npadd2p))
ctr.data = data.frame(date = date, x = y)
ctr.w = analyze.wavelet(ctr.data, "x", loess.span = 0, dt = 1, dj = 1/noTPts, lowerPeriod = 4, upperPeriod = 64, make.pval = T, n.sim = 100, verbose=F)
# calculate the average Phase/Angle for the 'center' for periods of 10-14 month
subset = which(ctr.w$Period >= (10) & ctr.w$Period <= (14))

ctr.angle = rep(0,length(y))
ictr = index

my.angle  = rep(0,length(y))

rad2deg = 180 / pi

for (i in 1:length(y)) {
	ctr.angle[i] = mean(ctr.w$Phase[subset,i]) * rad2deg
}

phase.coh = dhf.coh = dist = rep(0, nstates)

names(phase.coh) = names(dhf.coh) = names(dist)= state.df$name

rec.year = rec.state = array(NA,c(length(date),nstates))


pdf(file=paste0(mod_name[[1]],'.wavelet.by.state.new.pdf'),onefile=TRUE,width=15,height=13)

par(mfrow = c(5, 3), mar = c(2, 4, 2, 4))

for (ix in 1:nstates) {


	stateName = state.df$name[ix]
	deng = as.numeric(logState[, ix])

	# replace NA with '1'
	deng[is.na(deng)] <- log(1)
	deng = deng - mean(deng)
	deng = deng/sd(deng)
	# Let's padd the data with zero's until the next 2^n
	x = c(rep(0,npadd2m),deng,rep(0,npadd2p))

	my.data = data.frame(date = date, x = x, y=y)
	my.w = analyze.wavelet(my.data, "x", loess.span = 0, dt = 1, dj = 1/noTPts, lowerPeriod = 4, upperPeriod = 64, make.pval = T, n.sim = 100, verbose=F)

	for (i in 1:length(x)) my.angle[i] = mean(my.w$Phase[subset,i]) * rad2deg

	wt.image(my.w, color.key = "quantile", n.levels = 250, periodlab = "Period (Month)", show.date = T,
	          legend.params = list(lab = " ", mar = 6,lab.line=0.5),
	         graphics.reset = FALSE,main=stateName, date.format = "%F %T", timelab='') #

	my.rec = reconstruct(my.w, timelab = "Year", only.ridge = T, plot.rec = F,verbose=F)

	rec.year[,ix] =  reconstruct(my.w, only.ridge = T, plot.rec = F,verbose=F,sel.period=my.w$Period[subset])$series$x.r
    rec.state[,ix] = my.rec$series$x.r
	if (sum(is.nan(rec.year[,ix]))== length(rec.year[,ix])) rec.year[,ix] = NA
	plot(my.rec$series$x[(npadd2m+1):(npadd2m+1+noTPts)], type = "l", xlab = "", ylab = "log(incidence)", col = "blue",xaxt='no',ylim=range(my.rec$series$x,my.rec$series$x.r),main=stateName)
	lines(my.rec$series$x.r[(npadd2m+1):(npadd2m+1+noTPts)], col = "red")
	axis(1, at = startYear, labels = unique(years))
	legend('topright',legend=c('Data','Reconstructed'),text.col=c('blue','red'),bty='n')

	plot(ctr.angle[(npadd2m+1):(npadd2p+noTPts)],type='l',xlab='',ylab='Phase Angle (degrees)', col='blue',xaxt='no',ylim=c(-200,200),main=stateName)
	lines(my.angle[(npadd2m+1):(npadd2p+noTPts)],type='l',col='red',xlab='',ylab='')
	axis(1, at = startYear, labels = unique(years))
	legend('topright',legend=c(state.df$name[ictr],stateName),text.col=c('blue','red'),bty='n')


}

dev.off()

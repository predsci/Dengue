##
## Plot all available Dengue incidence using the DICE database 
## The list of countries is: c("BR", "MX", "LK", "TH") and "SG" which has weekly data
## Most of the parameters defined below are not used but are required when retrieving the data
## using the DICE functions.
## Output filename for plots is: 'dengue_heat_maps.pdf'
##

rm(list=ls())

require(DICE)
library(RColorBrewer)
library(fields)

## Set default values for all parameters - these are being reset when the script is called
## with differnt values


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
## Select RegState - will be done below  since we are looping over a few countries
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

filename = 'dengue_heat_maps.pdf'
pdf(file = filename, width = 10, height = 22)

par(mar = c(6, 9, 1, 4), mfrow=c(5,1)) 
#set.panel(5,1) # 2X2 matrix of plots


RegState.vec = c("BR", "MX", "LK", "TH")

for (istate in 1:length(RegState.vec)) {
	
	RegState = RegState.vec[istate]
	mod_name = c(NAME_2=RegState)
	cat("\n Processing ", RegState, '\n')
	# Get data. Default behavior is to return the specified season in $mydata as well as all available data in $all_years_epi.
	complete_mydata <- get.DICE.data(data_source = NULL, mod_level = mod_level, fit_level = fit_level,
                                 year = year, model = model, nperiodsFit = nfit, mod_name=mod_name,
                                 RegState = NULL, fit_names='all', isingle = isingle,
                                 db_opts=db_opts, disease = disease, epi_model = epi_model,
                                 method = method, all_years_flag=T, all_cad_clim=T, raw_col=NULL)

	mydata = complete_mydata$mydata # season for modeling
	all_years_epi = complete_mydata$all_years_epi


	data <- all_years_epi$fit$raw
	years <- all_years_epi$years
	months <- all_years_epi$months
	dates  <- all_years_epi$dates

	nts = dim(data)[1]
NA_ind = NULL
for (i in 1:nts) {
  if (all(is.na(data[i,]))) NA_ind = rbind(NA_ind, i)
  if (sum(data[i,]) == 0) NA_ind = rbind(NA_ind, i)
}
if (!is.null(NA_ind)) {
  data = data[-NA_ind,]
  dates  = dates[-NA_ind]
  years  = years[-NA_ind]
  months = months[-NA_ind]

}

	data <- as.matrix(data)
	tiny = 1 # 1e-06
	data <- data + tiny
	
	nx <- dim(data)[1]
	ny <- dim(data)[2]
	
	for (j in 1:ny) {

		data[1:nx, j] = log10(data[1:nx, j])
		#data[1:nx, j] = data[1:nx, j] - mean(data[1:nx, j])
		#data[1:nx, j] = data[1:nx, j]/sd(data[1:nx, j])
		#data[, j] = log10(data[, j])
	}

	# for (j in 1:ny) {

		# data[, j] = log(data[, j])
		# data[, j] = data[, j] - mean(data[, j])
		# data[, j] = data[, j]/sd(data[, j])
	# }


	rf <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) # make colors 'RdBu' 'Spectral'
	r <- rf(128)


	unq.years <- unique(years)
	n.unq <- length(unq.years)
	index <- rep(NA, n.unq)

	for (i in 1:n.unq) {
		index[i] = which(years == unq.years[i])[1]
	}
	every.xth.label = 4
	if (RegState == "MX") 
		every.xth.label = 5
	select.indexes = seq(1, nx, every.xth.label)
	select.x.labels = months[select.indexes]
	x.axis.labels = rep("", nx)
	x.axis.labels[select.indexes] = select.x.labels

	xlab = paste0(years[1], "-", years[nx])
	
	image(x = 1:nx, y = 1:ny, as.matrix(data), zlim =range(data) , ylab = "", xlab = xlab, axes = FALSE, col = r, main = mydata$model$name, xaxt = "n") #range(data)
	box()

	axis(1, at = index, label = years[index]   , las = 2, cex.axis = 0.9)
	axis(2, at = 1:ny , label = mydata$fit$name, las = 1, cex.axis = 0.9)
	image.plot(legend.only = TRUE, zlim = range(data), col = r, legend.width = 0.9, legend.shrink = 0.9, legend.mar = 2, nlevel = 128)

}

RegState = "SG"
cat("\n Processing ", RegState, "\n")
mod_name = c(NAME_2=RegState)
# Get data. Default behavior is to return the specified season in $mydata as well as all available data in $all_years_epi.
	# Get data. Default behavior is to return the specified season in $mydata as well as all available data in $all_years_epi.
	complete_mydata <- get.DICE.data(data_source = NULL, mod_level = mod_level, fit_level = 2,
                                 year = year, model = model, nperiodsFit = nfit, mod_name=mod_name,
                                 RegState = NULL, fit_names='all', isingle = isingle,
                                 db_opts=db_opts, disease = disease, epi_model = epi_model,
                                 method = method, all_years_flag=T, all_cad_clim=T, raw_col=NULL)
                                 
mydata = complete_mydata$mydata # season for modeling
all_years_epi = complete_mydata$all_years_epi

data <- all_years_epi$fit$epi
years <- all_years_epi$years
months <- all_years_epi$months

tiny = 1 #1e-06
data <- data + tiny

nx <- length(data)

data <- as.matrix(data)

data = log10(data)
#data = data - mean(data)
#data = data/sd(data)


unq.years <- unique(years)
n.unq <- length(unq.years)
index <- rep(NA, n.unq)

for (i in 1:n.unq) {
	index[i] = which(years == unq.years[i])[1]
}
xlab = paste0(years[1], "-", years[nx])

plot(x = 1:nx, data, ylab = "", xlab = xlab, col = "red", main = mydata$model$name, lwd = 3, ylim=range(data), type='l', xaxt='n')
axis(1, at = index, label = years[index], las = 2, cex.axis = 1, asp = (nx/ny))

dev.off()

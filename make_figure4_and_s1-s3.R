##
## Make Figure 4.  This code uses an Rds data file for the Mexico Data
## First page of pdf file is used for Figure 4.  
## Pages 2-4 are Figure S1, S2, S3 and S4
rm(list = ls())

library("Hmisc")
library("lubridate")
library("zoo")
library("RColorBrewer")
library("sp")
library("raster")
library("dplyr")
library("GADMTools")
library("rgeos")
library("fields")

##
## select country using ISO-2 code
##

country = "BR"
iso3 = "BRA"

horizon = 1:4

filename = "figure4_and_s1-s3.pdf"
if (country == "BR") {
	pdf(file = filename, width = 15, height = 8)
	par(mfrow = c(4, 7), mar = c(4.5, 3, 1, 0.5))
	plots_per_page = 4 * 7
} else if (country == "TH") {
	pdf(file = filename, width = 25, height = 12)
	par(mfrow = c(6, 13), mar = c(4.5, 3, 1, 0.5))
	plots_per_page = 6 * 13
} else if (country == "MX") {
	pdf(file = filename, width = 15, height = 8)
	par(mfrow = c(4, 5), mar = c(4.5, 3, 1, 0.5))
	plots_per_page = 4 * 5
} else {
	pdf(file = filename, width = 15, height = 8)
	par(mfrow = c(4, 5), mar = c(4.5, 3, 1, 0.5))
	plots_per_page = 4 * 5
}


my.horizon = 4

#country = "TH"
#iso3    = 'THA'
## Read the data 
DENVdat = readRDS("DICE_dengue.Rds")

# For Brazil this  is the state level data

if (country == "BR") {
	attr_ind = DENVdat$attr$ABBV_2 == country & DENVdat$attr$level == 4
	identifiers = DENVdat$attr$identifier[attr_ind]
	state_names = DENVdat$attr$NAME_4[attr_ind]

}

if (country == "TH") {
	attr_ind = DENVdat$attr$ABBV_2 == country & DENVdat$attr$level == 5
	identifiers = DENVdat$attr$identifier[attr_ind]
	state_names = DENVdat$attr$NAME_5[attr_ind]
}

if (country == "MX") {

	attr_ind = DENVdat$attr$ABBV_2 == country & DENVdat$attr$level == 3
	identifiers = DENVdat$attr$identifier[attr_ind]
	state_names = DENVdat$attr$NAME_3[attr_ind]

}

my_pop = DENVdat$SEDAC_pop[identifiers]
my_pop = round(my_pop)
# incidence
my_cols = c("year", "month", "ndays", identifiers)
dat = DENVdat$DENG_monthly[, my_cols]
# misc info

state_lat = DENVdat$attr$lat[attr_ind]
state_lon = DENVdat$attr$lon[attr_ind]
colnames(dat) = c("year", "month", "ndays", state_names)

# climate
cntr_cols = c("year", "day", identifiers)
sh = DENVdat$sh[, cntr_cols]
temperature = DENVdat$temperature[, cntr_cols]
precip = DENVdat$precip[, cntr_cols]


days_per_year = 365
days_per_week = 7
weeks_per_year = 52
month_per_year = 12

nstates = length(state_names)

NA_ind = apply(is.na(DENVdat$DENG_monthly[, identifiers]), 1, all)
dat = dat[!NA_ind, ]

days = dat[, "ndays"]

year = dat[, "year"]

month = dat[, "month"]

noTPts = length(days)

## Find the index of the start of each year
nyears = length(unique(year))
unqYear = unique(year)
startYear = rep(0, nyears)
for (i in 1:nyears) {
	startYear[i] = which(year == unqYear[i])[1]
}

## Take just the cases
stateData = dat[, -c(1, 2, 3)]

## rename columns using the full state names
colnames(stateData) = state_names

## convert to a matrix
stateData = as.matrix(stateData)


state.df = data.frame(name = rep(0, nstates), lon = rep(0, nstates), lat = rep(0, nstates), pop = rep(0, nstates))

state.df$name = state_names
state.df$lon = state_lon
state.df$lat = state_lat
state.df$pop = as.numeric(my_pop)



## Re-order the states by population
iorder = order(state.df$pop, decreasing = TRUE)

state.df = state.df[iorder, ]

## Reorder the stateData
mydf = stateData
for (i in 1:nstates) {
	stateName = colnames(stateData)[i]
	#stateName = substr(stateName,4,100)
	index = which(state.df[, "name"] == stateName)
	mydf[, index] = stateData[, i]
	colnames(mydf)[index] = stateName
}

stateData = mydf

##
## In the case of Mexico - we will remove any state that does not have cases for hald or more of the data time period
##


if (country == "MX") {
	ikeep = NULL
	for (istate in 1:nstates) {
		nzero = sum(stateData[, istate] == 0, na.rm = TRUE)
		if (nzero < (noTPts/2) || length(nzero) == 0) {
			ikeep = cbind(ikeep, istate)
		}
	}
	nstates = length(ikeep)
	stateData = stateData[, ikeep]
	state.df = state.df[ikeep, ]
}



year.min = min(year)
year.max = max(year)

##
## we will define three times: to-date = t_to, delivery-date = t_del and analysis date = t_an
## the t_to specifies that the current forecast will only use cases whose sympotm onset date is equal or less than t_to
## the delivery date specifies that the current forecast will use data up to and including t_del
## the analysis date specified when a given forecast was run
## for now lets set set the t_to to be the same as the delivery date 
## There are 13 years of data so let's use the previous ten to calculate averages 
## mean value for every month
##

nave = 10
nave1 = nave + 1

my.year.vec = unique(year)[nave1:nyears]

month_ave = unique(month)



##
## Define all the ARIMA models we want to test
##


p = 1
d = 0
q = 0
P = 3
D = 0
Q = 0

Period = 12

narima = length(p) * length(d) * length(q) * length(P) * length(D) * length(Q)

for (ihrz in 1:length(horizon)) {
	nplots = 1
	my.horizon = horizon[ihrz]

	cat("\n Horizon Prediction: ", my.horizon, "\n\n")


	deng.obsrv = array(NA, c((month_per_year * length(my.year.vec)), (4 + nstates)))
	deng.upper = array(NA, c((month_per_year * length(my.year.vec)), (4 + nstates), narima))
	deng.lower = array(NA, c((month_per_year * length(my.year.vec)), (4 + nstates), narima))
	deng.frcst = array(NA, c((month_per_year * length(my.year.vec)), (4 + nstates), narima))

	colnames(deng.obsrv) = c("year", "month", "frcstYear", "frcstMonth", state.df$name)

	dimnames(deng.frcst)[[2]] = dimnames(deng.upper)[[2]] = dimnames(deng.lower)[[2]] = c("year", "month", "frcstYear", "frcstMonth", state.df$name)

	arima.names = array(NA, narima)
	iarm = 1
	for (ip in 1:length(p)) {
		for (id in 1:length(d)) {
			for (iq in 1:length(q)) {
				for (iP in 1:length(P)) {
					for (iD in 1:length(D)) {
						for (iQ in 1:length(Q)) {
							arima.names[iarm] = paste0("(", p[ip], d[id], q[iq], ")(", P[iP], D[iD], Q[iQ], ")")
							iarm = iarm + 1
						}
					}
				}
			}
		}
	}

	dimnames(deng.frcst)[[3]] = arima.names


	icount = 1

	for (iyear in 1:length(my.year.vec)) {

		my.year = my.year.vec[iyear]
		cat("Processing Year: ", my.year, "\n")

		for (j in 1:month_per_year) {

			# This is the month that the forecast is made - last month of real data
			my.month = month_ave[j]
			my.row = which(year == my.year & month == my.month)

			cat("Processing Month: ", my.month, "\n")

			## Now do an ARIMF fit on the past ten year starting from now
			istart = my.row - 10 * month_per_year - my.horizon + 1
			istart = max(istart, 1)
			iend = my.row - my.horizon
			##
			## Do ARIMA fits and make a forecast nfrcst forward
##		

			## current data for all states - this is the observed number of cases
			
			cases.now = stateData[my.row, ]
			#cases.now = stateData[(my.row - my.horizon),]
			
			## This is when the forecast is made
			frcstMonth = month[(my.row - my.horizon)]
			frcstYear = year[(my.row - my.horizon)]
			cat("forecast made on", frcstMonth, frcstYear, "\n")

			iarm = 1
			for (ip in 1:length(p)) {

				for (id in 1:length(d)) {

					for (iq in 1:length(q)) {

						for (iP in 1:length(P)) {

							for (iD in 1:length(D)) {

								for (istate in 1:nstates) {

									## Let's fit the log of the data with a mean of zero and sd of one
									
									x = stateData[istart:iend, istate]
									x[is.na(x)] <- 0
									x = log(x + 1)
									mean.x = mean(x)
									x = x - mean.x
									sd.x = sd(x)
									x = x/sd.x
									fit <- try(arima(x, order = c(p[ip], d[id], q[iq]), seasonal = list(order = c(P[iP], D[iD], Q), period = 12), method = "ML"))
									if (isTRUE(class(fit) == "try-error")) 
										next
									result <- try(predict(fit, n.ahead = my.horizon, se.fit = TRUE))

									if (isTRUE(class(result) == "try-error")) 
										next
									if (length(result$pred) == my.horizon) {
										## Need to convert back from log to linear
										my.result = result$pred[my.horizon]

										my.result = exp(my.result * sd.x + mean.x)

										# This is an approximation.  We actually know that exp(sd(log[y])) != sd(y)
										my.sd = exp(fit$sigma2)
										lower = result$pred[my.horizon] * sd.x + mean.x - 1.96 * sqrt(fit$sigma2)
										upper = result$pred[my.horizon] * sd.x + mean.x + 1.96 * sqrt(fit$sigma2)

										deng.frcst[icount, state.df$name[istate], iarm] = my.result
										deng.upper[icount, state.df$name[istate], iarm] = exp(lower)
										deng.lower[icount, state.df$name[istate], iarm] = exp(upper)

									}

								}

								iarm = iarm + 1
							}
						}
					}
				}
			}

			deng.obsrv[icount, state.df$name] = cases.now

			deng.frcst[icount, "year", 1:narima] = my.year
			deng.upper[icount, "year", 1:narima] = my.year
			deng.lower[icount, "year", 1:narima] = my.year
			deng.obsrv[icount, "year"] = my.year

			deng.frcst[icount, "month", 1:narima] = my.month
			deng.upper[icount, "month", 1:narima] = my.month
			deng.lower[icount, "month", 1:narima] = my.month
			deng.obsrv[icount, "month"] = my.month

			deng.frcst[icount, "frcstYear", 1:narima] = frcstYear
			deng.upper[icount, "frcstYear", 1:narima] = frcstYear
			deng.lower[icount, "frcstYear", 1:narima] = frcstYear
			deng.obsrv[icount, "frcstYear"] = frcstYear

			deng.frcst[icount, "frcstMonth", 1:narima] = frcstMonth
			deng.upper[icount, "frcstMonth", 1:narima] = frcstMonth
			deng.lower[icount, "frcstMonth", 1:narima] = frcstMonth
			deng.obsrv[icount, "frcstMonth"] = frcstMonth

			icount = icount + 1

		}

	}


	## Now let's plot
	nts = dim(deng.obsrv)[1]

	for (istate in 1:nstates) {
		my.name = state.df$name[istate]
		upper = deng.upper[, state.df$name[istate], 1]
		lower = deng.lower[, state.df$name[istate], 1]
		obsrv = deng.obsrv[, state.df$name[istate]]
		frcst = deng.frcst[, state.df$name[istate], 1]

		ymin = 0
		ymax = max(upper, lower, obsrv, frcst, na.rm = TRUE)
		xlab = paste0(deng.obsrv[1, "year"], " - ", deng.obsrv[nts, "year"])
		plot(obsrv, type = "l", lwd = 1, col = "black", ylim = c(ymin, ymax), ylab = "Cases", xlab = xlab, xaxt = "n")
		lines(frcst, type = "l", col = "red", lwd = 3, xlab = "", ylab = "")
		axis(1, at = 1:nts, label = deng.obsrv[, "month"], las = 1, cex = 0.5)
		polygon(c(1:nts, rev(1:nts)), c(upper, rev(lower)), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
		legend("topleft", c(my.name, paste0(my.horizon, "- month ahead")), bty = "n")
		nplots = nplots + 1
	}
	if (nplots < plots_per_page) for (i in (nplots):plots_per_page) plot.new()
	
}

dev.off()

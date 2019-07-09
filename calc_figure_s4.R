## This is similar to arima.R but it also calculates the national level using both direct and
## aggregate methods

rm(list=ls())

require(DICE)

##
## Disease to model
##

disease = 'dengue'

##
## An integer or string-abbreviation selecting which data source should be used.  A list of data sources appears in DICE_data_sources.  If value is NULL, ## DICE will attempt to auto-choose a data source.
##

data_source=NULL #"Sri_weekly"
##
## Use the MySQl database (dengue) or not (flu)
##

db_opts=list(DICE_db="predsci", CDC_server=FALSE)

RegState = NULL

mod_name = c(NAME_2="SG")

fit_names = 'all'

mod_level = 2

fit_level = 2

nfit = nperiodsFit = 52

year = 2019

isingle = 1

epi_model = 2

model = 4

prior = 0

Temp = 1

method = 'mech'

complete_mydata <- get.DICE.data(data_source = NULL, mod_level = mod_level, fit_level = fit_level,
                                 year = year, model = model, nperiodsFit = nperiodsFit, mod_name=mod_name,
                                 RegState = RegState, fit_names=fit_names, isingle = isingle,
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
weeks  = all_years_epi$weeks

week_per_year  = 52
month_per_year = 12


# Population weighted centroid
latlon = data.frame(lat=mydata$fit$attr$sedac_lat, lon = mydata$fit$attr$sedac_lon)
pop=mydata$fit$pop
state_names = mydata$fit$name

nstates = mydata$fit$nregions

noTPts = length(dates)

## Find the index of the start of each year
nyears = length(unique(years))
unqYear = unique(years)
startYear = rep(0,nyears)
for(i in 1:nyears) {
	startYear[i] = which(years == unqYear[i])[1]
}

## Take just the cases
stateData = dat

state.df = data.frame(name=rep(0,nstates),lon=rep(0,nstates),lat=rep(0,nstates),pop=rep(0,nstates))

state.df$name = state_names
state.df$lon  = latlon$lon
state.df$lat  = latlon$lat
state.df$pop  = as.numeric(pop)

year.min = min(years)
year.max = max(years)

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
nave1 = nave+1

my.year.vec = unique(years)[nave1:nyears]

month_ave = unique(months)
week_ave = unique(weeks)[1:52]

deng.week.ave  = array(0,c(length(week_ave),nstates))
deng.month.ave = array(0,c(length(month_ave),nstates))
deng.year.ave  = array(0,c(length(week_ave),nstates))


colnames(deng.month.ave) = state.df$name
rownames(deng.month.ave) = month_ave


colnames(deng.week.ave) = colnames(deng.year.ave) = state.df$name
rownames(deng.week.ave) = rownames(deng.year.ave) = week_ave

nhorizon = 10

##
## Define all the ARIMA models we want to test
##


p = 1
d = 0
q = 0
P = 0:3
D = 0:1
Q = 0

narima = length(p) * length(d) * length(q) * length(P) * length(D) * length(Q)
month_per_year = 12

deng.frcst = array(NA, c((nhorizon*week_per_year*length(my.year.vec)),(5+nstates),narima))
deng.obsrv = array(NA, c((nhorizon*week_per_year*length(my.year.vec)),(5+nstates)))
deng.err   = array(NA, c((nhorizon*week_per_year*length(my.year.vec)),(5+nstates),narima))
deng.rel   = array(NA, c((nhorizon*week_per_year*length(my.year.vec)),(5+nstates),narima))
deng.null1 = array(NA, c((nhorizon*week_per_year*length(my.year.vec)),(5+nstates)))
deng.null2 = array(NA, c((nhorizon*week_per_year*length(my.year.vec)),(5+nstates)))

frcst.table = data.frame(current=week_ave,future=array(NA,c(length(week_ave),nhorizon)))

for (i in 2:(nhorizon+1)) {
        frcst.table[,i] = c(i:52, 1:(i-1))
}

colnames(deng.obsrv) = colnames(deng.null1) = colnames(deng.null2) = c('year','week','frcstYear','frcstWeek','nhorizon',state.df$name)

dimnames(deng.frcst)[[2]] = dimnames(deng.err)[[2]] = dimnames(deng.rel)[[2]] = c('year','week','frcstYear','frcstWeek','nhorizon',state.df$name)

arima.names = array(NA,narima)
iarm = 1
for (ip in 1:length(p)) {
	for (id in 1:length(d)) {
		for (iq in 1:length(q)) {
			for (iP in 1:length(P)) {
				for (iD in 1:length(D)) {
					for (iQ in 1:length(Q)) {
						arima.names[iarm] = paste0('(',p[ip],d[id],q[iq],')(',P[iP],D[iD],Q[iQ],')')
						iarm = iarm + 1
					}
				}
			}
		}
	}
}

dimnames(deng.frcst)[[3]]  = dimnames(deng.err)[[3]]  = dimnames(deng.rel)[[3]]  = arima.names

icount = 0

for (iyear in 1:length(my.year.vec)) {


	my.year = my.year.vec[iyear]
	cat('Processing Year: ',my.year,'\n')
	year_ave = unique(years)[iyear:(iyear + nave - 1)]
	
	index = which(years >= min(year_ave) & years <= max(year_ave))

	trnData = stateData[index]

	trnYear = years[index]

	trnWeek = weeks[index]

	deng.year.ave[1:week_per_year, 1] = round(mean(trnData,na.rm=TRUE))

	for (j in 1:length(week_ave)) {
		my.week = week_ave[j]
		jndex = as.numeric(which(trnWeek == my.week))
		deng.week.ave[j, 1] = round(mean(trnData[jndex], na.rm = TRUE))
	}

	for (j in 1:week_per_year) {

		# This is the month that the forecast is made - last month of real data
		my.week = week_ave[j]
		my.row = which(years == my.year & weeks == my.week)
		if (length(my.row) == 0) next
		cat('Processing Week: ', my.week,'\n')

		## Now do an ARIMF fit on the past ten year starting from now
		istart = my.row - 10 * week_per_year + 1
		iend   = my.row
##
## Do ARIMA fits and make a forecast nfrcst forward
##
		jstart = my.row+1
		jend = my.row + nhorizon
		jstart = min(noTPts, jstart)
		jend = min(noTPts, jend)

		## current data 
	    jnow = my.row
	    cases.now = stateData[jnow]

		frcstWeek = frcst.table[my.week,2:(nhorizon+1)]
		frcstYear = year[jstart:jend]

		nfrcst = length(frcstYear)
		frcstWeek = as.numeric(frcstWeek[1:nfrcst])

		week_an = rep(weeks[my.row], nfrcst)
		year_an = rep(years[my.row], nfrcst)
		horizon = 1:nfrcst

		iarm = 0
		for (ip in 1:length(p)) {

			for (id in 1:length(d)) {

				for (iq in 1:length(q)) {

					for (iP in 1:length(P)) {

						for (iD in 1:length(D)) {

							# increase ARIMA model number counter
							iarm = iarm + 1
							## Let's fit the log of the data with a mean of zero and sd of one

								x = stateData[istart:iend]
								x[is.na(x)] <- 0
								x = log(x + 1)
								mean.x = mean(x)
								x = x - mean.x
								sd.x = sd(x)
								x = x/sd.x
								fit <- try(arima(x, order = c(p[ip], d[id], q[iq]), seasonal = list(order = c(P[iP], D[iD], Q)), method = "ML"))
								if (isTRUE(class(fit) == "try-error"))
									next
								result <- try(predict(fit, n.ahead = nfrcst))
								if (isTRUE(class(result) == "try-error"))
									next
								if (length(result$pred) == length((icount + 1):(icount + nfrcst))) {
									## Need to convert back from log to linear
									my.result = result$pred[1:nfrcst]
									my.result = exp(my.result * sd.x + mean.x)
									deng.frcst[(icount + 1):(icount + nfrcst), state.df$name, iarm] = my.result

								}
						}
					}
				}
			}
		} ## End of loops over ARIMA models

		deng.frcst[(icount + 1):(icount + nfrcst), "year",1:narima]  = year_an
		deng.obsrv[(icount + 1):(icount + nfrcst), "year"] = year_an
		deng.err[(icount + 1):(icount + nfrcst), "year",  1:narima]  = year_an
		deng.rel[(icount + 1):(icount + nfrcst), "year",  1:narima]  = year_an
		deng.null1[(icount + 1):(icount + nfrcst), "year"]  = year_an
		deng.null2[(icount + 1):(icount + nfrcst), "year"] = year_an

		deng.frcst[(icount + 1):(icount + nfrcst), "week",1:narima]  = week_an
		deng.obsrv[(icount + 1):(icount + nfrcst), "week"]  = week_an
		deng.err[(icount + 1):(icount + nfrcst),   "week",1:narima]  = week_an
		deng.rel[(icount + 1):(icount + nfrcst),  "week", 1:narima]  = week_an
		deng.null1[(icount + 1):(icount + nfrcst), "week"]  = week_an
		deng.null2[(icount + 1):(icount + nfrcst), "week"]  = week_an

		deng.frcst[(icount + 1):(icount + nfrcst), "frcstYear",1:narima]  = frcstYear
		deng.obsrv[(icount + 1):(icount + nfrcst), "frcstYear"]  = frcstYear
		deng.err[(icount + 1):(icount + nfrcst), "frcstYear",  1:narima]  = frcstYear
		deng.rel[(icount + 1):(icount + nfrcst), "frcstYear", 1:narima ] = frcstYear
		deng.null1[(icount + 1):(icount + nfrcst), "frcstYear"]  = frcstYear
		deng.null2[(icount + 1):(icount + nfrcst), "frcstYear"] = frcstYear

		deng.frcst[(icount + 1):(icount + nfrcst), "frcstWeek",1:narima] = frcstWeek
		deng.obsrv[(icount + 1):(icount + nfrcst), "frcstWeek"] = frcstWeek
		deng.err[(icount + 1):(icount + nfrcst), "frcstWeek",  1:narima] = frcstWeek
		deng.rel[(icount + 1):(icount + nfrcst), "frcstWeek",  1:narima] = frcstWeek
		deng.null1[(icount + 1):(icount + nfrcst), "frcstWeek"] = frcstWeek
		deng.null2[(icount + 1):(icount + nfrcst), "frcstWeek"] = frcstWeek

		deng.frcst[(icount + 1):(icount + nfrcst), "nhorizon",1:narima] = 1:nfrcst
		deng.obsrv[(icount + 1):(icount + nfrcst), "nhorizon"] = 1:nfrcst
		deng.err[(icount + 1):(icount + nfrcst), "nhorizon",   1:narima] = 1:nfrcst
		deng.rel[(icount + 1):(icount + nfrcst), "nhorizon",   1:narima] = 1:nfrcst
		deng.null1[(icount + 1):(icount + nfrcst), "nhorizon"] = 1:nfrcst
		deng.null2[(icount + 1):(icount + nfrcst), "nhorizon"]  = 1:nfrcst

		mat = stateData[jstart:jend]
		mat[mat < 1] <- 1

		deng.obsrv[(icount + 1):(icount + nfrcst), state.df$name] = mat

		deng.err[(icount + 1):(icount + nfrcst), state.df$name,1:narima] = abs(as.numeric(deng.frcst[(icount + 1):(icount + nfrcst), state.df$name,1:narima]) - as.numeric(deng.obsrv[(icount +
			1):(icount + nfrcst), state.df$name]))

		deng.rel[(icount + 1):(icount + nfrcst), state.df$name,1:narima] = as.numeric(deng.err[(icount + 1):(icount + nfrcst), state.df$name,1:narima])/as.numeric(deng.obsrv[(icount + 1):(icount +
			nfrcst), state.df$name])


		deng.null1[(icount + 1):(icount + nfrcst), state.df$name] = abs(as.numeric(deng.year.ave[frcstWeek, ]) - stateData[jstart:jend])
		deng.null2[(icount + 1):(icount + nfrcst), state.df$name] = abs(as.numeric(deng.week.ave[frcstWeek, ]) - stateData[jstart:jend])

		icount = icount + nfrcst

	}

}

deng.err = deng.err[1:icount,,]
deng.rel = deng.rel[1:icount,,]
deng.frcst = deng.frcst[1:icount,,]
deng.obsrv = deng.obsrv[1:icount,]
deng.null1 = deng.null1[1:icount,]
deng.null2 = deng.null2[1:icount,]

## Now if we want to plot or average we can
## we can calculate the mean absolute error for predicting 1, 2, 3 and 4 month ahead

frcst.mae = array(NA, c(nhorizon, narima))
frcst.mre = array(NA, c(nhorizon, narima))
frcst.rel1 = array(NA, c(nhorizon, narima))
frcst.rel2 = array(NA, c(nhorizon, narima))

colnames(frcst.mae) = colnames(frcst.mre) = colnames(frcst.rel1) = colnames(frcst.rel2) = arima.names

for (i in 1:nhorizon) {

	index = which(as.numeric(deng.err[, "nhorizon", 1]) == i)

	my.state = state.df$name
	subset = deng.err[index, my.state, 1:narima]

	for (iarm in 1:narima) {
		frcst.mae[i, iarm] = mean(as.numeric(subset[, iarm]), na.rm = TRUE)
	}
	subset = deng.rel[index, my.state, 1:narima]

	for (iarm in 1:narima) {
		frcst.mre[i, iarm] = mean(subset[, iarm], na.rm = TRUE)
		if (is.nan(frcst.mre[i, iarm])) 
			frcst.mre[i, iarm] = NA
	}
	## Using NULL1
	subset = deng.null1[index, my.state]
	frcst.rel1[i, 1:narima] = mean(subset, na.rm = TRUE)
	frcst.rel1[i, 1:narima] = frcst.mae[i, 1:narima]/frcst.rel1[i, 1:narima]

	## Using NULL2
	subset = deng.null2[index, my.state]
	frcst.rel2[i, 1:narima] = mean(subset, na.rm = TRUE)
	frcst.rel2[i, 1:narima] = frcst.mae[i, 1:narima]/frcst.rel2[i, 1:narima]

}
print(frcst.rel2)
## write all the results - plotting will be done seperately because the code takes too long to run
results = list(frcst.mae=frcst.mae, frcst.mre = frcst.mre, frcst.rel1 = frcst.rel1, frcst.rel2 = frcst.rel2)

filename = paste0(mod_name[[1]],'.arima.RData')
save(results, file = filename)

###






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

mod_name = c(NAME_2="MX")

fit_names = 'all'

mod_level = 2

fit_level = 3

nfit = nperiodsFit = 52

year = 2017

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

## rename columns using the full state names
colnames(stateData) = state_names

## convert to a matrix
stateData = as.matrix(stateData)

state.df = data.frame(name=rep(0,nstates),lon=rep(0,nstates),lat=rep(0,nstates),pop=rep(0,nstates))

state.df$name = state_names
state.df$lon  = latlon$lon
state.df$lat  = latlon$lat
state.df$pop  = as.numeric(pop)

## Reorder the state.df and state Data by latitude - in increasing order

iorder = order(state.df$lat)
state.df = state.df[iorder,]

## Reorder the stateData
mydf = stateData
for (i in 1:nstates) {
	stateName = colnames(stateData)[i]
	#stateName = substr(stateName,4,100)
	index = which(state.df[,'name'] == stateName)
	mydf[,index] = stateData[,i]
	colnames(mydf)[index] = stateName
}

stateData = mydf



if (mod_name[[1]] == "MX") {
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

deng.month.ave = array(0,c(length(month_ave),nstates))
deng.year.ave  = array(0,c(length(month_ave),nstates))


colnames(deng.month.ave) = colnames(deng.year.ave) = state.df$name
rownames(deng.month.ave) = rownames(deng.year.ave) = month_ave

nhorizon = 4

##
## Define all the ARIMA models we want to test
##


p = 1
d = 0
q = 0
P = 0:3
D = 0:1
Q = 0

Period = 12

narima = length(p) * length(d) * length(q) * length(P) * length(D) * length(Q)
month_per_year = 12

deng.frcst = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nstates),narima))
deng.obsrv = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nstates)))
deng.err   = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nstates),narima))
deng.rel   = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nstates),narima))
deng.null1 = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nstates)))
deng.null2 = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nstates)))

##
## The national is direct and aggregate, more can be added

natl.names = c('Direct',"Aggregate")
nnatl = length(natl.names)


ndeng.frcst = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nnatl),narima))
ndeng.obsrv = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nnatl)))
ndeng.err   = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nnatl),narima))
ndeng.rel   = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nnatl),narima))
ndeng.null1 = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nnatl)))
ndeng.null2 = array(NA, c((nhorizon*month_per_year*length(my.year.vec)),(5+nnatl)))

frcst.table = data.frame(current=month_ave,future=array(NA,c(length(month_ave),nhorizon)))
frcst.table[,2] = c(2:12,1)
frcst.table[,3] = c(3:12,1,2)
frcst.table[,4] = c(4:12,1,2,3)
frcst.table[,5] = c(5:12,1,2,3,4)

colnames(deng.obsrv) = colnames(deng.null1) = colnames(deng.null2) = c('year','month','frcstYear','frcstMonth','nhorizon',state.df$name)

colnames(ndeng.obsrv) = colnames(ndeng.null1) = colnames(ndeng.null2) = c('year','month','frcstYear','frcstMonth','nhorizon', natl.names)

dimnames(deng.frcst)[[2]] = dimnames(deng.err)[[2]] = dimnames(deng.rel)[[2]] = c('year','month','frcstYear','frcstMonth','nhorizon',state.df$name)

dimnames(ndeng.frcst)[[2]] = dimnames(ndeng.err)[[2]] = dimnames(ndeng.rel)[[2]] = c('year','month','frcstYear','frcstMonth','nhorizon', natl.names)

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
dimnames(ndeng.frcst)[[3]] = dimnames(ndeng.err)[[3]] = dimnames(ndeng.rel)[[3]] = arima.names

icount = 0

for (iyear in 1:length(my.year.vec)) {


	my.year = my.year.vec[iyear]
	cat('Processing Year: ',my.year,'\n')
	year_ave = unique(years)[iyear:(iyear + nave - 1)]

	index = which(years >= min(year_ave) & years <= max(year_ave))

	trnData = stateData[index, ]

	trnYear = years[index]

	trnMonth = months[index]

	for (i in 1:nstates) {
		deng.year.ave[1:month_per_year, i] = round(mean(trnData[, i],na.rm=TRUE))
	}

	for (j in 1:length(month_ave)) {
		my.month = month_ave[j]
		jndex = as.numeric(which(trnMonth == my.month))
		for (i in 1:nstates) {
			deng.month.ave[j, i] = round(mean(trnData[jndex, i], na.rm = TRUE))
		}

	}

	for (j in 1:month_per_year) {

		# This is the month that the forecast is made - last month of real data
		my.month = month_ave[j]
		my.row = which(years == my.year & months == my.month)

		cat('Processing Month: ', my.month,'\n')

		## Now do an ARIMF fit on the past ten year starting from now
		istart = my.row - 10 * month_per_year + 1
		iend   = my.row
##
## Do ARIMA fits and make a forecast nfrcst forward
##
		jstart = my.row+1
		jend = my.row + 4
		jstart = min(noTPts, jstart)
		jend = min(noTPts, jend)

		## current data for all states
	    jnow = my.row
	    cases.now = stateData[jnow,]

		frcstMonth = frcst.table[my.month,2:5]
		frcstYear = year[jstart:jend]

		nfrcst = length(frcstYear)
		frcstMonth = as.numeric(frcstMonth[1:nfrcst])

		month_an = rep(months[my.row], nfrcst)
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
								result <- try(predict(fit, n.ahead = nfrcst))
								if (isTRUE(class(result) == "try-error"))
									next
								if (length(result$pred) == length((icount + 1):(icount + nfrcst))) {
									## Need to convert back from log to linear
									my.result = result$pred[1:nfrcst]
									my.result = exp(my.result * sd.x + mean.x)
									deng.frcst[(icount + 1):(icount + nfrcst), state.df$name[istate], iarm] = my.result

								}

							} # End of loop over states

							## Aggregate the state level data to National

							for (kcount in (icount + 1):(icount + nfrcst)){
								ndeng.frcst[kcount, "Aggregate", iarm] = sum(deng.frcst[kcount, (6:(5+nstates)), iarm],na.rm=TRUE)
							}

							## This fails when nfrcst = 1

							##ndeng.frcst[(icount + 1):(icount + nfrcst), "Aggregate", iarm] = rowSums(deng.frcst[(icount + 1):(icount + nfrcst), (6:(5+nstates)), iarm],na.rm=TRUE)

							## direct national

							 x = rowSums(stateData[istart:iend, ],na.rm=TRUE)
							 x[is.na(x)] <- 0
							 x = log(x + 1)
							 mean.x = mean(x)
							 x = x - mean.x
							 sd.x = sd(x)
							 x = x/sd.x
							 fit <- try(arima(x, order = c(p[ip], d[id], q[iq]), seasonal = list(order = c(P[iP], D[iD], Q), period = 12), method = "ML"))
							 if (isTRUE(class(fit) == "try-error"))
								next
							 result <- try(predict(fit, n.ahead = nfrcst, se.fit = TRUE))
							 if (isTRUE(class(result) == "try-error"))
									 next
							 if (length(result$pred) == length((icount + 1):(icount + nfrcst))) {
								## Need to convert back from log to linear
								my.result = result$pred[1:nfrcst]
								my.result = exp(my.result * sd.x + mean.x)
								ndeng.frcst[(icount + 1):(icount + nfrcst), 'Direct', iarm] = my.result

							}
						}
					}
				}
			}
		} ## End of loops over ARIMA models

		deng.frcst[(icount + 1):(icount + nfrcst), "year",1:narima] = ndeng.frcst[(icount + 1):(icount + nfrcst), "year",1:narima] = year_an
		deng.obsrv[(icount + 1):(icount + nfrcst), "year"] = ndeng.obsrv[(icount + 1):(icount + nfrcst), "year"] = year_an
		deng.err[(icount + 1):(icount + nfrcst), "year",  1:narima] = ndeng.err[(icount + 1):(icount + nfrcst), "year",  1:narima] = year_an
		deng.rel[(icount + 1):(icount + nfrcst), "year",  1:narima] = ndeng.rel[(icount + 1):(icount + nfrcst), "year",  1:narima] = year_an
		deng.null1[(icount + 1):(icount + nfrcst), "year"] = ndeng.null1[(icount + 1):(icount + nfrcst), "year"] = year_an
		deng.null2[(icount + 1):(icount + nfrcst), "year"] = ndeng.null2[(icount + 1):(icount + nfrcst), "year"] = year_an

		deng.frcst[(icount + 1):(icount + nfrcst), "month",1:narima] = ndeng.frcst[(icount + 1):(icount + nfrcst), "month",1:narima] = month_an
		deng.obsrv[(icount + 1):(icount + nfrcst), "month"] = ndeng.obsrv[(icount + 1):(icount + nfrcst), "month"] = month_an
		deng.err[(icount + 1):(icount + nfrcst),   "month",1:narima] = ndeng.err[(icount + 1):(icount + nfrcst),   "month",1:narima] = month_an
		deng.rel[(icount + 1):(icount + nfrcst),  "month", 1:narima] = ndeng.rel[(icount + 1):(icount + nfrcst),  "month", 1:narima] = month_an
		deng.null1[(icount + 1):(icount + nfrcst), "month"] = ndeng.null1[(icount + 1):(icount + nfrcst), "month"] = month_an
		deng.null2[(icount + 1):(icount + nfrcst), "month"] = ndeng.null2[(icount + 1):(icount + nfrcst), "month"] = month_an

		deng.frcst[(icount + 1):(icount + nfrcst), "frcstYear",1:narima] = ndeng.frcst[(icount + 1):(icount + nfrcst), "frcstYear",1:narima] = frcstYear
		deng.obsrv[(icount + 1):(icount + nfrcst), "frcstYear"] = ndeng.obsrv[(icount + 1):(icount + nfrcst), "frcstYear"] = frcstYear
		deng.err[(icount + 1):(icount + nfrcst), "frcstYear",  1:narima] = ndeng.err[(icount + 1):(icount + nfrcst), "frcstYear",  1:narima] = frcstYear
		deng.rel[(icount + 1):(icount + nfrcst), "frcstYear", 1:narima ] = ndeng.rel[(icount + 1):(icount + nfrcst), "frcstYear", 1:narima ] = frcstYear
		deng.null1[(icount + 1):(icount + nfrcst), "frcstYear"] = ndeng.null1[(icount + 1):(icount + nfrcst), "frcstYear"] = frcstYear
		deng.null2[(icount + 1):(icount + nfrcst), "frcstYear"] = ndeng.null2[(icount + 1):(icount + nfrcst), "frcstYear"] = frcstYear

		deng.frcst[(icount + 1):(icount + nfrcst), "frcstMonth",1:narima] = ndeng.frcst[(icount + 1):(icount + nfrcst), "frcstMonth",1:narima] = frcstMonth
		deng.obsrv[(icount + 1):(icount + nfrcst), "frcstMonth"] = ndeng.obsrv[(icount + 1):(icount + nfrcst), "frcstMonth"] = frcstMonth
		deng.err[(icount + 1):(icount + nfrcst), "frcstMonth",  1:narima] = ndeng.err[(icount + 1):(icount + nfrcst), "frcstMonth",  1:narima] = frcstMonth
		deng.rel[(icount + 1):(icount + nfrcst), "frcstMonth",  1:narima] = ndeng.rel[(icount + 1):(icount + nfrcst), "frcstMonth",  1:narima] = frcstMonth
		deng.null1[(icount + 1):(icount + nfrcst), "frcstMonth"] = ndeng.null1[(icount + 1):(icount + nfrcst), "frcstMonth"] = frcstMonth
		deng.null2[(icount + 1):(icount + nfrcst), "frcstMonth"] = ndeng.null2[(icount + 1):(icount + nfrcst), "frcstMonth"] = frcstMonth

		deng.frcst[(icount + 1):(icount + nfrcst), "nhorizon",1:narima] = ndeng.frcst[(icount + 1):(icount + nfrcst), "nhorizon",1:narima] = 1:nfrcst
		deng.obsrv[(icount + 1):(icount + nfrcst), "nhorizon"] = ndeng.obsrv[(icount + 1):(icount + nfrcst), "nhorizon"] = 1:nfrcst
		deng.err[(icount + 1):(icount + nfrcst), "nhorizon",   1:narima] = ndeng.err[(icount + 1):(icount + nfrcst), "nhorizon",   1:narima] = 1:nfrcst
		deng.rel[(icount + 1):(icount + nfrcst), "nhorizon",   1:narima] = ndeng.rel[(icount + 1):(icount + nfrcst), "nhorizon",   1:narima] = 1:nfrcst
		deng.null1[(icount + 1):(icount + nfrcst), "nhorizon"] = ndeng.null1[(icount + 1):(icount + nfrcst), "nhorizon"] = 1:nfrcst
		deng.null2[(icount + 1):(icount + nfrcst), "nhorizon"] = ndeng.null2[(icount + 1):(icount + nfrcst), "nhorizon"] = 1:nfrcst

		mat = stateData[jstart:jend, ]
		mat[mat < 1] <- 1

		deng.obsrv[(icount + 1):(icount + nfrcst), state.df$name] = mat

		for (kcount in (icount + 1):(icount + nfrcst)){
			ndeng.obsrv[kcount, "Direct"] = ndeng.obsrv[kcount, "Aggregate"] = sum(deng.obsrv[kcount,state.df$name],na.rm=TRUE) ##rowSums(mat,na.rm=TRUE)
		}

		deng.err[(icount + 1):(icount + nfrcst), state.df$name,1:narima] = abs(as.numeric(deng.frcst[(icount + 1):(icount + nfrcst), state.df$name,1:narima]) - as.numeric(deng.obsrv[(icount +
			1):(icount + nfrcst), state.df$name]))

		ndeng.err[(icount + 1):(icount + nfrcst), 'Direct',1:narima] = abs(as.numeric(ndeng.frcst[(icount + 1):(icount + nfrcst), "Direct",1:narima]) - as.numeric(ndeng.obsrv[(icount +
			1):(icount + nfrcst), "Direct"]))

		ndeng.err[(icount + 1):(icount + nfrcst), 'Aggregate',1:narima] = abs(as.numeric(ndeng.frcst[(icount + 1):(icount + nfrcst), "Aggregate",1:narima]) - as.numeric(ndeng.obsrv[(icount +
			1):(icount + nfrcst), "Aggregate"]))



		deng.rel[(icount + 1):(icount + nfrcst), state.df$name,1:narima] = as.numeric(deng.err[(icount + 1):(icount + nfrcst), state.df$name,1:narima])/as.numeric(deng.obsrv[(icount + 1):(icount +
			nfrcst), state.df$name])

		ndeng.rel[(icount + 1):(icount + nfrcst), "Direct",1:narima] = as.numeric(ndeng.err[(icount + 1):(icount + nfrcst), "Direct",1:narima])/as.numeric(ndeng.obsrv[(icount + 1):(icount +
			nfrcst), "Direct"])

		ndeng.rel[(icount + 1):(icount + nfrcst), "Aggregate",1:narima] = as.numeric(ndeng.err[(icount + 1):(icount + nfrcst), "Aggregate",1:narima])/as.numeric(ndeng.obsrv[(icount + 1):(icount +
			nfrcst), "Aggregate"])

		deng.null1[(icount + 1):(icount + nfrcst), state.df$name] = abs(as.numeric(deng.year.ave[frcstMonth, ] ) - stateData[jstart:jend, ])
		deng.null2[(icount + 1):(icount + nfrcst), state.df$name] = abs(as.numeric(deng.month.ave[frcstMonth, ]) - stateData[jstart:jend, ])

		if (jstart != jend) {

		ndeng.null1[(icount + 1):(icount + nfrcst), "Direct"] = abs(as.numeric(rowSums(deng.year.ave[frcstMonth, ] ,na.rm=TRUE)) - rowSums(stateData[jstart:jend, ],na.rm=TRUE))
		ndeng.null2[(icount + 1):(icount + nfrcst), "Direct"] = abs(as.numeric(rowSums(deng.month.ave[frcstMonth, ],na.rm=TRUE)) - rowSums(stateData[jstart:jend, ],na.rm=TRUE))

		ndeng.null1[(icount + 1):(icount + nfrcst), "Aggregate"] = abs(as.numeric(rowSums(deng.year.ave[frcstMonth, ] ,na.rm=TRUE)) - rowSums(stateData[jstart:jend, ],na.rm=TRUE))
		ndeng.null2[(icount + 1):(icount + nfrcst), "Aggregate"] = abs(as.numeric(rowSums(deng.month.ave[frcstMonth, ],na.rm=TRUE)) - rowSums(stateData[jstart:jend, ],na.rm=TRUE))
		} else {
		ndeng.null1[(icount + 1):(icount + nfrcst), "Direct"] = abs(as.numeric(sum(deng.year.ave[frcstMonth, ] ,na.rm=TRUE)) - sum(stateData[jstart, ],na.rm=TRUE))
		ndeng.null2[(icount + 1):(icount + nfrcst), "Direct"] = abs(as.numeric(sum(deng.month.ave[frcstMonth, ],na.rm=TRUE)) - sum(stateData[jstart, ],na.rm=TRUE))

		ndeng.null1[(icount + 1):(icount + nfrcst), "Aggregate"] = abs(as.numeric(sum(deng.year.ave[frcstMonth, ] ,na.rm=TRUE)) - sum(stateData[jstart, ],na.rm=TRUE))
		ndeng.null2[(icount + 1):(icount + nfrcst), "Aggregate"] = abs(as.numeric(sum(deng.month.ave[frcstMonth, ],na.rm=TRUE)) - sum(stateData[jstart, ],na.rm=TRUE))
		}

		icount = icount + nfrcst

	}

}

deng.err = deng.err[1:icount,,]
deng.rel = deng.rel[1:icount,,]
deng.frcst = deng.frcst[1:icount,,]
deng.obsrv = deng.obsrv[1:icount,]
deng.null1 = deng.null1[1:icount,]
deng.null2 = deng.null2[1:icount,]

ndeng.err = ndeng.err[1:icount,,]
ndeng.rel = ndeng.rel[1:icount,,]
ndeng.frcst = ndeng.frcst[1:icount,,]
ndeng.obsrv = ndeng.obsrv[1:icount,]
ndeng.null1 = ndeng.null1[1:icount,]
ndeng.null2 = ndeng.null2[1:icount,]


## Now if we want to plot or average we can
## we can calculate the mean absolute error for predicting 1, 2, 3 and 4 month ahead

frcst.mae = array(NA, c(nhorizon, nstates, narima))
frcst.mre = array(NA, c(nhorizon, nstates, narima))
frcst.rel1 = array(NA, c(nhorizon, nstates, narima))
frcst.rel2 = array(NA, c(nhorizon, nstates, narima))

dimnames(frcst.mae)[[2]] = dimnames(frcst.mre)[[2]] = dimnames(frcst.rel1)[[2]] = dimnames(frcst.rel2)[[2]] = state.df$name
dimnames(frcst.mae)[[3]] = dimnames(frcst.mre)[[3]] = dimnames(frcst.rel1)[[3]] = dimnames(frcst.rel2)[[3]] = arima.names

## Repeat on the National level

nfrcst.mae = array(NA, c(nhorizon, nnatl, narima))
nfrcst.mre = array(NA, c(nhorizon, nnatl, narima))
nfrcst.rel1 = array(NA, c(nhorizon, nnatl, narima))
nfrcst.rel2 = array(NA, c(nhorizon, nnatl, narima))

dimnames(nfrcst.mae)[[2]] = dimnames(nfrcst.mre)[[2]] = dimnames(nfrcst.rel1)[[2]] = dimnames(nfrcst.rel2)[[2]] = c("Direct","Aggregate")
dimnames(nfrcst.mae)[[3]] = dimnames(nfrcst.mre)[[3]] = dimnames(nfrcst.rel1)[[3]] = dimnames(nfrcst.rel2)[[3]] = arima.names

for (i in 1:nhorizon) {

	index = which(as.numeric(deng.err[, "nhorizon", 1]) == i)

	for (istate in 1:nstates) {

		my.state = state.df$name[istate]
		subset = deng.err[index, my.state, 1:narima]

		for (iarm in 1:narima) {
			frcst.mae[i, istate, iarm] = mean(as.numeric(subset[, iarm]), na.rm = TRUE)
		}
		subset = deng.rel[index, my.state, 1:narima]

		for (iarm in 1:narima) {
			frcst.mre[i, istate, iarm] = mean(subset[, iarm], na.rm = TRUE)
			if (is.nan(frcst.mre[i, istate, iarm]))
				frcst.mre[i, istate, iarm] = NA
		}
		## Using NULL1
		subset = deng.null1[index, my.state]
		frcst.rel1[i, istate, 1:narima] = mean(subset, na.rm = TRUE)
		frcst.rel1[i, istate, 1:narima] = frcst.mae[i, istate, 1:narima]/frcst.rel1[i, istate,1:narima]

		## Using NULL2
		subset = deng.null2[index, my.state]
		frcst.rel2[i, istate, 1:narima] = mean(subset, na.rm = TRUE)
		frcst.rel2[i, istate, 1:narima] = frcst.mae[i, istate, 1:narima]/frcst.rel2[i, istate,1:narima]

	}

	# repeat for the National direct and aggregate

	index = which(as.numeric(ndeng.err[, "nhorizon", 1]) == i)

	for (istate in 1:nnatl) {
		my.state = natl.names[istate]
		subset = ndeng.err[index, my.state, 1:narima]
		for (iarm in 1:narima) {
			nfrcst.mae[i, my.state, iarm] = mean(as.numeric(subset[, iarm]), na.rm = TRUE)
		}

		subset = ndeng.rel[index, my.state, 1:narima]

		for (iarm in 1:narima) {
			nfrcst.mre[i, my.state, iarm] = mean(subset[, iarm], na.rm = TRUE)
			if (is.nan(nfrcst.mre[i, my.state, iarm]))
				nfrcst.mre[i,my.state, iarm] = NA
		}

		## Using NULL1
		subset = ndeng.null1[index, my.state]
		nfrcst.rel1[i, my.state, 1:narima] = mean(subset, na.rm = TRUE)
		nfrcst.rel1[i, my.state, 1:narima] = nfrcst.mae[i, my.state, 1:narima]/nfrcst.rel1[i,my.state,1:narima]
		## Using NULL2
		subset = ndeng.null2[index, my.state]
		nfrcst.rel2[i, my.state, 1:narima] = mean(subset, na.rm = TRUE)
		nfrcst.rel2[i, my.state, 1:narima] = nfrcst.mae[i, my.state, 1:narima]/nfrcst.rel2[i, my.state,1:narima]


	}

}

## Now lets order the states by population - using the last year
iorder = order(state.df$pop,decreasing=FALSE)

frcst.mae  = frcst.mae[,iorder,]
frcst.mre  = frcst.mre[,iorder,]
frcst.rel1 = frcst.rel1[,iorder,]
frcst.rel2 = frcst.rel2[,iorder,]

## write all the results - plotting will be done seperately because the code takes too long to run
results = list(frcst.mae=frcst.mae, frcst.mre = frcst.mre, frcst.rel1 = frcst.rel1, frcst.rel2 = frcst.rel2,nfrcst.mae=nfrcst.mae, nfrcst.mre = nfrcst.mre, nfrcst.rel1 = nfrcst.rel1, nfrcst.rel2 = nfrcst.rel2)

filename = paste0(mod_name[[1]],'.arima.RData')
save(results, file = filename)

###






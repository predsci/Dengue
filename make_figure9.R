##
## CalCulation for Figure 9 - effect of delay in reporting Singapore weekly data
##

rm(list = ls())


library("Hmisc")
library("zoo")
library("RColorBrewer")

horizon = 1:10

nhorizon = length(horizon)

require(DICE)

##
## select country using ISO-2 code
##

pdf = TRUE

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

year.max = max(years)
year.min = year.max-2

## make a nice table
sg <- cbind(years,weeks,stateData)

##
## we will define three times: to-date = t_to, delivery-date = t_del and analysis date = t_an
## the t_to specifies that the current forecast will only use cases whose sympotm onset date is equal or less than t_to
## the delivery date specifies that the current forecast will use data up to and including t_del
## the analysis date specified when a given forecast was run
## for now lets set set the t_to to be the same as the delivery date 
## There are 13 years of data so let's use the previous ten to calculate averages 
## mean value for every week
##

nave = 10
nave1 = nave + 1

my.year.vec = unique(years)[nave1:nyears]

week_ave = unique(weeks)

##
## Define all the ARIMA models we want to test
##

p = 1
d = 0
q = 0
P = 3
D = 0
Q = 0

Period = 52

narima = length(p) * length(d) * length(q) * length(P) * length(D) * length(Q)

deng.obsrv = deng.null = array(NA, c((week_per_year * length(year.min:year.max)), 3))
deng.upper = array(NA, c((week_per_year * length(year.min:year.max)), 2 + nhorizon))
deng.lower = array(NA, c((week_per_year * length(year.min:year.max)), 2 + nhorizon))
deng.frcst = array(NA, c((week_per_year * length(year.min:year.max)), 2 + nhorizon))

colnames(deng.frcst) = colnames(deng.lower) = colnames(deng.upper) = c("year", "week", paste0("horizon", 1:nhorizon))
colnames(deng.obsrv) = colnames(deng.null) = c("year", "week", "cases")

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



for (ihrz in 1:length(horizon)) {

	my.horizon = horizon[ihrz]

	cat("\n Horizon Prediction: ", my.horizon, "\n\n")

	icount = 1

	for (my.year in year.min:year.max) {

		#my.year = my.year.vec[iyear]
		cat("Processing Year: ", my.year, "\n")

		for (j in 1:week_per_year) {
			## 
			if (my.year == year.min && j < 17) 
				next
				
			# This is the week that the forecast is made - last week of real data
			my.week = week_ave[j]
			my.row = which(years == my.year & weeks == my.week)

			if (length(my.row) == 0) 
				next

			# let's calcualte a NULL model for this week
			idx = which(weeks == my.week)
			## Find which idx == my.row and take just the ten before
			ilast = which(idx == my.row) 
			istart = ilast - 10
			iend   = ilast -1 
			idx = idx[istart:iend]
			
			cat("Processing Week: ", my.week, "\n")

			## Now do an ARIMF fit on the past ten year starting from now
			istart = my.row - 10 * week_per_year - my.horizon + 1
			istart = max(istart, 1)
			iend = my.row - my.horizon
			##
			## Do ARIMA fits and make a forecast nfrcst forward
##		

			## current data for all states - this is the observed number of cases
			
			cases.now = stateData[my.row]


			## This is when the forecast is made
			frcstWeek = weeks[(my.row - my.horizon)]
			frcstYear = years[(my.row - my.horizon)]
			cat("forecast made on", frcstWeek, frcstYear, "\n")

			iarm = 1
			for (ip in 1:length(p)) {

				for (id in 1:length(d)) {

					for (iq in 1:length(q)) {

						for (iP in 1:length(P)) {

							for (iD in 1:length(D)) {

								## Let's fit the log of the data with a mean of zero and sd of one
								
								x = stateData[istart:iend]
								x[is.na(x)] <- 0
								x = log(x + 1)
								mean.x = mean(x)
								x = x - mean.x
								sd.x = sd(x)
								x = x/sd.x
								fit <- try(arima(x, order = c(p[ip], d[id], q[iq]), seasonal = list(order = c(P[iP], D[iD], Q), period = 12), method = "ML"))
								if (isTRUE(class(fit) == "try-error")) {## Give it nother try with a different arima order
									fit <- try(arima(x, order = c(p[ip], d[id], q[iq]), seasonal = list(order = c(3, 1, 0), period = 12), method = "ML"))
									
								}
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

									deng.frcst[icount, (2 + ihrz)] = my.result
									deng.upper[icount, (2 + ihrz)] = exp(lower)
									deng.lower[icount, (2 + ihrz)] = exp(upper)

								}

								iarm = iarm + 1
							}
						}
					}
				}
			}

			deng.obsrv[icount, "cases"] = cases.now
			deng.null[ icount,"cases" ] = mean(stateData[idx])
			
			deng.frcst[icount, "year"] = my.year
			deng.upper[icount, "year"] = my.year
			deng.lower[icount, "year"] = my.year
			deng.obsrv[icount, "year"] = my.year
			deng.null[ icount, "year"] = my.year 
			
			deng.frcst[icount, "week"] = my.week
			deng.upper[icount, "week"] = my.week
			deng.lower[icount, "week"] = my.week
			deng.obsrv[icount, "week"] = my.week
			deng.null[ icount, "week"] = my.week 
			
			icount = icount + 1

		}

	}

}

## Trim
deng.null = deng.null[1:(icount-1),]
deng.obsrv = deng.obsrv[1:(icount-1),]
deng.lower = deng.lower[1:(icount-1),]
deng.upper = deng.upper[1:(icount-1),]
deng.frcst = deng.frcst[1:(icount-1),]


## Error analysis

null.mae  = deng.obsrv

frcst.mae = frcst.rel = deng.frcst

colnames(null.mae)  = c('year','week', 'ae')
colnames(frcst.mae) = c('year','week',paste0('horizon',1:nhorizon))
colnames(frcst.rel) = c('year','week',paste0('horizon',1:nhorizon))

null.mae[,'ae']  = abs(deng.null[,'cases']- deng.obsrv[,'cases'])

for (ihrz in 1:nhorizon) {
	frcst.mae[,(2+ihrz)] = abs(deng.frcst[,(2+ihrz)]-deng.obsrv[,'cases'])
	frcst.rel[,(2+ihrz)] = frcst.mae[,(2+ihrz)]/null.mae[,'ae']
}

##
## Average the error by season 
frcst.years = unique(frcst.rel[,'year'])
nfrcst.years = length(frcst.years)
mean.rel = array(NA,c(nfrcst.years,nhorizon))
colnames(mean.rel) = paste0('horizon',1:nhorizon)
rownames(mean.rel) = frcst.years

for (iyear in 1:nfrcst.years) {
	for (ihrz in 1:nhorizon) {
		idx = which(frcst.rel[,'year'] == frcst.years[iyear])
		mean.rel[iyear,ihrz] = mean(frcst.rel[idx,(2+ihrz)])
		mean.rel[iyear,ihrz] = round(mean.rel[iyear,ihrz],digits=2)
	}
}


results = list(deng.obsrv=deng.obsrv,deng.null=deng.null,deng.frcst=deng.frcst,deng.lower=deng.lower,deng.upper=deng.upper,null.mae=null.mae,frcst.mae=frcst.mae,frcst.rel=frcst.rel,mean.rel=mean.rel)

# save
filename = paste0(mod_name[[1]],'.arima.v',nhorizon,'.RData')
save(results,file=filename)
##
## Now let's plot
##

if (pdf == TRUE) {
	filename = "FIGURE9.PDF"
	pdf(file=filename,width=20,height=8)
} 

ncol = max(nhorizon/2,1)
par(mfrow = c(2,ncol), mar = c(5, 5, 2, 1))

colvec = c("red", "#386cb0", "#7fc97f", "#beaed4")


nts = dim(deng.obsrv)[1]

ymin = 0
ymax = max(deng.frcst[, 3:(2 + nhorizon)], deng.lower[, 3:(2 + nhorizon)], deng.upper[, 3:(2 + nhorizon)], deng.obsrv[, "cases"], na.rm = TRUE)

xlab = paste0(deng.obsrv[1, "year"], " - ", deng.obsrv[nts, "year"])

#	
# This is the mapping between prediction horizon and data delay
#	
delay = (horizon-1)

delay = paste0(delay,'-week delay')

delay[1] = 'no delay'

for (ihrz in 1:nhorizon) {

	plot(1:nts,deng.obsrv[, "cases"], type = "l", lwd = 1, col = "black", ylim = c(ymin, ymax), ylab = "Cases", xlab = xlab, xaxt = "n")

	upper = deng.upper[, (2 + ihrz)]
	lower = deng.lower[, (2 + ihrz)]
	frcst = deng.frcst[, (2 + ihrz)]

	mycol = colvec[1]
	if (nhorizon <= 4) mycol = colvec[ihrz]
		
	polygon(c(1:nts, rev(1:nts)), c(upper, rev(lower)), col = alpha(mycol, 0.2), border = FALSE)

	lines(deng.null[, "cases"], type = "l", col = 'black', lwd = 2, lty=2) #"darkgrey"

	lines(deng.obsrv[, "cases"], type = "l", col = "black", lwd = 2)

	lines(frcst, type = "l", col = alpha(mycol, 0.9), lwd = 3, xlab = "", ylab = "")

	legend("topleft", c(delay[ihrz]), bty = "n",text.col=alpha(mycol, 1.0),cex=1.2)
	#legend("topleft", c(delay[ihrz],mean.rel[,ihrz]), bty = "n",text.col=alpha(mycol, 1.0),cex=1.2)

	axis(1, at = 1:nts, label = deng.obsrv[, "week"], las = 1, cex = 0.5)
	legend('topright',legend=c('reported','hist. ave.'),lty=c(1,2),bty='n',lwd=2,cex=1.2)
}

if (pdf == TRUE)
	dev.off()

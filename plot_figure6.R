##
## plotting the results of an ARIMA analysis for a given country
## This version also analyses and plots the National results - direct and aggregate
## Make Figure 6


rm(list = ls())

library("Hmisc")
library("RColorBrewer")
require(DICE)

country = "BR"

if (country == "BR") {
	iso3 = "BRA"
	country.name = "BRAZIL"
} else if (country == "TH") {
	iso3 = "THA"
	country.name = "THAILAND"
} else if (country == "MX") {
	iso3 = "MEX"
	country.name = "MEXICO"
} else {
	iso3 = "NAN"
	country.name = "unknown"
}

complete_data=get.DICE.data(data_source=NULL, mod_level = 2, mod_name=c(NAME_2="BR"), fit_names='all', fit_level=3, year = 2015, nperiodsFit = 12, db_opts=list(DICE_db="predsci", CDC_server=FALSE), disease="dengue", all_years_flag=T, all_cad_clim=T)

mydata = complete_data$mydata

rain = all_years_clim$fit$precip

temp = all_years_clim$fit$temp

sh   = all_years_clim$fit$sh

names(rain) = names(temp) = names(sh) = mydata$fit$name


filename = paste0(country, ".arima.RData")

load(filename)

frcst.mae = results$frcst.mae

frcst.mre = results$frcst.mre

frcst.rel = results$frcst.rel2


## Repeat for national

nfrcst.mae = results$nfrcst.mae

nfrcst.mre = results$nfrcst.mre

nfrcst.rel = results$nfrcst.rel2

nhorizon = dim(frcst.mae)[1]

nstates = dim(frcst.mae)[2]

narima = dim(frcst.mae)[3]

arima.names = dimnames(frcst.mae)[[3]]

state.names = dimnames(frcst.mae)[[2]]

# Repeat for national

nnatl = dim(nfrcst.mae)[2]

natl.names = dimnames(nfrcst.mae)[[2]]


## Now load the results without the covar
filename = paste0(country, ".nocovar.arima.RData")

load(filename)

frcst.mae.nocovar = results$frcst.mae

frcst.mre.nocovar = results$frcst.mre

frcst.rel.nocovar = results$frcst.rel2

## Repeat for national

nfrcst.mae.nocovar = results$nfrcst.mae

nfrcst.mre.nocovar = results$nfrcst.mre

nfrcst.rel.nocovar = results$nfrcst.rel2


narima.nocovar = dim(frcst.mae.nocovar)[3]

arima.names.nocovar = dimnames(frcst.mae.nocovar)[[3]]


## calculate the ratio of MAE with covar to MAE without covar.
## The one with COVAR has lags of 0,1,2 - which is why dim(3) is 9 and not 3
##
## 
frcst.covar.nocovar = frcst.mae 
dim.vec = dim(frcst.covar.nocovar)

frcst.covar.nocovar[1:dim.vec[1], 1:dim.vec[2], 1:dim.vec[3]] = NA

for (j in 1:3) {
	for (k in 1:3) {
	frcst.covar.nocovar[,,((j-1)*3+k)] = frcst.mae[,,((j-1)*3+k)]/frcst.mae.nocovar[,,j]
	}
}

ratio = array(0, c(dim(frcst.covar.nocovar)[[2]],9))
rownames(ratio) = dimnames(frcst.covar.nocovar)[[2]]
colnames(ratio) = c(paste0("Lag",0:2), 'mean precip', 'sd precip', 'mean temp', 'sd temp', 'mean SH', 'sd SH')

for(i in 1:dim(frcst.covar.nocovar)[2]) {
	for (j in 1:3) ratio[i, j] = median(frcst.covar.nocovar[,i,6+j])
	## The order of the sattes is not the same need to find the right one
	k = which(mydata$fit$name == rownames(ratio)[i])
	ratio[i, (4)] = mean(rain[,k])
	ratio[i, (5)] = sd(  rain[,k])
	ratio[i, (6)] = mean(  temp[,k])
	ratio[i, (7)] = sd(temp[,k])
	ratio[i, (8)] = mean(  sh[,k])
	ratio[i, (9)] = sd(  sh[,k])
}

iorder = order(ratio[,1], decreasing = FALSE)
ratio = ratio[iorder,]
write.csv(ratio, file="BR.covar.info.csv")


## Begin plotting of state level results

narima.tot = narima+narima.nocovar

rf <- colorRampPalette(rev(brewer.pal(11, "Spectral"))) # make colors
colvec <- rf(narima.tot)

colmat <- array(NA, c(narima.tot, nhorizon))
rownames(colmat) = c(arima.names,arima.names.nocovar)
colnames(colmat) = paste0("horizon=", 1:nhorizon)


colmat[, 1:nhorizon] = rep(alpha(colvec, 0.75), nhorizon) ##c(alpha(colvec,0.8),alpha(colvec,0.8),alpha(colvec,0.8),alpha(colvec,0.8))

colvec.covar <- rf(narima)
colmat.covar <- array(NA, c(narima, nhorizon))
rownames(colmat.covar) = c(arima.names)
colnames(colmat.covar) = paste0("horizon=", 1:nhorizon)


colmat.covar[, 1:nhorizon] = rep(alpha(colvec.covar, 0.75), nhorizon) ##c(alpha(colvec,0.8),alpha(colvec,0.8),alpha(colvec,0.8),alpha(colvec,0.8))

cex = rep(1.5, narima.tot)

cex.axis = 0.8

frcst.mre[is.infinite(frcst.mre)] <- NA

nstates1 = nstates + 1
nnatl1 = nnatl + 1


filename = "figure6.pdf"

pdf(file = filename, width = 10, height = 12)

par(mfrow = c(1, 3), mar = c(5, 6.5, 2, 1))

xmin = round(min(frcst.mae, frcst.mae.nocovar, na.rm = TRUE), digits = 1)
xmax = round(max(frcst.mae, frcst.mae.nocovar, na.rm = TRUE), digits = 1)
xmin = floor(xmin)
xmax = round(xmax)

plot(x = frcst.mae[1, , 1], y = 1:nstates, col = colmat[1, 1], xlim = range(frcst.mae, na.rm = TRUE), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = expression("MAE of ARIMA Models"), 
	yaxt = "n", log = "x") #xlim=c(xmin,xmax)
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = frcst.mae[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
	for (j in 1:narima.nocovar) lines(x = frcst.mae.nocovar[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[(j+narima), i], type = "p", pch = 17, ylab = "", xlab = "", cex = cex[1])
	#for (j in 1:narima) lines(x = frcst.mae[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
}

abline(h = 1:nstates, col = alpha("grey", 0.3))
axis(2, at = 1:nstates, label = state.names, cex.axis = cex.axis, las = 2)


legend("topleft", legend =c( paste0(arima.names, 12),paste0(arima.names.nocovar, 12)), text.col = colvec, bty = "n")

xmin = round(min(frcst.mre, frcst.mre.nocovar, na.rm = TRUE), digits = 1)
xmax = round(max(frcst.mre, frcst.mre.nocovar, na.rm = TRUE), digits = 1)
#xmin = floor(xmin)
xmax = round(xmax)

tit = paste0("ARIMA ANALYSIS FOR: ", country.name)

plot(x = frcst.mre[1, , 1], y = 1:nstates, col = colmat[1, 1], xlim = range(frcst.mre, na.rm = TRUE), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = expression("MRAE of ARIMA Models"), 
	yaxt = "n", main = tit, log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = frcst.mre[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
	for (j in 1:narima.nocovar) lines(x = frcst.mre.nocovar[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[(j+narima), i], type = "p", pch = 17, ylab = "", xlab = "", cex = cex[1])	
}
abline(h = 1:nstates, col = alpha("grey", 0.3))
axis(2, at = 1:nstates, label = state.names, cex.axis = 0.8, las = 2)

xmin = round(min(frcst.rel, frcst.rel.nocovar, na.rm = TRUE), digits = 2)
xmax = round(max(frcst.rel, frcst.rel.nocovar, na.rm = TRUE), digits = 2)

# plot(x = frcst.rel[1, , 1], y = 1:nstates, col = colmat[1, 1], xlim = range(frcst.rel, na.rm = TRUE), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = expression("MAE(ARIMA)/MAE(NULL)"), 
	# yaxt = "n", log = "x")
# for (i in 1:nhorizon) {
	# for (j in 1:narima) lines(x = frcst.rel[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
	# for (j in 1:narima.nocovar) lines(x = frcst.rel.nocovar[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[(j + narima), i], type = "p", pch = 17, ylab = "", xlab = "", 
		# cex = cex[1])

# }
# abline(h = 1:nstates, col = alpha("grey", 0.3))
# abline(v = 1, col = alpha("black", 0.8), lty = 2)
# axis(2, at = 1:nstates, label = state.names, cex.axis = cex.axis, las = 2)



## ratio of MAE with COVAR to MAE without COVAR
##

plot(x = frcst.covar.nocovar[1,,1], y = 1:nstates, col = colmat.covar[1, 1], xlim = range(frcst.covar.nocovar, na.rm = TRUE), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = expression("MAE-COVAR/MAE-NOCOVAR"), 
	yaxt = "n", log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = frcst.covar.nocovar[i, ,j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat.covar[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
}
legend("topright", legend =paste0(arima.names, 12), text.col = colvec.covar, bty = "n")
abline(h = 1:nstates, col = alpha("grey", 0.3))
abline(v = 1, col = alpha("black", 0.8), lty = 2)
axis(2, at = 1:nstates, label = state.names, cex.axis = cex.axis, las = 2)

dev.off()


filename = paste0(country,'/',country,".natl.arima.pdf")

pdf(file = filename, width = 10, height = 12)

par(mfrow = c(1, 3), mar = c(5, 6.5, 2, 1))

xmin = round(min(nfrcst.mae, nfrcst.mae.nocovar, na.rm = TRUE), digits = 1)
xmax = round(max(nfrcst.mae, nfrcst.mae.nocovar, na.rm = TRUE), digits = 1)


if (country == "MEX") {
	xmin = 750
	xmax = 4000
}


plot(x = nfrcst.mae[1, , 1], y = 1:nnatl, col = colmat[1, 1], xlim = range(nfrcst.mae, nfrcst.mae.nocovar, na.rm = TRUE), ylim = c(1, nnatl1), type = "p", pch = 19, ylab = "", xlab = expression("MAE of ARIMA Models"), 
	yaxt = "n", log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = nfrcst.mae[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = 1.5)
	for (j in 1:narima.nocovar) lines(x = nfrcst.mae.nocovar[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[(j + narima), i], type = "p", pch = 17, ylab = "", xlab = "", 
		cex = cex[1])	
}

abline(h = 1:nnatl, col = alpha("grey", 0.3))
axis(2, at = 1:nnatl, label = natl.names, cex.axis = 1, las = 2)

legend("topleft", legend =c( paste0(arima.names, 12),paste0(arima.names.nocovar, 12)), text.col = colvec, bty = "n")

xmin = round(min(nfrcst.mre, nfrcst.mre.nocovar, na.rm = TRUE), digits = 2)
xmax = round(max(nfrcst.mre, nfrcst.mre.nocovar, na.rm = TRUE), digits = 2)

xmax = round(xmax)

if (country == "MEX") {
	xmin = 0.3
	xmax = 2.5
}



tit = paste0("NATIONAL LEVEL ARIMA ANALYSIS: ", country.name)

plot(x = nfrcst.mre[1, , 1], y = 1:nnatl, col = colmat[1, 1], xlim = range(nfrcst.mre, na.rm = TRUE), ylim = c(1, nnatl1), type = "p", pch = 19, ylab = "", xlab = expression("MRAE of ARIMA Models"), 
	yaxt = "n", main = tit, log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = nfrcst.mre[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = 1.5)
	for (j in 1:narima.nocovar) lines(x = nfrcst.mre.nocovar[i, , j], y = ((1:nnatl)+ (i - 1) * 0.2), col = colmat[(j + narima), i], type = "p", pch = 17, ylab = "", xlab = "", 
		cex = cex[1])
}

abline(h = 1:nnatl, col = alpha("grey", 0.3))
axis(2, at = 1:nnatl, label = natl.names, cex.axis = 1, las = 2)

xmin = round(min(nfrcst.rel, nfrcst.rel.nocovar, na.rm = TRUE), digits = 2)
xmax = round(max(nfrcst.rel, nfrcst.rel.nocovar, na.rm = TRUE), digits = 2)

if (country == "MEX") {

	xmin = 0.3
	xmax = 1.7

}
plot(x = nfrcst.rel[1, , 1], y = 1:nnatl, col = colmat[1, 1], xlim = range(nfrcst.rel, na.rm = TRUE), ylim = c(1, nnatl1), type = "p", pch = 19, ylab = "", xlab = expression("MAE(ARIMA)/MAE(NULL)"), 
	yaxt = "n", log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = nfrcst.rel[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = 1.5)
	for (j in 1:narima.nocovar) lines(x = nfrcst.rel.nocovar[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[(j + narima), i], type = "p", pch = 17, ylab = "", xlab = "", 
		cex = cex[1])	
}
abline(h = 1:nnatl, col = alpha("grey", 0.3))
abline(v = 1, col = alpha("black", 0.8), lty = 2)
axis(2, at = 1:nnatl, label = natl.names, cex.axis = 1, las = 2)

dev.off()


##
## Create tables with mean values of the Relative MAE of the different ARIMA models


nNULL = 1

stats = array(NA, c((nhorizon * nNULL), (narima.tot)))
colnames(stats) = c(arima.names, arima.names.nocovar)
rownames(stats) = rep(paste0("horizon=", 1:nhorizon), each = nNULL)

icount = 1
for (i in 1:nhorizon) {
	for (j in 1:narima) {
		stats[icount, (j)] = mean(frcst.rel[i, , j], na.rm = TRUE)
	}
	icount = icount + nNULL
}

icount = 1

for (i in 1:nhorizon) {
	for (j in 1:narima.nocovar) {
		stats[icount, (j+narima)] = mean(frcst.rel.nocovar[i, , j], na.rm = TRUE)
	}	
	icount = icount + nNULL
}

nstats = array(NA, c((nhorizon * nNULL * nnatl), ( narima.tot)))
colnames(nstats) = c(arima.names, arima.names.nocovar)
rownames(nstats) = rep(paste0("horizon=", 1:nhorizon), each = (nnatl * nNULL))

icount = 1
## Repeat for the National Aggregate and Direct
for (i in 1:nhorizon) {
	for (j in 1:narima) {
		nstats[icount, (j)] = mean(nfrcst.rel[i, "Direct", j], na.rm = TRUE)
		nstats[(icount + 1), (j)] = mean(nfrcst.rel[i, "Aggregate", j], na.rm = TRUE)

	}
	icount = icount + nNULL * 2
	
}


icount = 1
for (i in 1:nhorizon) {
	for (j in 1:narima.nocovar) {
		nstats[icount, (j + narima)] = mean(nfrcst.rel.nocovar[i, "Direct", j], na.rm = TRUE)
		nstats[(icount + 1), (j + narima)] = mean(nfrcst.rel.nocovar[i, "Aggregate", j], na.rm = TRUE)

	}
	icount = icount + nNULL * 2
	
}

tstats = rbind(stats, nstats)

filename = paste0(country,'/',country, ".arima.csv")
write.csv(tstats, file = filename)




##
## plotting the results of an ARIMA analysis for a given country
## This version also analyses and plots the National results - direct and aggregate
##


rm(list = ls())

library("Hmisc")
library("RColorBrewer")

country = "TH"

if (country == "BR") {
	iso3 = "BRA"
	country.name = "BRAZIL"
} else if (country == "TH") {
	iso3 = "THA"
	country.name = "THAILAND"
} else if (country == "MX") {
	iso3 = "MEX"
	country.name = "MEXICO"
} else if (country == "SG") {
	iso3 = "SGP"
	country.name = "Singapore"
} else {
	iso3 = ""
	country.name = "Unknown"	
}



filename = paste0(country, ".arima.v4.RData")

# for SG 
#filename = paste0(country, ".arima.v10.RData")

load(filename)


frcst.mae = results$frcst.mae

frcst.mre = results$frcst.mre

frcst.rel1 = results$frcst.rel1

frcst.rel2 = results$frcst.rel2


## Repeat for national

nfrcst.mae = results$nfrcst.mae

nfrcst.mre = results$nfrcst.mre

nfrcst.rel1 = results$nfrcst.rel1

nfrcst.rel2 = results$nfrcst.rel2

nhorizon = dim(frcst.mae)[1]

nhorizon = 4

nstates = dim(frcst.mae)[2]

narima = dim(frcst.mae)[3]

arima.names = dimnames(frcst.mae)[[3]]

state.names = dimnames(frcst.mae)[[2]]

# Repeat for national

nnatl = dim(nfrcst.mae)[2]

natl.names = dimnames(nfrcst.mae)[[2]]


## Begin plotting of state level results

rf <- colorRampPalette(rev(brewer.pal(11, "Spectral"))) # make colors
colvec <- rf(narima)

colmat <- array(NA, c(narima, nhorizon))
rownames(colmat) = arima.names
colnames(colmat) = paste0("horizon=", 1:nhorizon)

#colmat[,1:nhorizon] = c(alpha(colvec,0.8),alpha(colvec,0.6),alpha(colvec,0.4),alpha(colvec,0.2))

colmat[, 1:nhorizon] = rep(alpha(colvec, 0.75), nhorizon) ##c(alpha(colvec,0.8),alpha(colvec,0.8),alpha(colvec,0.8),alpha(colvec,0.8))

cex = rep(1.5, narima)

cex.axis = 0.9

frcst.mre[is.infinite(frcst.mre)] <- NA

nstates1 = nstates + 1
nnatl1 = nnatl + 1


filename = paste0(country, "arima.pdf")

pdf(file = filename, width = 10, height = 12)

par(mfrow = c(1, 3), mar = c(5, 6.5, 2, 1))

xmin = round(min(frcst.mae, na.rm = TRUE), digit = 1)
xmax = round(max(frcst.mae, na.rm = TRUE), digits = 1)
xmin = floor(xmin)
xmax = round(xmax)

plot(x = frcst.mae[1, , 1], y = 1:nstates, col = colmat[1, 1], xlim = range(frcst.mae, na.rm = TRUE), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = expression("MAE of ARIMA Models"), 
	yaxt = "n", log = "x") #xlim=c(xmin,xmax)
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = frcst.mae[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
}

abline(h = 1:nstates, col = alpha("grey", 0.3))
axis(2, at = 1:nstates, label = state.names, cex.axis = cex.axis, las = 2)

legend = expression(arima.names[12])
#legend("topleft", legend = paste0(arima.names, 12), text.col = colvec, bty = "n")

xmin = round(min(frcst.mre, na.rm = TRUE), digit = 1)
xmax = round(max(frcst.mre, na.rm = TRUE), digits = 1)
#xmin = floor(xmin)
xmax = round(xmax)



if (country == 'MX') {
	country.name = "Mexico"
}
tit = paste0("ARIMA ANALYSYS FOR: ", country.name)

plot(x = frcst.mre[1, , 1], y = 1:nstates, col = colmat[1, 1], xlim = range(frcst.mre, na.rm = TRUE), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = expression("MRAE of ARIMA Models"), 
	yaxt = "n", main = tit, log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = frcst.mre[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
}
abline(h = 1:nstates, col = alpha("grey", 0.3))
axis(2, at = 1:nstates, label = state.names, cex.axis = 0.8, las = 2)

xmin = round(min(frcst.rel1, na.rm = TRUE), digits = 2)
xmax = round(max(frcst.rel1, na.rm = TRUE), digits = 2)

plot(x = frcst.rel1[1, , 1], y = 1:nstates, col = colmat[1, 1], xlim = range(frcst.rel1, frcst.rel2, na.rm = TRUE), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = expression("MAE(ARIMA)/MAE(NULL)"), 
	yaxt = "n", log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) {
		#lines(x = frcst.rel1[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
		lines(x = frcst.rel2[i, , j], y = ((1:nstates) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
	}
}
abline(h = 1:nstates, col = alpha("grey", 0.3))
abline(v = 1, col = alpha("black", 0.8), lty = 2)
axis(2, at = 1:nstates, label = state.names, cex.axis = cex.axis, las = 2)
legend("topright", legend = paste0(arima.names, 12), text.col = colvec, bty = "n")

dev.off()

##
## Now plot seperately for each forecast month
##

cex = rep(1.2, narima)

cex.axis = 0.8
if (country == "BR") 
	cex.axis = 0.8
if (country == "TH") 
	cex.axis = 0.4

filename = paste0(country, "arima.v2.pdf")

pdf(file = filename, width = 10, height = 12)

par(mfcol = c(4, 3), mar = c(4.5, 6.5, 1, 1))

xmin = round(min(frcst.mae, na.rm = TRUE), digits = 1)
xmax = round(max(frcst.mae, na.rm = TRUE), digits = 1)
xmin = floor(xmin)
xmax = round(xmax)


for (i in 1:nhorizon) {
	par(font = 1)
	xlab = ""
	if (i == nhorizon) 
		xlab = expression("MAE of ARIMA Models")
	plot(x = frcst.mae[i, , 1], y = 1:nstates, col = colmat[1, 1], xlim = c(xmin, xmax), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = xlab, yaxt = "n", log = "x")
	mean.arima = rep(0, narima)
	for (j in 1:narima) {
		lines(x = frcst.mae[i, , j], y = 1:nstates, col = colmat[j, 1], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
		mean.arima[j] = mean(frcst.mae[i, , j], na.rm = TRUE)
		mean.arima[j] = round(mean.arima[j])
	}
	abline(h = 1:nstates, col = alpha("grey", 0.3))
	axis(2, at = 1:nstates, label = state.names, cex.axis = cex.axis, las = 2)
	legend("topleft", legend = paste0(i, "-month forecast"), text.col = "black", bty = "n", cex = 1.1)
	par(font = 2)
	legend("bottomright", legend = mean.arima, text.col = colmat[, 1], bty = "n", cex = 1)
}

xmin = round(min(frcst.mre, na.rm = TRUE), digit = 1)
xmax = round(max(frcst.mre, na.rm = TRUE), digits = 1)
#xmin = floor(xmin)
xmax = ceiling(xmax)

for (i in 1:nhorizon) {
	par(font = 1)
	xlab = ""
	if (i == nhorizon) 
		xlab = expression("MRAE of ARIMA Models")
	plot(x = frcst.mre[i, , 1], y = 1:nstates, col = colmat[1, 1], xlim = range(xmin, xmax), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = xlab, yaxt = "n", 
		log = "x")
	for (j in 1:narima) {
		lines(x = frcst.mre[i, , j], y = 1:nstates, col = colmat[j, 1], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
		mean.arima[j] = mean(frcst.mre[i, , j], na.rm = TRUE)
		mean.arima[j] = round(mean.arima[j], digits = 1)

	}
	abline(h = 1:nstates, col = alpha("grey", 0.3))
	axis(2, at = 1:nstates, label = state.names, cex.axis = cex.axis, las = 2)
	par(font = 2)
	legend("bottomright", legend = mean.arima, text.col = colmat[, 1], bty = "n", cex = 1)
}

xmin = round(min(frcst.rel1, na.rm = TRUE), digits = 2)
xmax = round(max(frcst.rel2, na.rm = TRUE), digits = 1)
xmax = ceiling(xmax)
print(xmax)
if (country == "TH") 
	xmax = 5

for (i in 1:nhorizon) {
	xlab = ""
	if (i == nhorizon) 
		xlab = expression("MAE(ARIMA)/MAE(NULL)")
	plot(x = frcst.rel1[1, , 1], y = 1:nstates, col = colmat[1, 1], xlim = c(xmin, xmax), ylim = c(1, nstates1), type = "p", pch = 19, ylab = "", xlab = xlab, yaxt = "n", log = "x")
	mean.arima1 = rep(0, narima)
	mean.arima2 = rep(0, narima)
	for (j in 1:narima) {
		lines(x = frcst.rel1[i, , j], y = 1:nstates, col = colmat[j, 1], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
		lines(x = frcst.rel2[i, , j], y = 1:nstates, col = colmat[j, 1], type = "p", pch = 17, ylab = "", xlab = "", cex = cex[1])
		mean.arima1[j] = mean(frcst.rel1[i, , j], na.rm = TRUE)
		mean.arima1[j] = round(mean.arima1[j], digits = 2)
		mean.arima2[j] = mean(frcst.rel2[i, , j], na.rm = TRUE)
		mean.arima2[j] = round(mean.arima2[j], digits = 2)
	}
	abline(h = 1:nstates, col = alpha("grey", 0.3))
	abline(v = 1, col = alpha("black", 0.8), lty = 2)
	axis(2, at = 1:nstates, label = state.names, cex.axis = cex.axis, las = 2)
	legend = paste0(arima.names, 12)
	if (i == 1) 
		legend("topright", legend = legend, text.col = colvec, bty = "n", cex = 1)
	legend = paste0(mean.arima1, "/", mean.arima2)
	par(font = 2)
	legend("bottomright", legend = legend, text.col = colmat[, 1], bty = "n", cex = 1)
}

dev.off()


##
## Plot the national direct and indirect
##


filename = paste0(country, ".natl.arima.pdf")

pdf(file = filename, width = 10, height = 12)

par(mfrow = c(1, 3), mar = c(5, 6.5, 2, 1))

xmin = round(min(nfrcst.mae, na.rm = TRUE), digit = 1)
xmax = round(max(nfrcst.mae, na.rm = TRUE), digits = 1)
#xmin = floor(xmin)
#xmax = round(xmax)

if (country.name == "MEXICO") {
	xmin = 750
	xmax = 4000
}


plot(x = nfrcst.mae[1, , 1], y = 1:nnatl, col = colmat[1, 1], xlim = range(nfrcst.mae, na.rm = TRUE), ylim = c(1, nnatl1), type = "p", pch = 19, ylab = "", xlab = expression("MAE of ARIMA Models"), 
	yaxt = "n", log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = nfrcst.mae[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = 1.5)
}

abline(h = 1:nnatl, col = alpha("grey", 0.3))
axis(2, at = 1:nnatl, label = natl.names, cex.axis = 1, las = 2)

legend = expression(arima.names[12])
legend("topleft", legend = paste0(arima.names, 12), text.col = colvec, bty = "n")

xmin = round(min(nfrcst.mre, na.rm = TRUE), digit = 2)
xmax = round(max(nfrcst.mre, na.rm = TRUE), digits = 2)
#xmin = floor(xmin)
xmax = round(xmax)

if (country.name == "MEXICO") {
	xmin = 0.3
	xmax = 2.5
}



tit = paste0("NATIONAL LEVEL ARIMA ANALYSYS: ", country.name)

plot(x = nfrcst.mre[1, , 1], y = 1:nnatl, col = colmat[1, 1], xlim = range(nfrcst.mre, na.rm = TRUE), ylim = c(1, nnatl1), type = "p", pch = 19, ylab = "", xlab = expression("MRAE of ARIMA Models"), 
	yaxt = "n", main = tit, log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = nfrcst.mre[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = 1.5)
}
abline(h = 1:nnatl, col = alpha("grey", 0.3))
axis(2, at = 1:nnatl, label = natl.names, cex.axis = 1, las = 2)

xmin = round(min(nfrcst.rel1, na.rm = TRUE), digits = 2)
xmax = round(max(nfrcst.rel1, na.rm = TRUE), digits = 2)

if (country.name == "MEXICO") {

	xmin = 0.3
	xmax = 1.7

}
plot(x = nfrcst.rel1[1, , 1], y = 1:nnatl, col = colmat[1, 1], xlim = range(nfrcst.rel1, nfrcst.rel2, na.rm = TRUE), ylim = c(1, nnatl1), type = "p", pch = 19, ylab = "", xlab = expression("MAE(ARIMA)/MAE(NULL)"), 
	yaxt = "n", log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) {
		lines(x = nfrcst.rel1[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 19, ylab = "", xlab = "", cex = 1.5)
		lines(x = nfrcst.rel2[i, , j], y = ((1:nnatl) + (i - 1) * 0.2), col = colmat[j, i], type = "p", pch = 17, ylab = "", xlab = "", cex = 1.5)
	}
}
abline(h = 1:nnatl, col = alpha("grey", 0.3))
abline(v = 1, col = alpha("black", 0.8), lty = 2)
axis(2, at = 1:nnatl, label = natl.names, cex.axis = 1, las = 2)

dev.off()


##
## Create tables with mean values of the Relative MAE of the different ARIMA models

nNULL = 2

stats = array(NA, c((nhorizon * nNULL), (1 + narima)))
colnames(stats) = c("NULL Model", arima.names)
rownames(stats) = rep(paste0("horizon=", 1:nhorizon), each = nNULL)

icount = 1
for (i in 1:nhorizon) {
	for (j in 1:narima) {
		stats[icount:(icount + nNULL - 1), 1] = 1:nNULL
		stats[icount, (j + 1)] = mean(frcst.rel1[i, , j], na.rm = TRUE)
		stats[(icount + 1), (j + 1)] = mean(frcst.rel2[i, , j], na.rm = TRUE)

	}
	icount = icount + nNULL
}

nstats = array(NA, c((nhorizon * nNULL * nnatl), (1 + narima)))
colnames(nstats) = c("NULL Model", arima.names)
rownames(nstats) = rep(paste0("horizon=", 1:nhorizon), each = (nnatl * nNULL))

icount = 1
## Repeat for the National Aggregate and Direct
for (i in 1:nhorizon) {
	for (j in 1:narima) {
		nstats[icount:(icount + nNULL * 2 - 1), 1] = paste0(rep(1:nNULL, each = 2), "-", natl.names)

		nstats[icount, (j + 1)] = mean(nfrcst.rel1[i, "Direct", j], na.rm = TRUE)
		nstats[(icount + 1), (j + 1)] = mean(nfrcst.rel1[i, "Aggregate", j], na.rm = TRUE)

		nstats[(icount + 2), (j + 1)] = mean(nfrcst.rel2[i, "Direct", j], na.rm = TRUE)
		nstats[(icount + 3), (j + 1)] = mean(nfrcst.rel2[i, "Aggregate", j], na.rm = TRUE)


	}
	icount = icount + nNULL * 2
}

tstats = rbind(stats, nstats)

filename = paste0(country, ".arima.csv")
write.csv(tstats, file = filename)





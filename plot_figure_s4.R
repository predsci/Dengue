##
## plotting the results of an ARIMA analysis for a given country
## This version also analyses and plots the National results - direct and aggregate
## Make Figure S4 - Singapore


rm(list = ls())

library("Hmisc")
library("RColorBrewer")

country = "SG"

iso3 = "SGP"
country.name = "Singapore"

filename = paste0(country, ".arima.RData")

load(filename)


frcst.mae = results$frcst.mae

frcst.mre = results$frcst.mre

frcst.rel1 = results$frcst.rel1

frcst.rel2 = results$frcst.rel2


nhorizon = dim(frcst.mae)[1]


nstates = 1

narima = dim(frcst.mae)[2]

arima.names = dimnames(frcst.mae)[[2]]

## Begin plotting of state level results

rf <- colorRampPalette(rev(brewer.pal(11, "Spectral"))) # make colors
colvec <- rf(narima)

colmat <- array(NA, c(narima, nhorizon))
rownames(colmat) = arima.names
colnames(colmat) = paste0("horizon=", 1:nhorizon)

#colmat[,1:nhorizon] = c(alpha(colvec,0.8),alpha(colvec,0.6),alpha(colvec,0.4),alpha(colvec,0.2))

colmat[, 1:nhorizon] = rep(colvec, nhorizon) ##c(alpha(colvec,0.8),alpha(colvec,0.8),alpha(colvec,0.8),alpha(colvec,0.8))

cex = rep(1.5, narima)

cex.axis = 0.9

frcst.mre[is.infinite(frcst.mre)] <- NA

filename = "figure_s4.pdf"

pdf(file = filename, width = 10, height = 12)

par(mfrow = c(1, 3), mar = c(5, 6.5, 2, 1))

xmin = round(min(frcst.mae, na.rm = TRUE), digit = 1)
xmax = round(max(frcst.mae, na.rm = TRUE), digits = 1)
xmin = floor(xmin)
xmax = round(xmax)

plot(x = frcst.mae[1, 1], y = 1, col = colvec[1], xlim = range(frcst.mae, na.rm = TRUE), ylim = c(1, nhorizon), type = "p", pch = 19, ylab = "", xlab = expression("MAE of ARIMA Models"), 
	yaxt = "n", log = "x") #xlim=c(xmin,xmax)
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = frcst.mae[i, j], y = i, col = colvec[j], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
}

abline(h = 1:nhorizon, col = alpha("grey", 0.3))
axis(2, at = 1:nhorizon, label = paste0(1:nhorizon, " weeks"), cex.axis = cex.axis, las = 2)

#legend = expression(arima.names[12])
#legend("topleft", legend = paste0(arima.names, 12), text.col = colvec, bty = "n")

xmin = round(min(frcst.mre, na.rm = TRUE), digit = 1)
xmax = round(max(frcst.mre, na.rm = TRUE), digits = 1)
#xmin = floor(xmin)
xmax = round(xmax)

xmin = 0.1
xmax = 1.2

tit = paste0("ARIMA ANALYSIS: ", country.name)

plot(x = frcst.mre[1, 1], y = 1, col = colvec[1], xlim = c(xmin,xmax), ylim = c(1, nhorizon), type = "p", pch = 19, ylab = "", xlab = expression("MRAE of ARIMA Models"), 
	yaxt = "n", main = tit, log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) lines(x = frcst.mre[i, j], y = i, col = colvec[j], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
}
abline(h = 1:nhorizon, col = alpha("grey", 0.3))
axis(2, at = 1:nhorizon, label = paste0(1:nhorizon, " weeks"), cex.axis = cex.axis, las = 2)

xmin = round(min(frcst.rel2, na.rm = TRUE), digits = 2)
xmax = round(max(frcst.rel2, na.rm = TRUE), digits = 2)
xmin = 0.15
xmax = 1.2
plot(x = frcst.rel2[1, 1], y = 1, col = colvec[1], xlim = c(xmin,xmax), ylim = c(1, nhorizon),  type = "p", pch = 19, ylab = "", xlab = expression("MAE(ARIMA)/MAE(NULL)"), 
	yaxt = "n", log = "x")
for (i in 1:nhorizon) {
	for (j in 1:narima) {		
		lines(x = frcst.rel2[i, j], y = i, col = colvec[j], type = "p", pch = 19, ylab = "", xlab = "", cex = cex[1])
	}
}

abline(h = 1:nhorizon, col = alpha("grey", 0.3))
axis(2, at = 1:nhorizon, label = paste0(1:nhorizon, " weeks"), cex.axis = cex.axis, las = 2)
legend("topleft", legend = paste0(arima.names, 52), text.col = colvec, bty = "n")

dev.off()





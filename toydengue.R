##
## Multi serotype vector-host SEIR model
##

rm(list=ls())

dyn.load("nserovechost.so") 

if (is.loaded("nserovechost")) {
	cat("\n Loaded n-serotype Vector-Host Code \n")
} else	{
	cat(" nserovechost.so has not loaded ",'\n')
	cat(" use the ./compile.pyf script to compile and then try to re-run R code")
	q(save="yes")
}
	
	
weeks_per_year = 52

days_per_year = 365


#
# Number of sertypes >= 1

ns = 4

#
# Human total population

NH = 2e7

#
# Number of years to simulate (converted to weeks)
#

nyears = 30
nweeks = nyears * weeks_per_year


# Host birth = death rate (1/60 to 1/70)

muH = 1./60. / days_per_year


#
# Vector death rate (1/42 - 1/8 day^-1)

muV = 1.0 / 14.0 

#
# One over Average latent period in host (4-7 days)

sigmaH = 1.0 / 5.0

#
# One over Average latent period in vector (7-14 days)

sigmaV = 1.0 / 10.0

#
# One over Average Infectious period in host - this is a vector of length ns (one over 4-12 days)

gamma =1.0 / 8.0 *  array(1,ns)


#
# Transmission rate host -> vector - this is a vector of length ns (70 year^-1)

beta = 70.0 / days_per_year * array(1,ns)

#
# Transmission rate vector -> host - this is a vector of length ns (70 year^-1)

alpha = 70.0 / days_per_year * array(1,ns)


#
# Per capita infection-induced mortality probabilities - this is a vector of length ns (low or zero)

rho = 0.0 * array(1, ns)

#
# The average number of mosquitoes per person (2)

k  = 2

#
# Seasonality parameter for mosquito birth rate: [1- Va * cos(2*pi*time)] (0.0-0.05)

Va = 0.0

#
# Enhanced immunity to other serotypes. Complete immunity eta = 0, no immunity eta = 1. Partial immunity 0 < eta_i < 1. This is a vector of length ns

eta = 0.5 * array(1, ns)

#
# One over the Average length of time spent in state of enhanced immunity to other serotypes this is a vector of length ns. (2-9 months)

delta =  1. / days_per_year * array(1, ns)

#
# The effect of increased susceptibility to infection with a second serotype chi_j > 1 this is a vector of length ns


chi = 1.0 * array(1, ns)

#
#  Loss of antibody-dependent enhancement ADE. The period od ADE is 1/omega

omega = 1.0 / days_per_year * array(1, ns)

#
# Seed for random number generator - currently not used 

iseed = 12345

## --------------------------------
## Start of Initial Conditions

#
# Initial fraction of population that is fully susceptible to all ns serotypes 

fracS0 = 0.29 

#
# Initial number of Humans exposed to each serotype. This is an array of length ns

sEp = 10 * 1:ns

#
# Initial number of infectious with serotype i.  This is an array of length ns

sIp = 10 * 1:ns

#
# Initial Number of vectors exposed to serotype i. This is an arrray of length ns

Ve = array(0, ns)

#
# Initial number of vectors infectious with serotype i. This is an arrray of length ns 

Vi = array(0, ns)

## End of Initial Conditions 
## --------------------------------


## calculate R0[i]
##

R0 = k * alpha * beta * sigmaH * sigmaV /(muV *(gamma + muH)*(sigmaH + muH)*(sigmaV + muV))

print(R0)

#
# Allocate an array for the time series of the ns serotypes

rtn = array(0, c(ns, nweeks))

out <- .Fortran("nserovechost",ns=as.integer(ns), nweeks=as.integer(nweeks), NH = as.double(NH), muH=as.double(muH), muV = as.double(muV), sigmaH = as.double(sigmaH), sigmaV = as.double(sigmaV),gamma=as.double(gamma), beta=as.double(beta), alpha=as.double(alpha),rho=as.double(rho),k=as.double(k), Va=as.double(Va), eta=as.double(eta), delta=as.double(delta), chi=as.double(chi), omega=as.double(omega),fracS0=as.double(fracS0),sEp = as.double(sEp), sIp=as.double(sIp), Ve = as.double(Ve), Vi = as.double(Vi), iseed=as.integer(iseed), rtn=as.double(rtn))


rtn <- matrix(out$rtn,nc=nweeks)

# 
# A time vector in weeks - used for plotting will be converted to years 
tps = 1:nweeks

rtnTot = array(0,nweeks)

tps = tps / weeks_per_year

for (i in 1:nweeks) 
	rtnTot[i] = sum(rtn[1:ns,i])
	
rownames(rtn) = paste0('den',1:ns)

colvec = rainbow(ns)

par(mar=c(5,5,2,5),mfrow=c(2,1))
plot(tps,rtn[1,],type='l',col=colvec[1],xlab = 'Time (years)',ylab='Incidence (n-serotypes)',ylim=c(0,max(rtn)))
for (i in 1:ns) lines(tps,rtn[i,] ,col=colvec[i])
legend('topleft',rownames(rtn),bty='n',text.col=colvec)
legend('topright','total',bty='n',text.col='black')
par(new=TRUE)
plot(tps, rtnTot, col='black', type='l', xaxt='n', yaxt='n', xlab='',ylab='',lwd=2)
axis(4)
mtext(text='Total Incidence', side=4, line =2)

# Now plot again, after removing the initial large peak
#

iburn = 10 # In years
index = which(tps > iburn)
tpsS <- tps[index]
rtnS <- rtn[,index]
rtnTotS <- rtnTot[index]




par(mar=c(5,5,2,5))
plot(tpsS,rtnS[1,],type='l',col=colvec[1],xlab = 'Time (years)',ylab='Incidence (n-serotypes)',ylim=c(0,max(rtnS)))
for (i in 1:ns) lines(tpsS,rtnS[i,] ,col=colvec[i])
legend('topleft',rownames(rtn),bty='n',text.col=colvec)
legend('topright','total',bty='n',text.col='black')
par(new=TRUE)
plot(tpsS, rtnTotS, col='black', type='l', xaxt='n', yaxt='n', xlab='',ylab='',lwd=2)
axis(4)
mtext(text='Total Incidence', side=4, line =2)





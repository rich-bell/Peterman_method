#	Peterman method
# Script of running the Kalman filter with dlm package
# Written by Jeremy Collie and Rich Bell on 2022-Dec-01
#	Code written in R 3.6.2

library(dlm)
# Step 1: read in the Bocaccio stock-recruitment data and reshape the arrays

boca<-read.csv("bocaccio_2017.csv",header = TRUE)
head(boca)

#  SSB in mt, R is age-0 Recruitment in thousands

# make some exploratory plots
plot(boca$SSB_t,boca$R,xlab="SSB (mt)",ylab="Recruitment (thousands)")

# Recruitment deviations are estimated for 1954-2015 though the assessment begins in 1892
# ------------------------------------------------------------------------------------------------
# Step 2: create the S and Y variables and fit the time-invariant Ricker model
S <- as.numeric(boca$SSB_t)
Y <- log(boca$R/boca$SSB_t)
plot(S,Y,xlab="SSB (mt)",ylab="log(R/S)")
#	The expectation in a classic Stock-recruit curve would be points following along a line with a negative slope.
ricker <- lm(Y~S)
summary(ricker)
abline(reg=ricker)

#---------------------------------------------------------------------------------
# Step 3: fit the Kalman filter
nyr <- length(Y)
# initialize parameters for dlm
init.mean.a1 <-coef(ricker)[1]
init.b1 <-coef(ricker)[2]
Va0<-100000
Vb0<-100000

mod=function(para)
	dlm(FF=matrix(c(1,0),nc=2),
	JFF=matrix(c(0,1),nc=2),
	X=S,V=(exp(para[1]))^2,
	GG=diag(2),
	W=diag(c(exp(para[2])^2,0)),
	m0=c(init.mean.a1,init.b1),
	C0=diag(c(Va0,Vb0)))
  var_par=dlmMLE(Y,parm = log(c(0.7,0.7)+10^-7),build = mod, control=list(maxit=10000),method='BFGS',hessian=T)

mod.filt <- dlmFilter(Y,mod.fun(var_par$par))
mod.smooth <- dlmSmooth(Y,mod.fun(var_par$par))
#  plot the a(t) values and corresponding confidence intervals
plot(boca$year,dropFirst(mod.smooth$s[,1]),xlab="",ylab="")
title(main="Bocaccio",xlab="Year",ylab="Smoothed a value")
# add the confidence interval to this plot
# mod.filt$mod$V contains the estimated observation error variance (0.486)
mod.filt$mod$V
# mod.filt$mod$VW contains the estimated process error covariance matrix (0.071)
mod.filt$mod$W
# the estimated signal-to-noise ratio is 0.15, which is on the low side
mod.filt$mod$W[1,1]/mod.filt$mod$V

#--------------------------------------------------------------------------------
# Step 4: Plot the smoothed a value and 95% confidence intervals
# The variances are extracted in terms of their Singular Value Decompositions
se.smo=do.call(rbind,lapply(dlmSvd2var(mod.smooth$U.S, mod.smooth$D.S),function(cov.S)sqrt(diag(cov.S))))
a.smo <- mod.smooth$s[-1,1]
pl <- a.smo + qnorm(0.05, sd = se.smo[-1,1])
pu <- a.smo + qnorm(0.95, sd = se.smo[-1,1])
matplot(boca$year,cbind(pl,pu),type = 'l', lty=2, col=2, xlab="",ylab="")
lines(boca$year,a.smo,type = 'o', pch = 20, col = "blue")
title(main="Bocaccio",xlab="Year",ylab="Smoothed a value")


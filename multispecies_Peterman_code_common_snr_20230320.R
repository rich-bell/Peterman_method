#	20230320
#	time-varying productivity estimates with state-space Ricker model, common signal to noise ratio

#	Rich Bell, Adrien Tableau, Jeremy Collie, Coilin Minto
#	rich.bell@tnc.org


rm(list=ls())
library(fields)
library(dlm)
library(TeachingDemos)



# dlm ---------------------------------------------------------------------
load(file="west_coast_input_SSB_R_20230309.RData")

head(df)	#	species clipped to just years with rec dev
head(df2)	#	species clipped to just years with rec dev, but NA's added so every time series goes back to 1954


#	Time-Varying Ricker model

###########################################
Va0=10^5
Vb0=10^5

# model all stocks together, estimated parameters: sd_obs for each stock and signal to noise ratio for all stocks together

#	initial parameter estimates, one for each species plus the common signal to noise ratio
init.dat3.par<-readRDS("dat3_init_par_20230309.RDS")
init.snr<-log(0.4)

stocks<-unique(df$stock)
use.stock<-stocks[c(1:31)]

df3=split(df2[df2$stock %in% use.stock,],df2$stock[df2$stock %in% use.stock],drop=T)
ns=length(df3)
df4=lapply(df3,function(stock)stock[!is.na(stock$SSB_t),])			#	removes NA
b0=sapply(df4,function(stock)summary(lm(stock$log.R.S~stock$SSB_t))[[4]][2,1])
a0=mapply(function(stock,b0)mean(stock$log.R.S[1:5]-b0*stock$SSB_t[1:5]),df4,b0)


mod3=function(para){
  dlm(FF=cbind(diag(ns),matrix(0,nc=ns,nr=ns)),
      JFF=cbind(matrix(0,nc=ns,nr=ns),diag(1:ns)),
      X=do.call(cbind,lapply(df3,function(stock)stock$SSB_t)),
      V=diag((exp(para[1:ns]))^2,nr=ns),
      GG=diag(2*ns),
      W=diag(c((exp(para[ns+1])*exp(para[1:ns]))^2,rep(0,ns))),
      m0=c(a0,b0),
      C0=diag(c(rep(Va0,ns),rep(Vb0,ns))))}


#	Can take an hour to run.
dat3=dlmMLE(
  y=do.call(cbind,lapply(df3,function(stock)stock$log.R.S)),
  parm = c(init.dat3.par[c(use.stock)],init.snr),
  build = mod3, 
  control=list(maxit=100000,trace=1),
  #method='BFGS',hessian=T,debug=TRUE)	#	debug mode commented out.  in initial fitting, used different methods. 
  method='BFGS',hessian=T)

#	output of the model are the estimated parameters which are the stand dev of V (observation error) and the signal to noise ratio.  Divide those two to get stand dev of W, process error.  From those can then build the time varying Ricker model
#	three col are estimated sd of V, calculated sd of W, and estiamted SNR across all species
cbind(sqrt(diag(mod3(dat3$par)$V)),sqrt(diag(mod3(dat3$par)$W)[1:ns]),sqrt(diag(mod3(dat3$par)$W)[1:ns])/sqrt(diag(mod3(dat3$par)$V)))
exp(dat3$par[ns+1])	#	snr


smo3=dlmSmooth(y=do.call(cbind,lapply(df3,function(stock)stock$log.R.S)),mod = mod3(dat3$par))
#smo3$s	#	this sort of replaces smo3, because it is just the smoothed terms
fil3=as.data.frame(dlmFilter(y=do.call(cbind,lapply(df3,function(stock)stock$log.R.S)),mod = mod3(dat3$par))$m)
se.smo3=do.call(rbind,lapply(dlmSvd2var(smo3$U.S, smo3$D.S),function(cov.S)sqrt(diag(cov.S))))


#	Need full smo3 to get se.smo3, but then need it back in data.frame for other components
smo3<-as.data.frame(smo3$s)
yr<-c(unique(df3[[1]]$year), 2022)

fil3$year=smo3$year=yr
iyear=sapply(split(df,df$stock,drop=T),function(x)x$year[1])	

out_inits=cbind.data.frame(fi=sapply(1:length(df3),function(x){
  fil3[fil3$year==iyear[x],x]
}),
sm=sapply(1:length(df3),function(x){
  smo3[smo3$year==iyear[x],x]
}),a0=a0)		#	}),a0=selection$a0)

out_inits$CVfi=100*(out_inits$fi-out_inits$a0)^2/out_inits$fi
out_inits$CVsm=100*(out_inits$sm-out_inits$a0)^2/out_inits$sm

# time series outputs
out3=lapply(1:ns,function(s){
#  tmp=cbind.data.frame(df3[[s]][,c('year','stock','log.R.S')],smo3$s[-1,c(s,s+ns)],se.smo3[-1,c(s,s+ns)])
  tmp=cbind.data.frame(df3[[s]][,c('year','stock','log.R.S')],smo3[-1,c(s,s+ns)],se.smo3[-1,c(s,s+ns)])
  names(tmp)[4:7]=c('a','b','sd_a','sd_b')
  return(tmp)})


dev.new()
par(mfrow=c(5,7),mar=c(2,2,2,2),family='serif')
for(i in 1:length(use.stock)){
plot(out3[[i]]$year, out3[[i]]$a,typ='l',main=use.stock[i])
}	#	end for loop

par(mfrow=c(5,7),mar=c(2,2,2,2),family='serif')
for(i in 1:length(use.stock)){
plot(smo3[,i],typ='l',main=use.stock[i])
}	#	end for loop


# Time Invariant ricker model stock by stock

out4=lapply(df3,function(st){
  tmp=cbind.data.frame(unique(st$stock),data.frame(matrix(as.vector(summary(lm(st$log.R.S~st$SSB_t))$coefficients[,1:2]),nr=1)))
  names(tmp)=c('stock','a','b','sd_a','sd_b')
  return(tmp)
})

# same thing but with the b parameter from model 3 to make it easy to compare

out5=mapply(function(st,o3){
  b=o3$b[nrow(o3)]
  sd_b=o3$sd_b[nrow(o3)]
  da=st[!is.na(st$log.R.S),]
  a=mean(da$log.R.S-b*da$SSB_t)
  sd_a=sd(da$log.R.S-b*da$SSB_t)/sqrt(length(da$log.R.S))
  tmp=cbind.data.frame(unique(st$stock),data.frame(a,b,sd_a,sd_b))
  names(tmp)[1]='stock'
  return(tmp)
},df3,out3,SIMPLIFY = F)


###############################3



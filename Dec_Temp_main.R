## Codes used by the temperature team -------------------------------------

## ------------------------------------------------------------------------
## Working directory (needs to be adapted)
Tdir = "data/"
OUTdir="sim_data/"
plotDir="plots/"
#setwd("/crop_failure/")


## ------------------------------------------------------------------------
## Packages used

library(parallel)
library(ncdf4)

## ------------------------------------------------------------------------
## Data

source("readextranc.R") #netcdf file functions

args=(commandArgs(TRUE))
print(args)
if(length(args)>0){
    varname =args[1]
    Lsim =as.integer(args[2])
    mo.start =as.integer(args[3])
    day.start =as.integer(args[4])
    nsim =as.integer(args[5])
    yymin=as.integer(args[6])
    yymax=as.integer(args[7])
    fileanalo=args[8]
    jobid=args[9]
    alpha.TN <<- as.numeric(args[10])
    alpha.cal <<- as.numeric(args[11])
}else{ ## Default option
    varname="tx"
    Lsim = 31 #number of days of the simulation
    mo.start = 12 #month of beginning
    day.start =01 #day of beginning
    nsim =1000 #number of simulations
    yymin=1950 #first year
    yymax=2018 #last year
    fileanalo="/output_medium.txt" #analogue files
    jobid="test"
    alpha.TN <<- 0.75 #choice of the weight on temperature
    alpha.cal <<- 6 # Weight on seasonal cycle (how close do you want to be from the calendar day)
}

args[1] = varname
args[2] = Lsim
args[3] = mo.start
args[4] = day.start
args[5] = nsim
args[6] = yymin
args[7] = yymax
args[8] = fileanalo
args[9] = jobid
args[10] = alpha.TN
args[11] = alpha.cal

## Sets days in calendar year
## Creates l.mo, l.da et moda (list of days in year)
l.mo=c(1:12)
l.da=c(31,28,31,30,31,30,31,31,30,31,30,31)
moda=c()
for(i in l.mo){
  for(j in 1:l.da[i]){
    moda=c(moda,i*100+j)
  }
}
## ------------------------------------------------------------------------
## Read input data
## Read analog file
##example: fileanalo="http://www-lscedods.cea.fr/estimr2/DASE_NK/ana_slp_surface_base_rms_NA_latest_-80.0_50.0_22.5_70.0_1_30_20.txt"
analo = read.table(fileanalo,header=TRUE)
date.a = analo$date

date.a.cal=match(as.integer(substr(analo$date,5,8)),moda)

## Read temperature  data for the region [45.5°N,51.5°N,-1.5°E,8°E] - EOBS dataset
filin = paste(Tdir,varname,"_FR_mean.nc",sep="")
nc = nc_open(filin)
varnc=nc$var[[varname]]
## Treatment of time
nctime = nc$dim[['time']]
time=nctime$vals
conv.time=caldat(time+julday(1,1,1950))
years=conv.time$year
months=conv.time$month
days=conv.time$day
TN_day=ncvar_get(nc,varid=varnc)
nc_close(nc)
Date=10000*years+100*months+days

TN=data.frame(Date,TN_day)


## ------------------------------------------------------------------------
## Computer parameters for the parallel calculation
## Please adapt to your own cluster
ncpus = detectCores()
print(paste(detectCores(),"cores detected"))
print(paste("Calcul sur",ncpus,"CPUs"))

## ------------------------------------------------------------------------
## Stochastic simulation with weights on distance to calendar day
## and weights on temperature rank.
## lanamax is the number of days around the target day from which you can pick analogues
"simu.extrHW1" = function(t.start=19500101,Lsim=90,alpha.cal = 0.5,
  alpha.TN = 0.1,lanamax=30)
  {
  I0 = which(analo$date == t.start)
  t0.cal = match(as.integer(substr(t.start,5,8)),moda)
  t0=t.start
  T.sim=c(TN[TN$Date == t0,2])
  t.sim=c(t0)
  ndum=c()
  for(i in 1:Lsim){
    I = which(analo$date == t0)
    I1 = I +1
    t1=analo$date[I1]
    ana.d1 = c(analo$date[I0+i],unlist(analo[I1,2:21]))
    ana.d1 = intersect(ana.d1,TN$Date)
    ## Weight on calendar day to respect the seasonal cycle
    ana.d1.cal = match(as.integer(substr(ana.d1,5,8)),moda)
    diff.ana.cal = abs(ana.d1.cal - (t0.cal+i))
    weight.cal = exp(-alpha.cal*diff.ana.cal/lanamax)
    sdum=sum(weight.cal,na.rm=TRUE)
    weight.cal = weight.cal/sdum
    ## Weight on temperature (to get warmer analogues)
    d1=ana.d1
    Idum = match(ana.d1,TN$Date)
    TN.d1 = TN[[2]][Idum] ## Temperatures analogues
    TN.d1.sort=sort(TN.d1,index.return=TRUE,decreasing=TRUE,
                    na.last=TRUE,method="radix")
    weight.TN = exp(-alpha.TN*c(1:length(TN.d1)))
    ## Number of analogues = length(TN.d1)
    sdum=sum(weight.TN,na.rm=TRUE)
    weight.TN[TN.d1.sort$ix] = weight.TN/sdum
    ## Product of weights for probabilities normalization
    weight.all=weight.cal * weight.TN
    weight.all[is.na(weight.all)]=0
    weight.all = weight.all/sum(weight.all,na.rm=TRUE)
    ndum=c(ndum,length(which(weight.all >= 0.05)))
    ## Random sampling
    d1.max = sample(d1,size=1,prob=weight.all)
    T.sim = c(T.sim,TN[TN$Date == d1.max,2])
    t.sim = c(t.sim, d1.max)
    t0=ifelse(idyn0,d1.max,t1)
  }
  return(cbind(t.sim,T.sim))
}
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
## Definition of the simulation function
fun.name = "simu.extrHW1"
print(paste("Applying",fun.name))
SIMU.FUNC = match.fun(fun.name)

## ------------------------------------------------------------------------
## Wrapper for stochastic simulations
"wrap.extrHW" = function(k)
  {
      XX = SIMU.FUNC(t.start=t.start,Lsim=Lsim,
                       alpha.TN=alpha.TN,alpha.cal=alpha.cal)
    return(XX)
  }
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Wrapper for mean temperatures calculations
"wrap.mean" = function(i)
  {
    mm = mean(Xdum[[i]][,2])
    return(mm)
  }
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Simulations for the reanalysis
"simu.yy" = function(idyn=TRUE,yymin=1950,yymax=2018,
                     mo0=mo.start,day0=day.start)
{
  idyn0 <<- idyn
  l.X.mean.dyn=list()
  l.T.mean.dyn=list()
  l.X.dyn = list()
  for(yy in yymin:yymax){
    print(paste("Processing",yy))
    t.start <<- (yy*100+mo0)*100+day0
    Xdum <<- mclapply(seq(1,nsim,by=1),wrap.extrHW,mc.cores=ncpus)
    l.X.dyn[[as.character(yy)]]=Xdum
    X.mean = mclapply(seq(1,nsim,by=1),wrap.mean,mc.cores=ncpus)
    l.X.mean.dyn[[as.character(yy)]]=unlist(X.mean)
    iTN=which(TN$Date==t.start)
    l.T.mean.dyn[[as.character(yy)]]=mean(TN[c(iTN:(iTN+Lsim)),2],
                                          na.rm=TRUE)
    }
    return(list(l.X.mean=l.X.mean.dyn,
                l.T.mean=l.T.mean.dyn,
                l.X=l.X.dyn,ymin=yymin,ymax=yymax))
}
## ------------------------------------------------------------------------
ncpus=1
simu.dyn=simu.yy(idyn=TRUE,yymin=yymin,yymax=yymax)
simu.sta=simu.yy(idyn=FALSE,yymin=yymin,yymax=yymax)

setwd(OUTdir)
fname=paste(varname,"-m",mo.start,"d",day.start,"_cal",alpha.cal,"_temp_",alpha.TN,"-",jobid,".Rdat",sep="")

save(file=fname,simu.dyn,alpha.cal,alpha.TN,simu.sta,args)

## ------------------------------------------------------------------------
## Compute the number of days between 0 and 10 degrees
## in December
setwd(OUTdir)
list.in = system("ls tx-m12d1_*-test.Rdat",intern=TRUE)

## Data input
i=1
l.alpha=c()
Sum_dyn=list()
Sum_stat=list()
for(filin in list.in){
  load(filin)
  l.alpha[[i]]=alpha.TN
  sum_dyn <- rbind()
  sum_stat=rbind()
  for(k in 1:length(unique(years))){
    dataYearDyn <- simu.dyn$l.X[[k]]
    dataYearSta <- simu.sta$l.X[[k]]
    tempSta <- c()
    tempDyn <- c()
    for(j in 1:10){
      dataDyn <- dataYearDyn[[j]]
      dataSta <- dataYearSta[[j]]
      tempDyn <- c(tempDyn, sum(dataDyn>0 & dataDyn<10))
      tempSta <- c(tempSta, sum(dataSta>0 & dataSta<10))
    }
    sum_dyn <- rbind(sum_dyn,tempDyn)
    sum_stat <- rbind(sum_stat,tempSta)
  }
  Sum_dyn[[i]] = sum_dyn
  Sum_stat[[i]] = sum_stat
  i=i+1
}
print(i)

## number of days for the observation set
sum_obs <- c()
for(i in unique(years)){
  data <- TN$TN_day[years == i & months == 12]
  sum_obs <- c(sum_obs,sum(data < 10 & data > 0))
}

## ------------------------------------------------------------------------
## Fit of a Beta-binomial distribution to the number
## of days between 0 and 10 degrees  in December

## Estimation of the parameters of the beta binomial distribution
n <- 31
m1 <- mean(sum_obs)
m2 <- mean((sum_obs)^2)

alpha <- (n*m1-m2)/(n*(m2/m1 - m1 -1) + m1)
beta <- (n-m1)*(n-m2/m1)/(n*(m2/m1-m1-1)+m1)

library('rmutil')
# rmutil provides the inverse CDF of beta binomial but has different parameters
s=alpha+beta
m=alpha/(alpha+beta)

l.probs=c(0.5,0.3,0.1,0.01,0.001,0.0001)
Ndum=qbetabinom(l.probs,size=n,m,s)

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## Codes of the figures of Pfleiderer et al.

## ------------------------------------------------------------------------
## Figure 4 of Pfleiderer et al.

## ------------------------------------------------------------------------
## Figure A1 of Pfleiderer et al.

setwd(plotDir)
list.in = list.files(Tdir)
list.in = paste0(Tdir,'/',list.in)
i=1
l.alpha=c(0.4,0.6)
Tobs=c()
Tsimdyn=list()
Tsimsta=list()
for(filin in list.in){
  load(filin)
  Tobs=unlist(simu.dyn$l.T.mean)
  Tsimdyn[[i]]=unlist(simu.dyn$l.X.mean)
  Tsimsta[[i]]=unlist(simu.sta$l.X.mean)
  i=i+1
}
sd.obs=sd(Tobs)
m.obs=mean(Tobs)
i0=which(l.alpha==0)
sd.dyn0=sd(Tsimdyn[[i0]])
m.dyn0=mean(Tsimdyn[[i0]])

proba.RTdyn=c()
proba.RTsta=c()
for(i in 1:length(l.alpha)){
  dum=1-pnorm(mean(Tsimdyn[[i]]),mean=m.obs,sd=sd.obs)
  proba.RTdyn=c(proba.RTdyn,dum)
  dum=1-pnorm(mean(Tsimsta[[i]]),mean=m.obs,sd=sd.obs)
  proba.RTsta=c(proba.RTsta,dum)
}

l.proba=1-pnorm(c(20:24),mean=m.obs,sd=sd.obs)
l.probs=c(0.5,0.7,0.9,0.99,0.99999)
Tdum=qnorm(l.probs,mean=m.obs,sd=sd.obs)

filout=paste("Tg_",args[1],"-m",args[4],".pdf",sep="")
pdf(filout)
default_marges=c(5.1,4.1,4.1,2.1)
par(mar=default_marges + c(0, 1, -1, 0),cex.axis = 1.5,cex.lab  = 1.5)
plot(c(0,0.8),c(0,18),type="n",xlab="",ylab="December Temperature (?C)",xlim=(c(0,1)),ylim=(c(0,18)),axes=FALSE,cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
abline(h=Tobs[["2015"]],lty=3)
boxplot(Tobs,add=TRUE,boxwex=0.2,at=0.2)
for(i in 1:length(l.alpha)){
  boxplot(Tsimdyn[[i]],at=l.alpha[i],axes=FALSE,add=TRUE,
          col="red",boxwex=0.2)
  boxplot(Tsimsta[[i]],at=l.alpha[i],axes=FALSE,add=TRUE,
          col="blue",boxwex=0.2)
}
labels_list=c('Observations','Z500','SLP')
axis(side=1,at=c(1,2,3)/5,labels=labels_list,cex.axis=1.2)
box()
dev.off()

## ------------------------------------------------------------------------
## Figure A2 of Pfleiderer et al.

## run HWgen_AnlmSa.R
idynO=TRUE
alphas.cal <- c(0,.1,.2,.5,.75,1,2,3,4,5,6,7,8,9,10)
score=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
k=1
for(qqi in alphas.cal){
  print(qqi)
  a=simu.extrHW1(19500101,1000,qqi,0.75,30)
  b=(trunc(a/100)-trunc(a/10000)*100)[,1]
  score[k]=(sum(b==12)+sum(b==11)+sum(b==1)+sum(b==2)-1)/1000
  k=k+1
}

pdf("alpha_calendar.pdf")
default_marges=c(5.1,4.1,4.1,2.1)
par(mar=default_marges + c(0, 1,-2, 0),cex.axis = 1.5,cex.lab  = 1.5,lwd=2)
plot(c(0,10),c(0,80),type="n",xlab="alpha calendar",ylab="% of winter days",axes=FALSE)
points(alphas.cal, score*100,type="b",pch=5)
abline(v=6,lty='dotted',col='red')

axis(side=1)
axis(side=2)
box()
dev.off()

## ------------------------------------------------------------------------
## Figure A3 of Pfleiderer et al.

png(filename="BoxPlot_Numer_Days.png",width=1500,height=1500,units = "px", pointsize = 40, bg = "white", res = NA)
default_marges=c(5.1,4.1,4.1,2.1)
par(mar=default_marges + c(0, 1, -1, 0),cex.axis = 1.5,cex.lab  = 1.5)
plot(c(-0.2,1.1),c(0,32),type="n",xlab="alpha",ylab="Number of days",axes=FALSE)
abline(h=sum_obs[66],lty=3,lwd=4) #66 correspond to year 2015

boxplot(sum_obs,add=TRUE,boxwex=0.2,at=-0.1,lwd=4)
for(i in 1:length(l.alpha)){
  boxplot(c(Sum_dyn[[i]]),at=l.alpha[i],lwd=4,axes=FALSE,add=TRUE,col="red",boxwex=0.2)
  boxplot(c(Sum_stat[[i]]),at=l.alpha[i],lwd=4,axes=FALSE,add=TRUE,col="blue",boxwex=0.2)
}

axis(side=1,at=l.alpha,lwd=4)
axis(side=2,lwd=4)
axis(side=4,at=Ndum,labels=l.probs,lwd=0.5)
dev.off()




## END of SCRIPT
## ------------------------------------------------------------------------

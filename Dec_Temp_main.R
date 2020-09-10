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



## END of SCRIPT
## ------------------------------------------------------------------------

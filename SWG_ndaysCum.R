## An extreme heatwave stochastic weather generator based on analogues
## and a simplified "importance sampling" algorithm
## Pascal Yiou (LSCE), June 2016, Feb 2018, June 2019
## This code is distributed freely and as is under a CeCill License.
## here some adaptations are added see Pfleiderer et al. ???

base_dir="/Users/peterpfleiderer/Projects/analogues/" ## Needs to be adapted
try(setwd(base_dir), silent=TRUE)

#base_dir = "C:/Users/jason/OneDrive/Documents/GitHub/xtrem_precip_analogues/" ## Needs to be adapted
#try(setwd(base_dir), silent=TRUE)

# test

source("scripts/readextranc.R") #netcdf file functions
library(parallel)
library(ncdf4)

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
  alpha.var <<- as.numeric(args[10])
  alpha.cal <<- as.numeric(args[11])
  nn <<- as.integer(args[12])
}else{ ## Default option
  varname="rr"
  Lsim = 30 #number of days of the simulation
  mo.start = 05 #month of beginning
  day.start = 01 #day of beginning
  nsim = 1 #number of simulations
  yymin=2016 #first year
  yymax=2016 #last year
  fileanalo="Data/analogues_slp_-50_30_30_70.txt" #analogue files
  jobid="springPrecip_v1"
  alpha.var <<- 0.5 #choice of the weight on temperature
  alpha.cal <<- 0.1 # Weight on seasonal cycle (how close do you want to be from the calendar day)
  nn <<- 3
}


args_tmp = c()
args_tmp[1] = varname
args_tmp[2] = Lsim
args_tmp[3] = mo.start
args_tmp[4] = day.start
args_tmp[5] = nsim
args_tmp[6] = yymin
args_tmp[7] = yymax
args_tmp[8] = fileanalo
args_tmp[9] = jobid
args_tmp[10] = alpha.var
args_tmp[11] = alpha.cal
args_tmp[12] = nn

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
fileanalo = 'Data/analogues.txt'

analo = read.table(fileanalo,header=TRUE)
date.a = analo$date

date.a.cal=match(as.integer(substr(analo$date,5,8)),moda)

## Read temperature  data for the region [45.5°N,51.5°N,-1.5°E,8°E] - EOBS dataset
filin = paste(base_dir,'Data/',varname,"_FR_mean.nc",sep="")
nc = nc_open(filin)
varnc=nc$var[[varname]]
## Treatment of time
nctime = nc$dim[['time']]
time=nctime$vals
conv.time=caldat(time+julday(1,1,1950))
years=conv.time$year
months=conv.time$month
days=conv.time$day
varIN_day=ncvar_get(nc,varid=varnc)
nc_close(nc)
Date=10000*years+100*months+days

varIN=data.frame(Date,varIN_day)



## ------------------------------------------------------------------------
## Computer parameters for the parallel calculation
## Please adapt to your own cluster

ncpus = detectCores()
print(paste(detectCores(),"cores detected"))

ncpus = 2
print(paste("Calcul sur",ncpus,"CPUs"))

## ------------------------------------------------------------------------
## Stochastic simulation with weights on distance to calendar day
## and weights on temperature rank.
## lanamax is the number of days around the target day from which you can pick analogues
"simu.extrHW1" = function(t.start=20030601,Lsim=90,alpha.cal = 0.5,
                          alpha.var = 0.1,lanamax=30, nn=3)
{
  I0 = which(analo$date == t.start)
  t0.cal = match(as.integer(substr(t.start,5,8)),moda)
  t0=t.start
  varIN.sim=c(varIN[varIN$Date == t0,2])
  t.sim=c(t0)
  Prob.sim = c(1)
  ndum=c()
  for(i in 1:as.integer(Lsim/nn)){

    I = which(analo$date == t0)
    I1 = I + 1
    t1 = analo$date[I1]

    #print(paste(c(t0,t1),sep=' - '))
    ana.d1 = c(analo$date[I0+1],unlist(analo[I1,2:21]))
    ana.d1 = intersect(ana.d1,varIN$Date)
    ## Weight on calendar day to respect the seasonal cycle
    ana.d1.cal = match(as.integer(substr(ana.d1,5,8)),moda)
    diff.ana.cal = abs(ana.d1.cal - (t0.cal+i))
    weight.cal = exp(-alpha.cal*diff.ana.cal/lanamax)
    sdum=sum(weight.cal,na.rm=TRUE)
    weight.cal = weight.cal/sdum
    ## Weight on temperature (to get warmer analogues)
    d1=ana.d1
    Idum = match(ana.d1,varIN$Date)
    varIN.d1 = varIN[[2]][Idum] ## Temperatures analogues
    for (j in 1:(nn-1)){
      varIN.d1 = varIN.d1 + varIN[[2]][Idum + j]
    }
    varIN.d1.sort=sort(varIN.d1,index.return=TRUE,decreasing=TRUE,
                    na.last=TRUE,method="radix")
    weight.varIN = exp(-alpha.var*c(1:length(varIN.d1)))
    ## Number of analogues = length(varIN.d1)
    sdum=sum(weight.varIN,na.rm=TRUE)
    weight.varIN[varIN.d1.sort$ix] = weight.varIN/sdum
    ## Product of weights for probabilities normalization
    weight.all=weight.cal * weight.varIN
    weight.all[is.na(weight.all)]=0
    weight.all = weight.all/sum(weight.all,na.rm=TRUE)
    ndum=c(ndum,length(which(weight.all >= 0.05)))
    ## Random sampling
    id.max = sample(1:length(d1),size=1,prob=weight.all)
    Prob.sim = c(Prob.sim, (1/length(weight.all)) * (1/weight.all[id.max]), as.integer(1:(nn-1) /100) +1)
    #print(paste(c(Prob.sim,1/length(weight.all),weight.all[id.max]),sep=' - '))
    d1.max = d1[id.max]
    date.mod = which(varIN$Date == d1.max)
    for (j in 0:(nn-1)){
      t.sim = c(t.sim, varIN$Date[date.mod + j])
      varIN.sim = c(varIN.sim,varIN[date.mod + j,2])
    }
    #print(paste(c(i,t0,t1,varIN$Date[date.mod]), sep=' - '))
    t0=ifelse(idyn0,varIN$Date[date.mod + j],t1)
  }
  return(cbind(t.sim,varIN.sim,Prob.sim))
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
                       alpha.var=alpha.var,alpha.cal=alpha.cal,nn=nn)
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

"wrap.sum" = function(i)
  {
  mm = sum(Xdum[[i]][,2])
  return(mm)
}

"wrap.prod" = function(i)
{
  mm = prod(Xdum[[i]][,3])
  return(mm)
}
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Simulations for the reanalysis
"simu.yy" = function(idyn=TRUE,yymin=1948,yymax=2018,
                     mo0=mo.start,day0=day.start)
{
  idyn0 <<- idyn
  l.X.mean.dyn=list()
  l.X.sum.dyn=list()
  l.varIN.mean.dyn=list()
  l.varIN.sum.dyn=list()
  l.Prob.dyn=list()
  l.X.dyn = list()
  for(yy in yymin:yymax){
    print(paste("Processing",yy))
    t.start <<- (yy*100+mo0)*100+day0
    Xdum <<- mclapply(seq(1,nsim,by=1),wrap.extrHW,mc.cores=ncpus)
    l.X.dyn[[as.character(yy)]]=Xdum
    # mean
    X.mean = mclapply(seq(1,nsim,by=1),wrap.mean,mc.cores=ncpus)
    l.X.mean.dyn[[as.character(yy)]]=unlist(X.mean)
    # sum
    X.sum = mclapply(seq(1,nsim,by=1),wrap.mean,mc.cores=ncpus)
    l.X.sum.dyn[[as.character(yy)]]=unlist(X.sum)
    # prob
    Prob = mclapply(seq(1,nsim,by=1),wrap.prod,mc.cores=ncpus)
    l.Prob.dyn[[as.character(yy)]]=unlist(Prob)

    ivarIN=which(varIN$Date==t.start)
    l.varIN.mean.dyn[[as.character(yy)]]=mean(varIN[c(ivarIN:(ivarIN+Lsim)),2],na.rm=TRUE)
    l.varIN.sum.dyn[[as.character(yy)]]=sum(varIN[c(ivarIN:(ivarIN+Lsim)),2],na.rm=TRUE)
  }
    #return(l.X.dyn)
    return(list(l.X.mean=l.X.mean.dyn, l.X.sum=l.X.sum.dyn, l.varIN.mean=l.varIN.mean.dyn,l.varIN.sum.dyn, l.X=l.X.dyn,ymin=yymin,ymax=yymax, Prob.sim=l.Prob.dyn))
}
## ------------------------------------------------------------------------

simu.dyn=simu.yy(idyn=TRUE,yymin=yymin,yymax=yymax)
simu.sta=simu.yy(idyn=FALSE,yymin=yymin,yymax=yymax)

fname=paste("sim_data/",varname,"-m",mo.start,"d",day.start,"nDay",Lsim,"_cal",alpha.cal,"_precip_",alpha.var,"_N",nsim,"-aggregNcum",nn,'-',jobid,".Rdat",sep="")

save(file=fname,simu.dyn,alpha.cal,alpha.var,simu.sta,args_tmp)

q("no")
## END of SCRIPT
## ------------------------------------------------------------------------

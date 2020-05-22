#Trace une carte "longitude-latitude" d'un champ, avec la ligne de continents
# qui correspond
"image.cont" <-
function(lonF,latF,champ,titre="",Ichange=numeric(0),
         zlev=seq(min(champ,na.rm=TRUE),max(champ,na.rm=TRUE),length=11),
         transpose=TRUE,mar=c(5,5,5,6),legend=TRUE,xlab="Longitude",
         ylab="Latitude",satur=FALSE,paquet="fields",add=FALSE,
         lonrange=range(lonF),latrange=range(latF))
{
#library(paquet)
                                        #zlev=seq(-0.85,0.6,length=11)
col10=rainbow(length(zlev)-1,start=0,end=2/6)
#par( mar=c(10,5,5,5))
par( mar=mar)
if(satur){
champ[champ<=min(zlev)]=min(zlev)
champ[champ>=max(zlev)]=max(zlev)
}
if(length(Ichange)>0) champ[Ichange]=1
if (transpose)
  dum=t(matrix(champ,length(latF),length(lonF)))
else
  dum=matrix(champ,length(lonF),length(latF))
latF.sort=sort(latF,index.return=TRUE)
lonF.sort=sort(lonF,index.return=TRUE)
plot(lonrange,latrange,type="n",xlab=xlab,ylab=ylab,xlim=lonrange,
     ylim=latrange)
image(sort(lonF),sort(latF),dum[lonF.sort$ix,latF.sort$ix],
      col=col10[length(col10):1],
       xlab=xlab,ylab=ylab,main=titre,breaks=zlev,add=TRUE)
if(paquet=="fields"){
  library(fields)
  world(add=TRUE)
#  world(xlim=range(lonF),ylim=range(latF),add=TRUE)
}
if(paquet=="maps"){
  library(maps)
  map(add=TRUE)
}
if(legend) image.plot(dum[,length(latF):1],col=col10[length(col10):1],
       legend.only=TRUE,zlim=range(zlev))
}

#Trace une carte "longitude-latitude" de contours d'un champ, avec la ligne
#de continents qui correspond
"image.cont.c" <-
function(lonF,latF,champ,titre="",Ichange=numeric(0),
         zlev=pretty(champ,10),
#         zlev=seq(min(champ,na.rm=TRUE),max(champ,na.rm=TRUE),length=11),
         transpose=TRUE,mar=c(5,5,5,6),xlab="Longitude",
         ylab="Latitude",col="blue",add=FALSE,paquet="fields",lty=1)
{
#zlev=seq(-0.85,0.6,length=11)
col10=rainbow(length(zlev)-1,start=0,end=2/6)
#par( mar=c(10,5,5,5))
par( mar=mar)
if(length(Ichange)>0) champ[Ichange]=1
if (transpose)
  dum=t(matrix(champ,length(latF),length(lonF)))
else
  dum=matrix(champ,length(lonF),length(latF))
latF.sort=sort(latF,index.return=TRUE)
nlev=length(zlev)
contour(lonF,sort(latF),dum[,latF.sort$ix],
       xlab=xlab,ylab=ylab,main=titre,col=col,add=add,nlevels=nlev,
      levels=zlev,lty=lty)
if(paquet=="fields"){
  library(fields)
  world(xlim=range(lonF),ylim=range(latF),add=TRUE)}
if(paquet=="maps"){
  library(maps)
  map(add=TRUE)
}
#world(xlim=range(lonF),ylim=range(latF),add=TRUE)
}

# Fonction generale pour tracer une carte avec des points de couleur
trace.po <- function(titre="",xcontrib,lon=lon,lat=lat,legende=TRUE,
         qq=NULL,col=NULL,lonlim=c(-10,30),latlim=c(35,60),
         mar=c(3,3,3,8))
{
#layout(matrix(1:2,1,2),widths=c(4,1))
library(maps)
xmoy=xcontrib
if(is.null(col)) col6=rainbow(6,start=4/6,end=0)
else col6=col
if(is.null(qq)) qq=quantile(xmoy,seq(0,1,length=length(col6)),na.rm=TRUE)
l.col=c()
for(i in 1:length(xmoy)) l.col=c(l.col,which(xmoy[i]>=qq[1:(length(qq)-1)] &
      xmoy[i]<qq[2:length(qq)]))
par(mar=mar)
plot(lonlim,latlim,type="n",xlab="",ylab="",main=titre)
#world(add=TRUE)
map(add=TRUE)
points(lon,lat,cex=3,col="black",bg=col6[l.col],pch=21)
if (legende){
  image.plot(xmoy,col=col6,legend.only=TRUE,zlim=range(qq),legend.mar=8)
}
}
# Fin de la fonction trace.po

# Fonction generale pour tracer une carte avec des points de couleur
# sur des valeurs discretes
trace.po2 <- function(titre="",Xm,lon=lon,lat=lat,ncolor=5,col=NULL)
{
library(maps)
xmoy=Xm
if(is.null(col))col6=rainbow(ncolor,start=4/6,end=0)
else col6=col
plot(seq(-10,30,length=10),seq(35,60,length=10),type="n",
       xlab="Longitude",ylab="Latitude",main=titre)
#world(add=TRUE)
map(add=TRUE)
points(lon,lat,cex=2,col=col6[xmoy],bg=col6[xmoy],pch=21)
}
# Fin de la fonction trace.po2

#Trace une carte d'un champ discret, avec ligne de continents
"image.cont.disc" <-
function(lonF,latF,champ,titre="",col=c(1:4),pref="",leg=NULL)
{
opar=par()
nf=layout(matrix(1:2,1,2),c(4,1),TRUE)
par(mar=c(5,4,4,0))
dum=t(matrix(champ,length(latF),length(lonF)))
image(lonF,latF,dum,col=col,xlab="Longitude",ylab="Latitude",main=titre)
world(xlim=range(lonF),ylim=range(latF),add=TRUE)
par(mar=c(5,0,4,1))
plot(c(0:1),c(0:1),type="n",axes=FALSE,xlab="",ylab="")
if(is.null(leg)) leg=paste(pref,c(1:length(col)))
legend(0,0.5,leg,col4)
layout(matrix(1,1,1))
}

# Trace de plusieurs cartes (mchamp)
"image.mfield" <-
function(lonF,latF,mchamp,mchamp2=NULL,titre="",
         nplot=4,rng=range(mchamp[,1:nplot]))
{
col10=rainbow(10,start=0,end=2/6)
col10=col10[10:1]
dimchamp=dim(mchamp)
nreg=nplot
n1=ceiling(nplot/2)
layout(matrix(1:(2*n1),2,n1))
for (i in 1:nreg){
  if(rng==-1){
    breaks=seq(min(mchamp[,i]),max(mchamp[,i]),length=11)
  }
  else{
    breaks=seq(rng[1],rng[2],length=11)
  }
  image(lonF,latF,t(matrix(mchamp[,i],length(latF),length(lonF))),
        breaks=breaks,
        xlab="Longitude",ylab="Latitude",main=paste(titre,i),col=col10)
  if(!is.null(mchamp2))
    contour(lonF,latF,t(matrix(mchamp2[,i],length(latF),length(lonF))),
            add=TRUE)
  world(xlim=range(lonF),ylim=range(latF),add=TRUE)
}
layout(matrix(1,1,1))
}

"map_circle" <-function(Xdat,IF,infoKNMI)
{
  plot(seq(-10,30,length=10),seq(35,60,length=10),type="n",
       xlab="Longitude",ylab="Latitude",main="")
  world(add=TRUE)
  for(i in 1:length(IF)){
    lon.sta=as.numeric(strsplit(as.character(infoKNMI$LON[IF[i]]),":")[[1]])
    lat.sta=as.numeric(strsplit(as.character(infoKNMI$LAT[IF[i]]),":")[[1]])
    if(dat.SEA.can[ireg+(i-1)*5+1]-dat.SEA.tot[ireg+(i-1)*5+1]>0) col="blue"
    else col="red"
    cex=50*abs(dat.SEA.can[ireg+(i-1)*5+1]-dat.SEA.tot[ireg+(i-1)*5+1])
    points(lon.sta[1]+0.01*lon.sta[2],lat.sta[1]+0.01*lat.sta[2],
        cex=cex,col=col)
  }
  legend(-10,50,legend=paste(c(2,4,6),"%"),pch=21,pt.cex=0.5*c(2,4,6))
  title(paste("Nb Precip days % diff, months",mon.1,"-",mon.2,", reg",ireg))
}


###Trace une carte "longitude-latitude" de contours d'un champ d'anomalies,
## (donc centre en 0) avec la ligne de continents qui correspond
## Les isolignes <0 sont en pointilles et le 0 en gras
## P. Yiou, Oct. 2019.
"image.cont.c.ano" <-
  function(lonF,latF,champ,titre="",Ichange=numeric(0),
           zlev=pretty(champ,10),
           transpose=TRUE,mar=c(5,5,5,6),xlab="Longitude",
           ylab="Latitude",colo="blue",add=FALSE,paquet="maps",lty=1)
  {
    par( mar=mar)
    if(length(Ichange)>0) champ[Ichange]=1
    if (transpose){
      dum=t(matrix(champ,length(latF),length(lonF)))
    }else{
      dum=matrix(champ,length(lonF),length(latF))
    }
    latF.sort=sort(latF,index.return=TRUE)
    nlev=length(zlev)
    llty=rep(1,nlev)
    llty[zlev<0]=2
    llwd=rep(1,nlev)
    llwd[zlev==0]=2
    colo=rep(1,nlev)
    colo[zlev<0]='#189cad'
    colo[zlev>0]='#7418ad'
    colo[zlev==0]='#060e69'
    contour(lonF,sort(latF),dum[,latF.sort$ix],
            xlab=xlab,ylab=ylab,main=titre,col=colo,add=add,nlevels=nlev,
            levels=zlev,lty=llty,lwd=llwd)
    if(paquet=="fields"){
      library(fields)
      world(xlim=range(lonF),ylim=range(latF),add=TRUE, col='red')}
    if(paquet=="maps"){
      library(maps)
      map(col='#3d3d3d', add=TRUE)
    }
  } 

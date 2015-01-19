library(clim.pact)

now <- function(format="dd.mm.yyyy") {
  cmon<-c("Jan","Feb","Mar","Apr","May","Jun",
          "Jul","Aug","Sep","Oct","Nov","Dec")
  imon <- 1:12
  now <- date()
  day.now<-as.numeric(substr(now,9,10))
  mon.now <- substr(now,5,7)
  mon.now <- imon[is.element(cmon,mon.now)]
  year.now <- as.integer(substr(now,21,24))
  if (format=="dd.mm.yyyy") now <- paste(day.now,mon.now,year.now,sep=".")
  if (format=="dd.mon.yyyy") now <- paste(day.now,substr(now,5,7),year.now)
  if (format=="c(day,mon,year)")  now <- c(day.now,mon.now,year.now)
  now
}



###########################################################################

#Innlasting av månedsverdier:

KDVH4DS <-
function(StNr=18700,fom="01.01.1950",tom="now",param=c("TAM","RR"),month=TRUE,
         method="mean",verbose=FALSE,mean.from.daily=FALSE,min.data.points=20,whole=TRUE){

  if (length(param)==1 & param[1]=="TAM" & (!month)) param <- c("TAM","RR")  # For daily data objects
  if (month) param <- param[1]
  data(dvh.station.list,envir=environment())
#  dnmi.meta <- read.table("data/dvh.station.list",header=TRUE,as.is=TRUE)
  if (is.character(StNr)) {
    imatch <- is.element(lower.case(substr(dnmi.meta$Navn,1,nchar(StNr))),lower.case(StNr))
    if (sum(imatch)==1) StNr<- dnmi.meta$Stnr[imatch]
    if (sum(imatch)>1) {
      print(dnmi.meta$Navn[imatch])
      iloc <- as.numeric(readline(paste("Which of these ( 1 -",sum(imatch),")? ")))
      StNr<- dnmi.meta$Stnr[imatch][iloc]
    }
    if (sum(imatch)==0) {
      print(dnmi.meta$Navn)
      stop("Could not find requested station")
    }
  }
  alt<- as.numeric(dnmi.meta$Hoh[dnmi.meta$Stnr==StNr][1])
  lon<- dnmi.meta$Lon[dnmi.meta$Stnr==StNr][1]
  lat<- dnmi.meta$Lat[dnmi.meta$Stnr==StNr][1]
  location <- dnmi.meta$Navn[dnmi.meta$Stnr==StNr][1]

  if (length(param)==1) {
    ele <- switch(param,
                  'TAM'=101,'TAX'=112,'TAN'=122,'RR'=601,'SA'=711)
    unit <- switch(param,
                  'TAM'='deg C','TAX'='deg C','TAN'='deg C','RR'='mm/month')
  } else {
    ele <- rep(NA,2); unit <- rep(NA,2)
    for (i in 1:2) {
      ele[i] <- switch(param[i],
                    'TAM'=101,'TAX'=112,'TAN'=122,'RR'=601,'PRN'=401)
      unit[i] <- switch(param[i],
                    'TAM'='deg C','TAX'='deg C','TAN'='deg C','RR'='mm/month','PRN'='hPa')
      }
  } 
  fom.c <- c(as.numeric(substr(fom,1,2)),as.numeric(substr(fom,4,5)),
             as.numeric(substr(fom,7,10)))
  if (tom=="now") tom <- now()
  t1<-now(format="c(day,mon,year)")    
  e.rows.mon <- t1[2] + 12*t1[3] - (fom.c[2]+12*fom.c[3])
  e.rows.day <- julday(t1[1],t1[2],t1[3]) - julday(fom.c[1],fom.c[2],fom.c[3])
  
  if ((month) & (method=="mean") & (!mean.from.daily)) {
    Filnavn <- "http://typhoon.oslo.dnmi.no/metnopub/production/metno?re=15&ct=text/plain&del=space&ddel=dot&nod=NA&split=1"
  } else
    Filnavn <- "http://typhoon.oslo.dnmi.no/metnopub/production/metno?re=14&ct=text/plain&del=space&ddel=dot&nod=NA&split=1"
  
  for (i in 1:length(param)) Filnavn <- paste(Filnavn,"&p=",param[i],sep="")
  Filnavn <- paste(Filnavn,"&fd=",fom,"&td=",tom,"&s=",sep="")
  Filnavn <- paste(Filnavn,StNr,sep="")

  if (verbose) print(Filnavn)

  if ((month) & (method=="mean")) {
    #print("read data from URL")
    Datasett <- read.table(Filnavn,header=TRUE)
    #print(summary(Datasett))
    #print("Extract dates")
    if (!whole) selected <- rep(TRUE,length( Datasett$Year)) else
                selected <- (Datasett$Year + (Datasett$Month-0.5)/12) < (t1[3] + (t1[2]-0.5)/12)
    mm <- Datasett$Month[selected]
    yy <- Datasett$Year[selected]
    dat <- eval(parse(text=paste("Datasett$",param,"[selected]",sep="")))
    #print("make station.series object")
    #x11(); plot(yy+(mm-0.5)/12,dat,type="l"); print(rbind(Datasett$Dato,yy,mm,dat)); print(table(mm)); print(table(yy)); print(class(mm)); print(class(yy))
    kdvh <- station.obj(x=dat,yy=yy,mm=mm,obs.name=param,unit=unit,
                        station= StNr,lat=lat,lon=lon,alt=alt,
                        location=location,country="Norway",ele=ele)
    #print("kdvh made")
    #print(kdvh$val[1:3,])
  } else {   
    Datasett <- read.table(Filnavn,header=TRUE)
    day <- Datasett$Day
    mon <- Datasett$Month
    year <-Datasett$Year
    
    ny <- max(year,na.rm=TRUE)-min(year,na.rm=TRUE)+1
    yy <- sort(rep(seq(min(year,na.rm=TRUE),max(year,na.rm=TRUE),by=1),12))
    dat<-rep(NA,length(yy))
    mm <- rep(1:12,ny)
 #print(c(length(yy),length(mm),ny,e.rows.mon))
  
   if (month) {
     for (i in 1:e.rows.mon){
       ii <- is.element(mon,mm[i]) & is.element(year,yy[i])
       if (sum(ii)>min.data.points) {
          dat[i] <- eval(parse(text=paste("round(",method,
                         "(as.numeric(Datasett$",param[1],"[ii]),na.rm=TRUE),1)")))
       }
#    print(c(i,mm[i],yy[i],sum(ii),dat[i]))
     }
     if (now("c(day,mon,year)")[1] < 28) dat[yy==now("c(day,mon,year)")[3] & mm==now("c(day,mon,year)")[2]] <- NA
       kdvh <- station.obj(x=dat,yy=yy,mm=mm,obs.name=param,unit=unit,
                           station= StNr,lat=lat,lon=lon,alt=alt,
                           location=location,country="Norway",ele=ele)    
     } else {
      data.1 <- eval(parse(text=paste("Datasett$",param[1],sep="")))
      if (length(param>=2)) {
        if (!is.na(param[2])) data.2 <- eval(parse(text=paste("Datasett$",param[2],sep=""))) else data.2 <- rep(NA,length(data.1))
      } else data.2 <- rep(NA,length(data.1))
      kdvh <- station.obj.dm(t2m=data.1,precip=data.2,
                             yy=year,mm=mon,dd=day,obs.name=param,unit=unit,
                             station= StNr,lat=lat,lon=lon,alt=alt,
                             location=location,country="Norway",ele=ele)
     }  
   }

   kdvh$utm.e <- dnmi.meta$"Utm_e"[dnmi.meta$Stnr==StNr][1]
   kdvh$utm.n <- dnmi.meta$"Utm_n"[dnmi.meta$Stnr==StNr][1]
   kdvh$utm.sone <- dnmi.meta$"Utm_sone"[dnmi.meta$Stnr==StNr][1]
   kdvh$kommune <- dnmi.meta$Kommune[dnmi.meta$Stnr==StNr][1]
   kdvh$fylke <- dnmi.meta$Fylke[dnmi.meta$Stnr==StNr][1]
   kdvh$type <- dnmi.meta$Type[dnmi.meta$Stnr==StNr][1]
   kdvh$type.beskr <- dnmi.meta$"Type_beskrivelse"[dnmi.meta$Stnr==StNr][1]
   kdvh$fnr <- dnmi.meta$Fnr[dnmi.meta$Stnr==StNr][1]
   kdvh$knr <- dnmi.meta$Knr[dnmi.meta$Stnr==StNr][1]
   kdvh$URL <- Filnavn
   kdvh$method <- method
   invisible(kdvh)
}

avail.kvdh <- function(yr0=1950,mo0=1,print=TRUE) {
  data(dvh.station.list,envir=environment())
#  dnmi.meta <- read.table("data/dvh.station.list",header=TRUE,as.is=TRUE)
  iavail <- (dnmi.meta$yy0 <= yr0) 
  result <- list(Navn=dnmi.meta$Navn[iavail],Stnr=dnmi.meta$Stnr[iavail],
                 yy0=dnmi.meta$yy0[iavail],Type=dnmi.meta$Type[iavail])
  if (print) print(cbind(result$Navn,result$Stnr,result$yy0,result$Type))
  invisible(result)
}




untar <- function(path="/disk1/sesong/data/tar_nc",untar.path="/home/rasmusb/ECFCSTS") {
  file.list <- list.files(path=path,pattern=".tar")
  #print(file.list)
  wd0 <- getwd()
  setwd(untar.path)
  for (file in file.list) {
    command <- paste("tar xvf ",path,"/",file,sep="")
    #print(command)
    system(command)
  }

  file.list <- list.files(path="/disk1/sesong/data/nc_prdato",pattern=".nc")
  print(file.list)
  for (file in file.list) system(paste("ln -s /disk1/sesong/data/nc_prdato/",file,sep=""))

  setwd(wd0)
}

product <- function(x) {
  product <- exp(sum(log(x)))
  product
}


# path="/opdata/sesong/"
# path="/disk1/sesong/data/smma/"
# path="/home/rasmusb/ECFCSTS/smma/"


reliab.diag <- function(obs,ts,y0,fp=NULL,CDF=FALSE,res=300,type="png256",loc,a,
                        yy=NULL,mm=NULL,dd=NULL,FCST=NULL,test=FALSE,lag=1) {

  if (dev.cur() ==1) newFig()
  
  if (is.null(FCST)) {
    FCST <- list(vname=obs$obs.name)
    class(FCST) <- "monthly forecasts"
  }
  if (is.null(yy)) yy <- plotStation(obs,what="t")$yy; dev.off()
  if (is.null(mm)) mm <- plotStation(obs,what="t")$mm; dev.off()
  if (class(obs)[1]=="station") obs <- c(t(obs$val))
  print("reliab.diag:")
  if (class(FCST)=="monthly.forecasts") {
    print("monthly forecasts")
    i1 <- is.finite(obs)
    i2 <- is.finite(obs)
  } else {
    print("seasonal forecasts")
    print(c(length(yy),length(mm),NA,length(obs),NA,length(a$yy),length(a$mm)))
    print(paste("Observations: time ", min(yy*100+mm),"-",max(yy*100+mm),"nt=",length(yy)))
    print(paste("fcst: time ", min(a$yy*100+a$mm),"-",max(a$yy*100+a$mm),"nt=",length(a$yy)))
    i1 <- is.element(yy*100+mm,a$yy*100+a$mm) 
    i2 <- is.element(a$yy*100+a$mm,yy*100+mm)
  }
  if (sum(i1)==0 | sum(i2)==0) {
    #print(rbind(yy*100+mm,a$yy*100+a$mm))
    stop('reliab.diag: No matching dates!')
  } else print(paste("No. matching dates is",sum(i1),"=",sum(i2)))
  if (sum(i1) != sum(i2)) stop('reliab.diag: Messed up chronology!')
  #print(dim(ts)); print(length(obs))
  ts <- ts[i2,]; obs <- obs[i1]
  good <- is.finite(obs)
  ts <- ts[good,]; obs <- obs[good]
  if (length(obs) != dim(ts)[1]) stop(paste("reliab.diag: vectors must have same lengths",
                                             length(obs), dim(ts)[1]))
  nf <- length(obs)
  p <- rep(NA,nf)
  print(c(floor(min(obs,na.rm=TRUE)),ceiling(max(obs,na.rm=TRUE))))
  q <- seq(floor(min(c(ts,obs),na.rm=TRUE)),ceiling(max(c(ts,obs),na.rm=TRUE)),length=300)
  print(paste("Forecaste probabilities: nf=",nf))

  for (i in 1:nf) {
    cdf <- pnorm(q,mean=mean(ts[i,]),sd=sd(ts[i,]))
    edf <- empirical.ranking(ts[i,])
    if (CDF) p[i] <- max(cdf[q < y0],na.rm=TRUE) else
             p[i] <- max(edf$P[edf$x < y0],na.rm=TRUE)
  }

  if (sum(is.finite(p))==0) {print(summary(ts)); stop("reliab.diag: invalid probabilities")}
  
  if (is.null(fp)) fp <- seq(0,1,length=5)

  nb <- length(fp)
  fp[nb] <- 1.001
  orf <- rep(0,nb)
  print(paste("obs rel. probabilities: nf=",nb))
  print(paste("i  fp[i]   sum(iii)   sum(obs[iii]<=",y0,")"))
  for (i in 2:nb) {
    iii <- (p >= fp[i-1]) & (p < fp[i]) 
    if (sum(iii)>0) orf[i] <- sum(obs[iii] < y0)/sum(iii) else orf[i] <- 0
    print(paste("obs.rel.freq.: i=",i,"fp=",fp[i],"sum(iii)=",sum(iii),
                " range:",round(min(obs[iii],na.rm=TRUE),2),",",round(median(obs[iii],na.rm=TRUE),2),
                ",",round(max(obs[iii],na.rm=TRUE),2),
                " sum(obs[iii]) <y0:", sum(obs[iii] < y0)," orf=",round(orf[i],2)))
  }

  if (test) {
     print("---TEST---")
     print("p, obs, y0, fp, orf:")
     print(summary(p)); print(summary(obs)); print(y0); print(fp); print(orf)
     cols <- c("darkred","red","pink","green","darkgreen","lightblue","blue",
               "darkblue","wheat","yellow","magenta","cyan")
     bitmap(paste("reliab.diag-test-",loc,FCST$vname,"-",y0,"-",res,".",
                  substr(type,1,3),sep=""),type=type,res=res)
     par(mfcol=c(2,1))
     plot(p,ylim=c(-0.1,1),main="Probabilities")
     for (i in 1:nb) lines(c(1,length(p)), rep(fp[i],2),lty=2,col="steelblue")
     for (i in seq(nb,1,by=-1)) {
       iii <- (p <= fp[i])
       points((1:length(p))[iii],p[iii],pch=20,cex=0.5,col=cols[i])
       text(i*length(p)/12,-0.075,fp[i],col=cols[i],cex=1.2)
     }

     plot(c(1,length(obs)),range(c(ts,obs),na.rm=TRUE),type="n",main="Forecast & obs.")
     grid()
     for (i in 1:length(ts[1,])) points(ts[,i],col="grey",cex=0.7)
     points(obs,col="red")
     for (i in seq(nb,1,by=-1)) {
       iii <- (p <= fp[i])
       points((1:length(p))[iii],obs[iii],pch=20,cex=0.5,col=cols[i])
     }  

     
     dev.off()
     
     print("TEST - p summary:")
     print(summary(p))
  }

  h <- hist(p)
  #print("Brier score:")
  bs <- round(BS(obs,ts,y0,CDF=CDF),2)
  #print("Graphics:")
  #print(rbind(fp,orf))

  print(paste("Filename: reliab.diag",loc,FCST$vname,"-",y0,"_",lag,"-",res,".",
               substr(type,1,3),sep=""))
#  bitmap(paste("reliab.diag",loc,FCST$vname,"-",y0,"_",lag,"-",res,".",
#               substr(type,1,3),sep=""),type=type,res=res)
  x11()
  plot(c(0,1),c(0,1),type="l",lwd=3,col="grey70",main=paste("Reliability diagram for",FCST$vname),
       xlab="Forecast probability, yi",ylab="Observed relative frequency",
       sub=paste("Pr(x <",y0,"), Brier Score= ",bs))
  grid()
  lines(fp,orf,lty=2)
  points(fp,orf,pch=19)
  if (class(FCST)=="monthly.forecasts") mtext(side=4,paste("Monthly forecasts,",length(ts[1,]),
                                                "ensemble members & ",length(ts[,1]),"forecasts"))

  #print("Insert:")
  n <- length(h$breaks)*3
  hx <- rep(NA,n+1); hy <- hx
  ii <- 0
  for(i in seq(1,n,by=3)) {
    ii <- ii+1
    hy[i] <- 0
    hy[i+1:2] <-   0.25 *rep(h$density[ii],2)/max(h$density)
    hx[i:(i+1)] <- 0.30 * h$breaks[ii]/max(h$breaks)
    hx[i+2] <-     0.30 * h$breaks[ii+1]/max(h$breaks)
  }
  hy[n+1] <- 0; hx[n+1] <- 0.30

  if (mean(orf[1:5],na.rm=TRUE) < 1.1* mean(fp[1:5],na.rm=TRUE)) {
    x1 <- 0; x2 <- 0.3; y1 <- 0.7;y2 <- 1.0
  } else {
    x1 <- 0.7; x2 <- 1.0; y1 <- 0;y2 <- 0.3
  }
  hy <- hy + y1; hx <- hx + x1
      
  polygon(c(x1,x2,x2,x1,x1),c(y1,y1,y2,y2,y1),
          border="black",col="grey80",lwd=2)
  lines(hx,hy)
  dev.copy2eps(file=paste("reliab.diag",loc,FCST$vname,"-",y0,"_",lag,"-",res,".eps",sep=""))
#  dev.off()

  reliability <- list(fp=fp,orf=orf,bs=bs,h=h)
  class(reliability) <- "reliability"
  invisible(reliability)
}


# Ny routine: season.mean REB 13.02.2006
seasonal.mean <- function(x) {
  seasonal.mean <- c(ma.filt(x,3)[-1],NA) # REB 13.02.2006
  seasonal.mean
}


getecfcst <- function(vname="T2M",path="/opdata/sesong/",local=TRUE,
                      verbose=FALSE,precip.local="~/data/sesong/") {

  if (path=="/opdata/sesong/") name.id <- "smma"
  if (substr(path,1,26)=="/klimadata/rasmusb/monthly") name.id <- "MonthlyFcst"
   if ( (local) & file.exists(precip.local) ) {
     #print("Use local acrhive too!")
     op.data.file <- list.files(path=path,pattern=".nc")
     op.data.file <- op.data.file[grep(name.id,op.data.file)]
     #print(op.data.file)
     local.data.file <- list.files(path=precip.local,pattern=vname)
     local.data.file <- local.data.file[grep(name.id,local.data.file)]
     cp.files <- op.data.file[!is.element(op.data.file,local.data.file)]
     for (cp.file in cp.files) {
       command <- paste("system('cp -f ",path,cp.file," ",precip.local,".')",sep="")
       print(command)
       eval(parse(text=command))
     }
     path <- precip.local
   }
   print(paste("Data path=",path))
#   vname2 <- switch(vname,"T2M"="2T","MSL"="MSL","TP"="TP",
#                    "2t"="2T","p2t"="2T","tp"="TP","msl"="MSL",
#                    "t2m"="2T","slp"="MSL","temp"="2T","pre"="TP","prec"="TP",
#                    "precip"="TP","167.171"="2T")
   vname2 <- switch(vname,"T2M"="2T","MSL"="MSL","TP"="TP",
                    "2t"="2T","p2t"="2T","tp"="TP","msl"="MSL",
                    "t2m"="2T","slp"="MSL","temp"="2T","pre"="TP","prec"="TP",
                    "precip"="TP","167.171"="2T","2ta"="2T","v2ta"="2T")
   file.list <- list.files(path=path,pattern=".nc")
   file.list <- file.list[c(grep(paste(".",vname2,".",sep=""),file.list),
                            grep(paste(".",vname,".",sep=""),file.list))]
   file.list <- file.list[grep(name.id,file.list)]; file.list <- sort(file.list)
   print(vname2)
   print(file.list)
   nf <- length(file.list)
   for (ifcst in 1:nf) {
     dots <- instring(".",file.list[ifcst])
     fcst <- substr(file.list[ifcst],dots[2],dots[3]-1)
     print(file.list[ifcst])
     ncid<-open.ncdf(paste(path,file.list[ifcst],sep=""))
     nv <- ncid$nvars
     cdfvars <- rep("-",nv)
     for (i in 1:nv) cdfvars[i] <- ncid$var[[i]]$name
     miss <-  ncid$var[[1]]$missval
     nd <- ncid$var[[1]]$ndims
     cdfdims <- rep("-",nd)
     for (i in 1:nd) cdfdims[i] <- ncid$var[[1]]$dim[[i]]$name
     #print(cdfdims)

     ilon <- grep("lon",lower.case(cdfdims))
     ilat <- grep("lat",lower.case(cdfdims))
     itim <- grep("tim",lower.case(cdfdims))
     inum <- grep("num",lower.case(cdfdims))
     ifcp <- grep("fcperiod",lower.case(cdfdims))
     lon <- get.var.ncdf(ncid,cdfdims[ilon])
     lat <- get.var.ncdf(ncid,cdfdims[ilat])

     if (ifcst==1) { lon0 <- lon; lat0 <- lat }
     num <- get.var.ncdf(ncid,cdfdims[inum])
     tim <- get.var.ncdf(ncid,cdfdims[itim])
     if (length(ifcp)>0) fcperiod <- get.var.ncdf(ncid,cdfdims[ifcp])
     
     attr(lon,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[ilon],"$units",sep="")))
     attr(lat,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[ilat],"$units",sep="")))
     attr(tim,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[itim],"$units",sep="")))
     attr(num,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[inum],"$units",sep="")))
     t.unit <- attr(tim,"unit")

     torg <- substr(t.unit,regexpr("since",t.unit)+6,nchar(t.unit))
     t.unit <- substr(t.unit,1,regexpr("since",t.unit)-2)
  
     if (!is.null(torg)) {
       print(paste("torg=",torg))
       yy0 <- datestr2num(torg)[1]
       mm0 <- datestr2num(torg)[2]
       dd0 <- datestr2num(torg)[3]
     } 

     if (verbose) print(paste("Time unit:",lower.case(t.unit)))
     mmddyy <- caldat(tim/24 + julday(mm0,dd0,yy0))
     mm <- mmddyy$month
     yy <- mmddyy$year
     dd <- mmddyy$day
     tim <- tim/24
     t.unit <- "day"
     obj.type <- "daily.field.object"

     daysayear<- 365.25
     nx <- length(lon); ny <- length(lat); nt <- length(tim)
     v1 <- cdfvars[1]
     data <- get.var.ncdf(ncid,v1)

     arv <- att.get.ncdf(ncid, cdfvars[1], 'scale_factor')
     if( arv$hasatt ) scal <- arv$value else scal <- 1
     arv <- att.get.ncdf(ncid, cdfvars[1], 'add_offset')
     if( arv$hasatt ) offs <- arv$value else offs <- 0
     if (as.numeric(R.Version()$major) < 2) { 
       print(paste("The host is",Sys.info()[4],"running R version",R.Version()$major,"-",
                    R.Version()$minor,": scaling by x",scal,"and adding",offs))
       data <- data * scal
       data <- data + offs
     }
     
     close.ncdf(ncid)

     if ( (length(lon) != length(lon0)) | (length(lat) != length(lat0)) |
          (lon[1] != lon0[1]) | (lat[1] != lat0[1]) |
          (lon[length(lon)] != lon0[length(lon0)]) | (lat[length(lat)] != lat0[length(lat0)]) ) {
       print("-- Detecting different grids! ---")
       print("Apply interpolations")
       print(dim(data))
       dims <- dim(data)
       nt <- length(tim); nx <- length(lon0); ny <- length(lat0); nm <- length(num)
       new.dims <- dims; new.dims[1] <- nx; new.dims[2] <- ny
       #print(new.dims);
       vector.length <-  product(new.dims)
       lonxy <- rep(lon,length(lat)); latxy <- sort(rep(lat,length(lon)))
       #print("Check!")
       dat <- rep(NA,vector.length); print(vector.length); print(length(dat))
       if (length(dat) != vector.length) {
         print("something strange here - this should not happen!")
         dat <- c(dat,NA); print(vector.length); print(length(dat))
       }
       dim(dat) <- new.dims
       print(paste("nm=",nm,"length(fcperiod)=",length(fcperiod),"nt=",nt))
       for (im in 1:nm) {
         for (ifc in 1:length(fcperiod)) {
           for (it in 1:nt) {
             map <- data[,,im,ifc,it] 
             #print(c(length(lon),length(lat),NA,dim(map),NA,product(dim(map)),length(lonxy),length(latxy)))
             dat[,,im,ifc,it] <- interp(lonxy,latxy,map,lon0,lat0)$z 
           }
         }
       }
       lon <- lon0; lat <- lat0; data <- dat; rm(dat)
     }
     
     #if ((v1 == "2t") | (v1 == "p2t")) data <- data - 273

     lon[lon > 180] <- lon[lon > 180]-360

     if (verbose) print("Sort longs and lats")
     x.srt <- order(lon)
     y.srt <- order(lat)
     lon <- lon[x.srt]
     lat <- lat[y.srt]
     if (verbose) print(dim(data))
     print(v1)
     v1 <- switch(lower.case(v1),"2t"="t2m","p2t"="t2m","tp"="tp","msl"="slp",
                  "t2m"="t2m","slp"="slp","temp"="t2m","pre"="tp","prec"="tp",
                  "precip"="tp","167.171"="t2m","228.173"="tp","2ta"="t2m",
                  "v2ta"="t2m","tpara"="tp","t2a"="t2m")
     print(paste("Variable name:",v1))

     if (length(cdfdims)==4) {
       command <- paste("ecfcst <- list(",v1,"= data,lon=lon,lat=lat,num=num,tim=tim)",sep="")
       print(command)
       data <- data[x.srt,y.srt,,]
       eval(parse(text=command))
       ecfcst$yy <- yy;  ecfcst$mm <- mm; ecfcst$dd <- dd
       #print(summary(fcst))

       if (ifcst == 1) {
         if (verbose) print(paste("FCST <- list(fcst",fcst," = ecfcst)",sep=""))
         eval(parse(text=paste("FCST <- list(fcst",fcst," = ecfcst, N=nf, nx=nx, ny=ny ,vname=v1)",sep=""))) 
         if (verbose) print(summary(FCST))
         FCST$lon <- lon; FCST$lat <- lat; FCST$num <- num
       }  else {
         eval(parse(text=paste("FCST$fcst",fcst," <- ecfcst",sep="")))
         if (verbose) print(summary(FCST))
       }
   } else {
       print("5D data matrix!"); print(dim(data))
       command <- paste("ecfcst <- list(",v1,"= data,lon=lon,lat=lat,num=num,tim=fcperiod)",sep="")
       print(command)
       data5D <- data
       for (it in 1:nt) {
         data <- data5D[x.srt,y.srt,,,it]
         eval(parse(text=command))
         ecfcst$yy <- yy;  ecfcst$mm <- mm; ecfcst$dd <- dd
         #print(summary(fcst))

       if ( (ifcst == 1) & (it==1) ) {
         if (verbose) print(paste("FCST <- list(fcst",fcst," = ecfcst)",sep=""))
         eval(parse(text=paste("FCST <- list(fcst",fcst," = ecfcst, N=nf, nx=nx, ny=ny ,vname=v1)",sep=""))) 
         if (verbose) print(summary(FCST))
         FCST$lon <- lon; FCST$lat <- lat; FCST$num <- num
       }  else {
         eval(parse(text=paste("FCST$fcst",fcst," <- ecfcst",sep="")))
         if (verbose) print(summary(FCST))
       }
     }
   }
  }  # End loop

  if (as.numeric(R.Version()$major) < 2) { 
       print(paste("The host is",Sys.info()[4],"running R version",R.Version()$major,"-",
                    R.Version()$minor,": scaling by x",scal,"and adding",offs))
  } else {
       print(paste("The host is",Sys.info()[4],"running R version",R.Version()$major,"-",
                    R.Version()$minor,": No further scaling."))
  }

  FCST$dat.att <- cdfcont(paste(path,file.list[1],sep=""))
  class(FCST) <- "Ensemble seasonal forecast"
  save(file="getecfcst.rda",FCST)
  invisible(FCST)

}  # end function




ensemble.timeseries <- function(FCST,fcst.list,n.ens,obs,lag,verbose=FALSE) {

  if (verbose) print("ensemble.timeseries - start")
  cmon <- c("j","f","m","a","m","j","j","a","s","o","n","d")
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny
  dat <- rep(NA,nf*nx*ny); dim(dat) <- c(nf,ny,nx)
  tim <- 1:nf; yy <- rep(NA,nf); mm <- yy; dd <- yy
  id.x <- matrix(rep(FCST$vname,ny*nx),ny,nx); id.t <- rep(FCST$vname,nf)
  ts <- rep(NA,nf*n.ens); dim(ts) <- c(nf,n.ens)
  mon3 <- rep(" ",nf)
  
  for (inum in 1:n.ens) {
    for (ifcst in 1:nf) {
      cline <- paste("FCST$",fcst.list[ifcst],"$",FCST$vname,sep="")
      yy[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$yy[",lag+1,"]",sep="")))   # REB 15.02.2006: '+1' since first field is month 0.
      mm[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$mm[",lag+1,"]",sep="")))
      dd[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$dd[",lag+1,"]",sep="")))
      fcst <- eval(parse(text=cline))
      months <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$mm[",lag+1,":",lag+3,"]",sep="")))
      months[months ==13] <- 1
      mon3[ifcst] <- paste(cmon[months[1]],cmon[months[2]],cmon[months[3]],sep="")
      #print(c(inum,lag,lag+2,NA,dim(t(fcst[,,inum,lag])),NA,dim(dat[ifcst,,])))
      dat[ifcst,,] <- 0.3333* t(fcst[,,inum,lag+1] + fcst[,,inum,lag+2] + fcst[,,inum,lag+3])
    }  
 
    if (lower.case(FCST$vname)=="p2t") dat <- dat - 273
    if (lower.case(FCST$vname)=="tp")  dat <- dat * 3600*24*30 * 1000
    if (verbose) print(summary(c(dat)))
                                        # Needs changing to 28, 29, 30 or 31 days... XXX

    print(FCST$obj.type)
    FCST$obj.type <- "monthly.field.object"
    attr(tim,"units") <- "month"
    print(attributes(dat))
    field  <- list(dat=dat,lon=FCST$lon,lat=FCST$lat,tim=tim,
                   v.name=FCST$vname,id.x=id.x,id.t=id.t,
                   yy=yy,mm=mm,dd=dd,n.fld=1,
                   id.lon=rep(FCST$vname,FCST$nx),id.lat=rep(FCST$vname,FCST$ny),
                   attributes=NULL)
    class(field) <- c("field",FCST$obj.type)
    print(class(field))
    
    if (verbose) print(summary(c(dat)))
    st <- plotField(field,lon=obs$lon,lat=obs$lat,what="abs",col="grey40",add=TRUE)
#    class(st) <- c("station","monthly.station.record")
    print(summary(st)); print(class(st))
    a <- plotStation(st,l.anom=FALSE,what="t",add=TRUE,trend=FALSE,std.lev=FALSE,col="pink")
    print(c(inum,length(ts[,inum]),length(a$value),dim(ts)))
    ts[,inum] <- a$value
  }

  x11()

  clim.lag <- lagStation(obs,lag=lag)
  clim.lag <- rep(colMeans(clim.lag$val[is.element(clim.lag$yy,1961:1990),],na.rm=TRUE),3)
  clim.lag<-ma.filt(clim.lag,3)
  clim.lag <- clim.lag[13:24]
  
  b <-  plotStation(obs,l.anom=FALSE,what="t",trend=FALSE,std.lev=FALSE,col="blue")
  dev.off()

  ylab <- paste(obs$obs.name," (",obs$unit,")")
  if (verbose) print("Time interval for ensemble")  
  if (verbose) print(range(a$yy*100+a$mm))
  if (verbose) print("ensemble.timeseries - finished")
  i1 <- is.element(a$yymm,b$yymm); i2 <- is.element(b$yymm,a$yymm)
  y2 <- filter(b$value,rep(1,3)/3,sides=1)[i2]
  r <- rep(NA,n.ens); rmse=r; R2 <- r; bias=r; scale=r
  for (inum in 1:n.ens) {
    y1 <- ts[i1,inum]; ok <- is.finite(y1) & is.finite(y2) 
    r[inum] <- round(cor(y1[ok],y2[ok]),2)
    rmse[inum] <- round(sqrt(sum( (y1[ok]-y2[ok])^2 ))/nf,2)
    R2[inum] <- r[inum]^2
    bias[inum] <- mean(y2[ok]) - mean(y1[ok])
    scale[inum] <- sd(y2[ok])
  }
  
  Clim <- rep(NA,length(ts[,1]))
  for (im in 1:12) {
    ii <- is.element(a$mm,im)
    Clim[ii] <- clim.lag[im]
  }

  ensemble.timeseries <- list(ts=ts,a=a,months=months,mon3=mon3,location=obs$location,clim=Clim,
                              year=a$yy,month=a$mm,bias=bias,scale=scale,
                              ylab=ylab,r=r,rmse=rmse,R2=R2,model=rep("IFS",n.ens),run=1:n.ens,nmnts=3)
  class(ensemble.timeseries) <- "EC.ensemble"
  invisible(ensemble.timeseries)
}



adjust.ts <- function(obs,test=FALSE,fill=TRUE) {
  klima <- obs
  klima$val <- obs$val - anomaly.station(obs)$val
  bitmap(file="adjust.ts-test1.png",type="png256")
  Klima <- plotStation(klima,l.anom=FALSE,what="t",add=FALSE,trend=FALSE,std.lev=FALSE,col="blue")
  yy <- Klima$yy; mm <- Klima$mm; klima <- Klima$value
  klima.87.01 <- obs 
  klima.87.01$val <-  obs$val - anomaly.station(obs,period=c(1987,2001))$val
  klima.87.01 <- plotStation(klima.87.01,l.anom=FALSE,what="t",add=TRUE,trend=FALSE,std.lev=FALSE,col="blue")$value
  dev.off()
  
  if (test) {
    bitmap(paste("getecfcst.adjust.ts-test_",obs$station,".png",sep=""),type="png256")
    plot(klima[1:12],type="l",lwd=3)
    lines(klima.87.01[1:12],col="red",lty=2,lwd=2)
    for (ii in 1:49) points(obs$val[ii,],col="grey",pch=20)
    for (ii in 1:49) points(obs$val[ii,]-anomaly.station(obs)$val[ii,])
    lines(obs$val[1,]-anomaly.station(obs,period=c(1987,2001))$val[1,],col="blue")
    lines(mm[1:12],Klima$value[1:12],lty=2,col="darkblue")
    while (dev.cur() > 1) dev.off()
  }

  klima <- seasonal.mean(klima)
  klima.87.01 <- seasonal.mean(klima.87.01)

  #print(c(sum(is.finite(klima)),sum(is.finite(klima.87.01))))
  nl <- length(klima.87.01)
#print(length(klima.87.01[(nl-11):nl])); print(length(klima.87.01[(nl-23):(nl-12)]))

# Fill in climatological values where there are NAs.
#  print(rbind(klima[(nl-11):nl],klima[(nl-23):(nl-12)]))

  if (fill) {
#    klima.87.01[(nl-11):nl] <- klima.87.01[(nl-23):(nl-12)]        # REB 17.03.2006
#    klima[(nl-11):nl] <- klima[(nl-23):(nl-12)]                    # REB 17.03.2006

# New lines: REB 17.03.2006
    fill.in <- seq(1,length(klima),by=1)[!is.finite(klima)]
    if (length(fill.in)>0) {
      for (i in 1:length(fill.in)) klima[fill.in[i]] <- mean(klima[is.element(mm,mm[fill.in[i]])],na.rm=TRUE)
    }
    fill.in <- seq(1,length(klima.87.01),by=1)[!is.finite(klima.87.01)]
    if (length(fill.in)>0) {
      for (i in 1:length(fill.in)) klima.87.01[fill.in[i]] <- mean(klima.87.01[is.element(mm,mm[fill.in[i]])],na.rm=TRUE)
    }
  }
# End of new lines: REB 17.03.2006
  
  results <- list(klima=klima,klima.87.01=klima.87.01,mm=mm,yy=yy)
  invisible(results)
}



make.epsfig <- function(loc,obs,obs.ts,yy,mm,FCST,res=300,type="png256",a,ts,ii2,
                        y,ens.mean,mon3,
                        label=TRUE,fcst.list,months,verbose=FALSE,n.ens) {

  if (verbose) print("make.epsfig - start")
  cmon <- c("j","f","m","a","m","j","j","a","s","o","n","d")
  cmonth <- c("januar","februar","mars","april","mai","juni",
            "juli","august","september","oktober","november","desember")
  vname=switch(FCST$vname,"p2t"="2m temperatur (C)","t2m"="2m temperatur (C)",
               "tp"="Nedbor","msl"="Navniva lufttrykk")
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny
  year <- a$yy[nf]
  if (label) main <- paste(obs$location," dato:",cmonth[months[1]],year,
                           "til folgende",cmonth[months[3]]) else
             main <- " "


  bitmap(paste("kvalitet",loc,FCST$vname,"-",res,".",substr(type,1,3),sep=""),type=type,res=res) 
                                
  par(cex.main=0.8)
  plot(range(c(a$yy,a$yy+1)),range(c(y[is.finite(y)],c(ts[is.finite(ts)])),na.rm=TRUE)*1.05,
       type="n",xlab="Dato",ylab=paste(vname,"anomali"),main=main,
       sub=paste("Interpolert til ",obs$lat,"N/",obs$lon,"E. (station #",loc,")",sep=""))
  
  grid()
  col1 <- rep("grey75",nf); col2 <- rep("grey80",nf)
  col1[is.element(mon3,mon3[nf])] <- "grey35"
  col2[is.element(mon3,mon3[nf])] <- "grey70"
  
  for (i in 1:nf) {
    lines(rep(a$yy[i]+(a$mm[i]-0.5)/12,2),range(ts[i,])*1.005,col=col1[i],lwd=5)
    lines(rep(a$yy[i]+(a$mm[i]-0.5)/12,2),range(ts[i,]),col=col2[i],lwd=3)
    points(rep(a$yy[i]+(a$mm[i]-0.5)/12,n.ens),ts[i,],pch=20,cex=0.5,col="grey30")
    text(a$yy[i]+(a$mm[i]-0.5)/12,min(c(y[i],c(ts[i,])),na.rm=TRUE)-0.3,mon3[i],cex=0.7,srt=90)
    points(a$yy[i]+(a$mm[i]-0.5)/12,mean(ts[i,]),pch=20,cex=2)
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=1.25,col="grey10")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=1.00,col="grey20")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=0.85,col="grey30")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=0.70,col="grey40")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=0.60,col="grey50")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=0.50,col="grey60")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=0.40,col="grey70")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=0.30,col="grey80")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=0.20,col="grey90")
    points(a$yy[i]+(a$mm[i]-0.5)/12+0.01,mean(ts[i,])+0.01,pch=20,cex=0.10,col="white")
  }
  now <- length(y); now2 <- length(obs.ts$value)
  print(paste("Anomaly on ",yy[now]*100+mm[now],"is",y[now]," ( total value is",round(obs.ts$value[now2],2),
                            100*obs.ts$yy[now2]+obs.ts$mm[now2],")"))
                                                                                
  points(yy+(mm-0.5)/12, y,pch=20,cex=2.1,col="red")
  points(yy+(mm-0.5)/12, y,pch=21,cex=1.9,col="darkred")

  I1 <- is.element(100*yy+mm,100*a$yy+a$mm)
  I2 <- is.element(100*a$yy+a$mm,100*yy+mm)
  eval <- list(fcst=ts[I2,],obs=y[I1],year=yy[I1],month=mm[I1])
  save(file=paste("data/make.epsfig",loc,FCST$vname,"-",res,".",substr(type,1,3),
         ".rda",sep=""),eval)

# Most recent forecast:

  i <- nf
  lines(rep(a$yy[i]+(a$mm[i]-0.5)/12,2),range(ts[i,])*1.005,col="red",lwd=5)
  lines(rep(a$yy[i]+(a$mm[i]-0.5)/12,2),range(ts[i,]),col="pink",lwd=3)
  points(rep(a$yy[i]+(a$mm[i]-0.5)/12,n.ens),ts[i,],pch=20,cex=0.5,col="red")
  points(a$yy[i]+(a$mm[i]-0.5)/12,mean(ts[i,]),pch=20,col="steelblue",cex=1.5)
  points(a$yy[i]+(a$mm[i]-0.5)/12,mean(ts[i,]),pch=4,col="darkred",cex=1.25)
  points(a$yy[i]+(a$mm[i]-0.5)/12,mean(ts[i,]),pch=21,col="darkred",cex=1.5)

  text(a$yy[i]+(a$mm[i]-0.5)/12,min(c(y[i],c(ts[i,])),na.rm=TRUE)-0.3,mon3[i],
       cex=0.7,srt=90,font=2)
# Legends...

  ens.mean <- rowMeans(ts[ii2,]); good <- is.finite(y) & is.finite(ens.mean)
#  if (verbose) print(c(dim(ts),NA,length(y))); print(c(length(rowMeans(ts[ii2,])),length(colMeans(ts[ii2,])))); print(rbind(y,ens.mean,yy,mm))

#  cor.score <- round(cor(y[good],ens.mean[good]),2)   # Misleading because of recent high level...
  cor.score <- round(sum(y[good]*ens.mean[good])/sqrt(sum(y[good]^2)*sum(ens.mean[good]^2)),2)
  text(mean(a$yy),max(c(y,c(ts)),na.rm=TRUE),paste("korrelasjon m. ens.middel=",cor.score))
  if (label) mtext(side=4,paste("(ecfcst2eps - getecfcst.R) - ",date()))
 dev.off()
 if (verbose) print("make.epsfig - finished")
}



#--------------------------------------------------------------------------

#ecfcst2eps <- function(FCST,lag=1,loc=18700,n.ens=40,label=TRUE,verbose=FALSE,ensemble=NULL,
#                       test=FALSE,mean.from.daily=FALSE) {
ecfcst2eps <- function(FCST,lag=1,loc=18700,n.ens=51,label=TRUE,verbose=FALSE,ensemble=NULL,
                       test=FALSE,mean.from.daily=FALSE) {

  if (verbose) print("-----------------------------------------------------------------------------")
  if (verbose) print(paste("ecfcst2eps: lag=",lag,"  loc=",loc,"  n.ens=",n.ens,
                           "  label=",label))
  bitmap(file="scratch.png",type="png256")
  cmon <- c("januar","februar","mars","april","mai","juni",
            "juli","august","september","oktober","november","desember")
  if (is.null(FCST)) load("getecfcst.rda")
  fcst.list <- rownames(summary(FCST))
  fcst.list <- fcst.list[grep("fcst.",fcst.list)]

  ele <- switch(FCST$vname,"t2m"="TAM","tp"="RR","slp"="POM")
  obs <- KDVH4DS(loc,param=ele,mean.from.daily=mean.from.daily)
  obs.ts <- plotStation(obs,l.anom=FALSE,what="t",trend=FALSE,std.lev=FALSE)
  plotStation(obs,l.anom=TRUE,what="t",trend=FALSE,std.lev=FALSE)
  # if (verbose) {print("Monthly values:"); print(rbind(obs.ts$yy,obs.ts$mm,round(obs.ts$value,1)))}
  obs.ts$value <- seasonal.mean(obs.ts$value)
  
  if (is.null(ensemble)) ensemble <- ensemble.timeseries(FCST,fcst.list,n.ens=n.ens,obs=obs,lag=lag)
  a <- ensemble$a; ts <- ensemble$ts; months <- ensemble$months; mon3 <- ensemble$mon3
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny

  if (verbose) print(summary(a))

  # 1987-2001 Climatology - extend to the whole year (into the future)...
  if (verbose) print("Climatologies")
  ii <- is.element(obs.ts$yy*100 + obs.ts$mm, a$yy*100 + a$mm)
  ii2 <- is.element(a$yy*100 + a$mm,obs.ts$yy*100 + obs.ts$mm)
  klima <- adjust.ts(obs)$klima
  klima.87.01 <- adjust.ts(obs)$klima.87.01

  while (dev.cur() > 1) dev.off()
  if (test) {bitmap(file="ecfcst2eps-test1a.png",type="png256"); plot(klima,type="l"); while (dev.cur() > 1) dev.off()}
  if (test) {bitmap(file="ecfcst2eps-test2a.png",type="png256"); plot(klima.87.01,type="l"); while (dev.cur() > 1) dev.off()}

#  if (verbose) {print("Seasonal values:"); print(rbind(obs.ts$yy,obs.ts$mm,round(obs.ts$value,1),klima))}
  
# Observations -> convert to anomalies wrt 1961-90:
  if (verbose) print("Estmate anomalies")
  y <- obs.ts$value[ii] - klima[ii]; yy <- obs.ts$yy[ii]; mm <- obs.ts$mm[ii]
  if (verbose) {print("Observations:"); print(summary(y)); print(length(y)); print(range(yy*100+mm))}
  now <- length(y); now2 <- length(obs.ts$value)
  print(paste("Anomaly on ",yy[now]*100+mm[now],"is",y[now]," ( total value is",round(obs.ts$value[now2],2),
                            100*obs.ts$yy[now2]+obs.ts$mm[now2],")"))
#  if (verbose) {print("Seasonal anomalies:"); print(rbind(obs.ts$yy,obs.ts$mm,round(obs.ts$value,1)))}

  # 'Add another year of climatological values'
  if (verbose) print("Extend the climatology to the future...")
  iv <- (length(obs.ts$mm)-11):length(obs.ts$mm)

  if (verbose) print(table(obs.ts$mm[iv]))
  while (sum(as.numeric(table(obs.ts$mm[iv]))>1)) {
    if (verbose) print("Detected time gaps (ecfcst2eps)!")
    iv <- iv[-1]
  } 
    
#  print("iv:"); print(iv)
  mm.klima <- c(obs.ts$mm, obs.ts$mm[iv])
  yy.klima <- c(obs.ts$yy,obs.ts$yy[iv]+1)
  if (test) {bitmap(file="ecfcst2eps-test3.png",type="png256"); plot(yy.klima+(mm.klima-0.5)/12)}

# 'Add another year of climatological values'
#  if (verbose) print("Climatology:"); print(table(yy.klima*100+mm.klima))
  klima <- c(klima,klima[iv])
  klima.87.01 <- c(klima.87.01,klima.87.01[iv])
  if (test) {points(yy.klima+(mm.klima-0.5)/12,col="red",pch=20); while (dev.cur() > 1) dev.off()}
  if (test) {bitmap(file="ecfcst2eps-test1b.png",type="png256");
             plot(klima,type="l"); while (dev.cur() > 1) dev.off()}
  if (test) {bitmap(file="ecfcst2eps-test2b.png",type="png256");
             plot(klima.87.01,type="l"); while (dev.cur() > 1) dev.off()}
  
  v <- is.element(yy.klima*100 + mm.klima, a$yy*100 + a$mm)
  vi <- is.element(a$yy*100 + a$mm,yy.klima*100 + mm.klima)
  
# For debugging purposes
#  if (verbose) {print("ts[,]:"); print(summary(c(ts))); print(summary(c(yy))); print(summary(c(mm)))}
#  if (verbose) {print("a$yy/a$mm:"); print(summary(c(a$yy))); print(summary(c(a$mm)))}
#  if (verbose) {print("a$yy/a$mm:"); print(rbind(a$yy[vi],a$mm[vi]))}
#  if (verbose) {print("Vector lengths:"); print(c(length(ts[vi,1]),length(ts[,1]),length(klima[v]),length(klima)))}
#  if (verbose) {print("yy.klima/mm.klima:"); print(summary(c(yy.klima))); print(summary(c( mm.klima)))}
#  if (verbose) {print(rbind(yy.klima[v],mm.klima[v]))}
  
  for (i in 1:n.ens) ts[vi,i] <- ts[vi,i] - klima[v] + klima.87.01[v]
  print(paste("Correction of climatology: - klima[v] + klima.87.01[v]= ",
              round(klima.87.01[v][sum(v)] -  klima[v][sum(v)],2)))
  #print(summary(c(ts)))


  if (verbose) print("Figure...")
  make.epsfig(loc,obs,obs.ts,yy,mm,FCST,res=75,type="png256",a,ts,ii2,y,
              ens.mean,mon3,label,fcst.list,months,verbose,n.ens)
  make.epsfig(loc,obs,obs.ts,yy,mm,FCST,res=300,type="png256",a,ts,ii2,y,
              ens.mean,mon3,label,fcst.list,months,verbose,n.ens)
  

 #print(rbind(y,yy,mm))
  if (verbose) print("-----------------------------------------------------------------------------")

 vname=switch(FCST$vname,"p2t"="2m temperatur (C)","t2m"="2m temperatur (C)",
                         "tp"="Nedbor","msl"="Navniva lufttrykk")

 while (dev.cur() > 1) dev.off()
   class(ensemble)<- "EC.ensemble"
   results <- list(obs=y,ens=ts,yy=a$yy,mm=a$mm,num=n.ens,location=obs$location,
                   lag=lag,vname=vname,klima.61.90= klima, klima.87.01=klima.87.01,
                   lon=obs$lon, lat=obs$lat,ensemble=ensemble)
   save(file=paste("ecfcst2eps_",obs$station,".rda",sep=""),results)
   invisible(results)
}


# Evaluation: -----------------------------------------------------------------------------------------------
  

#evaluate <- function(FCST,lag=1,loc=18700,n.ens=40,label=TRUE,verbose=FALSE,ensemble=NULL,
#                       test=FALSE,mean.from.daily=FALSE,CDF=FALSE) {
evaluate <- function(FCST,lag=1,loc=18700,n.ens=51,label=TRUE,verbose=FALSE,ensemble=NULL,
                       test=FALSE,mean.from.daily=FALSE,CDF=FALSE) {
  if (verbose) print("-----------------------------------------------------------------------------")
  if (verbose) print(paste("evaluate: lag=",lag,"  loc=",loc,"  n.ens=",n.ens,"  label=",label))
  bitmap(file="scratch.png",type="png256")
  cmon <- c("januar","februar","mars","april","mai","juni",
            "juli","august","september","oktober","november","desember")
  if (is.null(FCST)) load("getecfcst.rda")
  fcst.list <- rownames(summary(FCST))
  fcst.list <- fcst.list[grep("fcst.",fcst.list)]

  ele <- switch(FCST$vname,"t2m"="TAM","tp"="RR","slp"="POM")
  breaks <- switch(FCST$vname,"t2m"=seq(-3,3,by=1),"tp"=seq(-50,50,by=10),"slp"=seq(-50,50,by=10))
  obs <- KDVH4DS(loc,param=ele,mean.from.daily=mean.from.daily)
  obs.ts <- plotStation(obs,l.anom=FALSE,what="t",trend=FALSE,std.lev=FALSE)
  plotStation(obs,l.anom=TRUE,what="t",trend=FALSE,std.lev=FALSE)
  obs.ts$value <- seasonal.mean(obs.ts$value)
  
  if (is.null(ensemble)) ensemble <- ensemble.timeseries(FCST,fcst.list,n.ens=n.ens,obs=obs,lag=lag)
  a <- ensemble$a; ts <- ensemble$ts; months <- ensemble$months; mon3 <- ensemble$mon3
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny

  if (verbose) print(summary(a))

  # 1987-2001 Climatology - extend to the whole year (into the future)...
  if (verbose) print("Climatologies")
  ii <- is.element(obs.ts$yy*100 + obs.ts$mm, a$yy*100 + a$mm)
  ii2 <- is.element(a$yy*100 + a$mm,obs.ts$yy*100 + obs.ts$mm)
  klima <- adjust.ts(obs)$klima
  klima.87.01 <- adjust.ts(obs)$klima.87.01

# Observations -> convert to anomalies wrt 1961-90:
  if (verbose) print("Estmate anomalies")
  y <- obs.ts$value[ii] - klima[ii]; yy <- obs.ts$yy[ii]; mm <- obs.ts$mm[ii]
  if (verbose) {print("Observations:"); print(summary(y)); print(length(y)); print(range(yy*100+mm))}
  now <- length(y); now2 <- length(obs.ts$value)
  if (verbose) print(paste("Anomaly on ",yy[now]*100+mm[now],"is",y[now]," ( total value is",round(obs.ts$value[now2],2),
                            100*obs.ts$yy[now2]+obs.ts$mm[now2],")"))

  # 'Add another year of climatological values'
  if (verbose) print("Extend the climatology to the future...")
  iv <- (length(obs.ts$mm)-11):length(obs.ts$mm)

  if (verbose) print(table(obs.ts$mm[iv]))
  while (sum(as.numeric(table(obs.ts$mm[iv]))>1)) {
    if (verbose) print("Detected time gaps (ecfcst2eps)!")
    iv <- iv[-1]
  } 
    
#  print("iv:"); print(iv)
  mm.klima <- c(obs.ts$mm, obs.ts$mm[iv])
  yy.klima <- c(obs.ts$yy,obs.ts$yy[iv]+1)
  if (test) {bitmap(file="evaluate-test.png",type="png256"); plot(yy.klima+(mm.klima-0.5)/12); dev.off()}

# 'Add another year of climatological values'
#  if (verbose) print("Climatology:"); print(table(yy.klima*100+mm.klima))
  klima <- c(klima,klima[iv])
  klima.87.01 <- c(klima.87.01,klima.87.01[iv])
  
  v <- is.element(yy.klima*100 + mm.klima, a$yy*100 + a$mm)
  vi <- is.element(a$yy*100 + a$mm,yy.klima*100 + mm.klima)
  
  for (i in 1:n.ens) ts[vi,i] <- ts[vi,i] - klima[v] + klima.87.01[v]
  if (verbose) print(paste("Correction of climatology: - klima[v] + klima.87.01[v]= ",
                            round(klima.87.01[v][sum(v)] -  klima[v][sum(v)],2)))

  if (verbose) print("Figure...")
  reliab.diag(y,ts,0,res=75,type="png256",loc=loc,a=a,yy=yy,mm=mm,FCST=FCST,CDF=CDF)
  reliab.diag(y,ts,0,res=300,type="png256",loc=loc,a=a,yy=yy,mm=mm,FCST=FCST,CDF=CDF)
  if (verbose) print("-----------------------------------------------------------------------------")
}


  
# Usikkerhet: -----------------------------------------------------------------------------------------------




  
make.histfig <- function (loc,obs,obs.ts,yy,mm,FCST,res=300,type="png256",a,ts,ii2,y,ens.mean,mon3,label=TRUE,
                          fcst.list,months,verbose=FALSE,n.ens,clim.3) {
  
  cmon <- c("januar","februar","mars","april","mai","juni",
            "juli","august","september","oktober","november","desember")
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny
  year <- a$yy[nf]
  if (label) main <- paste(obs$location," dato:",cmon[months[1]],year,"til folgende",cmon[months[3]]) else
             main <- " "
  vname=switch(FCST$vname,"p2t"="2m temperatur (C)","t2m"="2m temperatur (C)",
               "tp"="Nedbor","msl"="Navniva lufttrykk")

  if (lower.case(FCST$vname)=="p2t" | lower.case(FCST$vname)=="t2m") dx<-0.5
                                                               else  dx<-10  
 
  print(paste("Variable=",lower.case(FCST$vname),vname," dx=",dx))
  x<-seq(floor(min(c(y,c(ts)),na.rm=TRUE))-dx,ceiling(max(c(y,c(ts))+dx,na.rm=TRUE)),by=dx) 
  print(range(x)); print(range(ts[nf,]))
  print(c(floor(min(c(y,c(ts)),na.rm=TRUE)),ceiling(max(c(y,c(ts)),na.rm=TRUE))))
  
  imon <- is.element(mm,months[1])
  h <- hist(ts[nf,],breaks=x)
  print(summary(ts[nf,]))
  h$density = h$density/(sum(h$density*dx))
  pdf <- dnorm(x,mean=mean(y[imon],na.rm=TRUE),sd=sd(y[imon],na.rm=TRUE))
  y.max <- max(c(h$density),pdf)
  
  while (dev.cur() > 1) dev.off()
  if (lower.case(FCST$vname)=="p2t") {xlab <- "Temperatur avvik fra normalen"; air<-0.1} else
                                     {xlab <- "Nedbor avvik fra normalen"; air <- 0.175}
  
  bitmap(paste("usikkerhet",loc,FCST$vname,"-",res,".png",sep=""),type="png256",res=300)  
  plot(range(h$mids),c(0,y.max),type="n",
      xlab="Temperatur avvik fra normalen",ylab="Sannsynlighets tetthet",
      sub=paste("Interpolert til ",obs$lat,"N/",obs$lon,"E. (obs. fordeling for ",
                min(yy)," - ",max(yy),")",sep=""),col="darkblue", main=main)
  grid()

# REB 16.10.2006  x <- seq(floor(min(c(y,c(ts)),na.rm=TRUE))-1,ceiling(max(y,c(ts),na.rm=TRUE))+1,by=0.05)
  polygon(c(x,x[1]),c(pdf,0),lwd=2,col="grey95",border="grey80")
  lines(rep(mean(y[imon],na.rm=TRUE)+0.430819*sd(y[imon],na.rm=TRUE),2),c(0,1),lty=3,col="grey60")
  lines(rep(mean(y[imon],na.rm=TRUE)-0.430819*sd(y[imon],na.rm=TRUE),2),c(0,1),lty=3,col="grey60")

  ens.mean <- mean(ts[nf,])
  lines(rep(ens.mean,2),c(y.max*0.01,y.max-y.max*0.01),lwd=3,col="grey70")
  lines(rep(ens.mean,2),c(y.max*0.01,y.max-y.max*0.01),lwd=2,col="darkblue")
  lines(h$mids,h$density,lwd=6,col="darkblue")
  lines(h$mids,h$density,lwd=4,col="blue")

  for (i in 1:n.ens) lines(rep(ts[nf,i],2),c(-1,1)*y.max*0.01) 

  if (lower.case(FCST$vname)=="p2t" | lower.case(FCST$vname)=="t2m") upper.axis<-seq(-10,10,by=2) else
                                                                     upper.axis<-seq(-150,150,by=30)

  for (x in upper.axis) {
    lines(rep(x,2),y.max +c(-1,1)*y.max*0.01,lwd=2,col="grey30")
    text(x,y.max-0.02*y.max,round(x+mean(clim.3[imon],na.rm=TRUE),1),col="grey30",font=2)
  }
  if (label) mtext(side=4,paste("(ecfcst2hist - getecfcst.R)  - ",date()))
  while (dev.cur() > 1) dev.off()
}

make.probfig <- function (loc,obs,obs.ts,yy,mm,FCST,res=300,type="png256",a,ts,ii2,y,ens.mean,mon3,label,
                          fcst.list,months,verbose=FALSE,n.ens,clim.3) {
  
  cmon <- c("januar","februar","mars","april","mai","juni",
            "juli","august","september","oktober","november","desember")
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny
  year <- a$yy[nf]
  if (label) main <- paste(obs$location," dato:",cmon[months[1]],year,"til folgende",cmon[months[3]]) else
             main <- " "
  vname=switch(FCST$vname,"p2t"="2m temperatur (C)","t2m"="2m temperatur (C)",
               "tp"="Nedbor","msl"="Navniva lufttrykk")

  if (lower.case(FCST$vname)=="p2t" | lower.case(FCST$vname)=="t2m") dx<-0.5
                                                               else  dx<-5  
 
  print(paste("Variable=",lower.case(FCST$vname),vname," dx=",dx))
  x<-seq(floor(min(c(y,c(ts)),na.rm=TRUE))-5,ceiling(max(c(y,c(ts)),na.rm=TRUE))+5,by=dx) 
  print(range(x)); print(range(ts[nf,]))
  print(c(floor(min(c(y,c(ts)),na.rm=TRUE)),ceiling(max(c(y,c(ts)),na.rm=TRUE))))

  edf <- empirical.ranking(ts[nf,])
  while (dev.cur() > 1) dev.off()
  if (lower.case(FCST$vname)=="p2t") {xlab <- "Temperatur avvik fra normalen"; air<-0.1} else
                                     {xlab <- "Nedbor avvik fra normalen"; air <- 0.1}
  
  bitmap(paste("sannsynlighet",loc,FCST$vname,"-",res,".png",sep=""),type="png256",res=300)  
  plot(range(x),c(0,1),type="n",
      xlab="Temperatur avvik fra normalen",ylab="Pr(X<x)",
      sub=paste("Interpolert til ",obs$lat,"N/",obs$lon,"E. (obs. fordeling for ",
                min(yy)," - ",max(yy),")",sep=""),col="darkblue", main=main)
  grid()

# REB 16.10.2006  x <- seq(floor(min(c(y,c(ts)),na.rm=TRUE))-1,ceiling(max(y,c(ts),na.rm=TRUE))+1,by=0.05)
  imon <- is.element(mm,months[1])
  pdf <- pnorm(x,mean=mean(y[imon],na.rm=TRUE),sd=sd(y[imon],na.rm=TRUE))
  polygon(c(x,max(x),min(x)),c(pdf,0,0),lwd=2,col="grey95",border="grey80")
  lines(rep(mean(y[imon],na.rm=TRUE)+0.430819*sd(y[imon],na.rm=TRUE),2),c(0,1),lty=3,col="grey60")
  lines(rep(mean(y[imon],na.rm=TRUE)-0.430819*sd(y[imon],na.rm=TRUE),2),c(0,1),lty=3,col="grey60")

  ens.mean <- mean(ts[nf,])
  lines(rep(ens.mean,2),c(0,1),lwd=3,col="grey70")
  lines(rep(ens.mean,2),c(0,1),lwd=2,col="darkblue")
  
  lines(edf$x,edf$P,lwd=6,col="darkblue")
  lines(edf$x,edf$P,lwd=4,col="blue")
  
  if (lower.case(FCST$vname)=="p2t" | lower.case(FCST$vname)=="t2m") {
    for (i in 1:n.ens) {
      lines(rep(ts[nf,i],2),c(-0.01,0.01))
    }
  } else {
    for (i in 1:n.ens) {
      lines(rep(ts[nf,i],2),c(-0.01,0.01))
    }
  }
  if (lower.case(FCST$vname)=="p2t" | lower.case(FCST$vname)=="t2m") upper.axis<-seq(-10,10,by=2) else
                                                                     upper.axis<-seq(-150,150,by=30)

  #print(summary(clim.3[imon]))
  for (x in upper.axis) {
    lines(rep(x,2),c(0.97,1),lwd=2,col="grey30")
    text(x,0.95,round(x+mean(clim.3[imon],na.rm=TRUE),1),col="grey30",font=2)
  }
  if (label) mtext(side=4,paste("(ecfcst2hist - getecfcst.R)  - ",date()))
  while (dev.cur() > 1) dev.off()
}



#ecfcst2hist <- function(FCST,lag=1,loc=18700,n.ens=40,empirical=FALSE,label=TRUE,verbose=FALSE,
#                        ensemble=NULL,test=FALSE,mean.from.daily=FALSE,Prob=FALSE) {
ecfcst2hist <- function(FCST,lag=1,loc=18700,n.ens=51,empirical=FALSE,label=TRUE,verbose=FALSE,
                        ensemble=NULL,test=FALSE,mean.from.daily=FALSE,Prob=FALSE) {
  
  if (verbose) print("***-----------------------------------------------------------------------***")
  if (verbose)   print(paste("ecfcst2hist: lag=",lag,"  loc=",loc,"  n.ens=",n.ens,"  label=",label))  

  bitmap("scratch.png",type="png256")
  cmon <- c("januar","februar","mars","april","mai","juni",
            "juli","august","september","oktober","november","desember")

  if (is.null(FCST)) load("getecfcst.rda")
  fcst.list <- rownames(summary(FCST))
  fcst.list <- fcst.list[grep("fcst.",fcst.list)]
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny
  ele <- switch(FCST$vname,"t2m"="TAM","tp"="RR","slp"="POM")
  obs <- KDVH4DS(loc,param=ele,,mean.from.daily=mean.from.daily)
  obs.ts <- plotStation(obs,l.anom=FALSE,what="t",trend=FALSE,std.lev=FALSE)
  plotStation(obs,l.anom=TRUE,what="t",trend=FALSE,std.lev=FALSE)
  obs.ts$value <- seasonal.mean(obs.ts$value)
  clim <- colMeans(obs$val[(obs$yy >= 1961) & (obs$yy <= 1990),])
  clim.3 <- 0.3333 * (clim + c(clim[2:12],clim[1]) + c(clim[3:12],clim[1:2]))

  if (is.null(ensemble)) ensemble <- ensemble.timeseries(FCST,fcst.list,n.ens=n.ens,obs=obs,lag=lag)
  a <- ensemble$a; ts <- ensemble$ts; months <- ensemble$months; mon3 <- ensemble$mon3
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny
  
#print(summary(a))
  ii <- is.element(obs.ts$yy*100 + obs.ts$mm, a$yy*100 + a$mm)
  ii2 <- is.element(a$yy*100 + a$mm,obs.ts$yy*100 + obs.ts$mm)
  klima <- adjust.ts(obs)$klima
  klima.87.01 <- adjust.ts(obs)$klima.87.01
  while (dev.cur() > 1) dev.off()
  if (test) {bitmap(file="ecfcst2eps-test1a.png",type="png256"); plot(klima,type="l"); while (dev.cur() > 1) dev.off()}
  if (test) {bitmap(file="ecfcst2eps-test2a.png",type="png256"); plot(klima.87.01,type="l"); while (dev.cur() > 1) dev.off()}
  
# Observations -> convert to anomalies wrt 1961-90:
  if (verbose) print("Estmate anomalies")
  y <- obs.ts$value - klima; yy <- obs.ts$yy; mm <- obs.ts$mm
  if (verbose) {print("Observations:"); print(summary(y)); print(length(y)); print(range(yy*100+mm))}
  now <- length(y); now2 <- length(obs.ts$value)
  print(paste("Anomaly on ",yy[now]*100+mm[now],"is",y[now]," ( total value is",round(obs.ts$value[now2],2),
                            100*obs.ts$yy[now2]+obs.ts$mm[now2],")"))

  # 'Add another year of climatological values'
  if (verbose) print("Extend the climatology to the future...")
  iv <- (length(obs.ts$mm)-11):length(obs.ts$mm)

  if (verbose) print(table(obs.ts$mm[iv]))
  while (sum(as.numeric(table(obs.ts$mm[iv]))>1)) {
    if (verbose) print("Detected time gaps (ecfcst2hist)!")
    iv <- iv[-1]
  }

  print("iv:"); print(iv)
  mm.klima <- c(obs.ts$mm,obs.ts$mm[iv])
  yy.klima  <- c(obs.ts$yy,obs.ts$yy[iv]+1)
  if (test) {bitmap(file="ecfcst2eps-test3.png",type="png256"); plot(yy.klima+(mm.klima-0.5)/12)}

  #print(table(mm.klima))
  #print(table(yy.klima))

  # 'Add another year of climatological values'
#  if (verbose) print("Climatology:"); print(table(yy.klima*100+mm.klima))
  klima <- c(klima,klima[iv])
  klima.87.01 <- c(klima.87.01,klima.87.01[iv])
  if (test) {points(yy.klima+(mm.klima-0.5)/12,col="red",pch=20); while (dev.cur() > 1) dev.off()}
  if (test) {bitmap(file="ecfcst2eps-test1b.png",type="png256"); plot(klima,type="l"); while (dev.cur() > 1) dev.off()}
  if (test) {bitmap(file="ecfcst2eps-test2b.png",type="png256"); plot(klima.87.01,type="l"); while (dev.cur() > 1) dev.off()}
  
  #if (verbose) print("Climatology:"); print(table(yy.klima*100+mm.klima))
  klima <- c(klima,klima[iv])
  klima.87.01 <- c(klima.87.01,klima.87.01[iv])
   
  v <- is.element(yy.klima*100 + mm.klima, a$yy*100 + a$mm)
  vi <- is.element(a$yy*100 + a$mm,yy.klima*100 + mm.klima)
  
# For debugging purposes
#  if (verbose) {print("ts[,]:"); print(summary(c(ts))); print(summary(c(yy))); print(summary(c(mm)))}
#  if (verbose) {print("a$yy/a$mm:"); print(summary(c(a$yy))); print(summary(c(a$mm)))}
#  if (verbose) {print("a$yy/a$mm:"); print(rbind(a$yy[vi],a$mm[vi]))}
#  if (verbose) {print("Vector lengths:"); print(c(length(ts[vi,1]),length(ts[,1]),length(klima[v]),length(klima)))}
#  if (verbose) {print("yy.klima/mm.klima:"); print(summary(c(yy.klima))); print(summary(c( mm.klima)))}
#  if (verbose) {print(rbind(yy.klima[v],mm.klima[v]))}

  for (i in 1:n.ens) ts[vi,i] <- ts[vi,i] - klima[v] + klima.87.01[v]
  print(paste("Correction of climatology: - klima[v] + klima.87.01[v]= ",
              round(klima.87.01[v][sum(v)] -  klima[v][sum(v)],2)))
  #print(summary(c(ts)))

  if (!Prob) {
    make.histfig(loc,obs,obs.ts,yy,mm,FCST,res=75,type="png256",a,ts,ii2,y,ens.mean,mon3,label,
                 fcst.list,months,verbose,n.ens,clim.3)
    make.histfig(loc,obs,obs.ts,yy,mm,FCST,res=300,type="png256",a,ts,ii2,y,ens.mean,mon3,label,
               fcst.list,months,verbose,n.ens,clim.3)
  } else {
    make.probfig(loc,obs,obs.ts,yy,mm,FCST,res=75,type="png256",a,ts,ii2,y,ens.mean,mon3,label,
                 fcst.list,months,verbose,n.ens,clim.3)
    make.probfig(loc,obs,obs.ts,yy,mm,FCST,res=300,type="png256",a,ts,ii2,y,ens.mean,mon3,label,
               fcst.list,months,verbose,n.ens,clim.3)
  }
  if (verbose) print("***-----------------------------------------------------------------------***")

  vname=switch(FCST$vname,"p2t"="2m temperatur (C)","t2m"="2m temperatur (C)",
                         "tp"="Nedbor","msl"="Navniva lufttrykk")

  results <- list(obs=y,ens=ts,yy=a$yy,mm=a$mm,num=n.ens,location=obs$location,lag=lag,vname=vname,
                  klima.61.90= klima, klima.87.01=klima.87.01, lon=obs$lon, lat=obs$lat, ensemble=ensemble)
  save(file=paste("ecfcst2hist_",obs$station,".rda",sep=""),results)
  invisible(results)
}




#ecfcst2map <- function(FCST,lag=1,n.ens=40,label=TRUE,res=300) {
ecfcst2map <- function(FCST,lag=1,n.ens=51,label=TRUE,res=300) {

  #print("ecfcst2map")
  while (dev.cur() > 1) dev.off()
  bitmap(file="scratch.png",type="png256")
  
  cmon <- c("Jan","Feb","Mar","Apr","May","Jun",
            "Jul","Aug","Sep","Oct","Nov","Dec")
  if (is.null(FCST)) load("getecfcst.rda")
  fcst.list <- rownames(summary(FCST))
  fcst.list <- fcst.list[grep("fcst.",fcst.list)]
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny

  dat <- rep(NA,nf*nx*ny); dim(dat) <- c(nf,ny,nx)
  ensemble <- rep(NA,n.ens*nx*ny); dim(ensemble) <- c(nx,ny,n.ens)
  good <- rep(0,nx*ny); dim(good) <- c(ny,nx)
  tim <- 1:nf; yy <- rep(NA,nf); mm <- yy; dd <- yy
  id.x <- matrix(rep(FCST$vname,ny*nx),ny,nx); id.t <- rep(FCST$vname,nf)
  #print(dim(dat))

  #print("Loop over fcsts")
  for (ifcst in 1:nf) {
      cline <- paste("FCST$",fcst.list[ifcst],"$",FCST$vname,sep="")
      #print(cline)
      yy[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$yy[",lag+1,"]",sep="")))
      mm[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$mm[",lag+1,"]",sep="")))
      dd[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$dd[",lag+1,"]",sep="")))
      #print(c(dd[ifcst],mm[ifcst],yy[ifcst]))
      fcst <- eval(parse(text=cline))
      #print(c(length(fcst),NA,dim(fcst)))
      months <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$mm[",lag+1,":",lag+3,"]",sep="")))
      #print(months)
      year <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$yy[",lag+1,"]",sep="")))
      #print(year)
   
      dat[ifcst,,] <- 0; good[,] <- 0
      ensemble[,,] <- 0
      #print(paste("Ensemble mean",ifcst))
      for (inum in 1:n.ens) {
         Map <- 0.3333* t(fcst[,,inum,lag+1] + fcst[,,inum,lag+2] + fcst[,,inum,lag+3])
         dat[ifcst,,] <- dat[ifcst,,] + Map
         ensemble[,,inum] <- 0.3333* fcst[,,inum,lag+1] + fcst[,,inum,lag+2] + fcst[,,inum,lag+3]
         #contour(FCST$lon,FCST$lat,t(Map)); addland()
         good <- good + is.finite(Map)
      }
      #good[good == 0] <- NA
      #dat[ifcst,,] <- dat[ifcst,,]/good
    }
  # dat <- dat/40
  dat <- dat/n.ens

  months <- eval(parse(text=paste("FCST$",fcst.list[nf],"$mm[",lag+1,":",lag+3,"]",sep="")))
  print("last months");  print(months)

  print("Make a field object")
  #if (FCST$vname=="p2t") dat <- dat - 273
  field  <- list(dat=dat,lon=FCST$lon,lat=FCST$lat,tim=tim,
                 v.name=FCST$vname,id.x=id.x,id.t=id.t,
                 yy=yy,mm=mm,dd=dd,n.fld=1,yy0=min(yy),
                 id.lon=rep(FCST$vname,FCST$nx),id.lat=rep(FCST$vname,FCST$ny),
                 attributes=NULL)
  class(field) <- c("field","monthly.field.object")
  attr(field$tim,"units") <- "month"

  #print("Get ERA40")
  print(lower.case(FCST$vname))
  elem <- switch(lower.case(FCST$vname),"p2t"="t2m","t2m"="t2m","tp"="prec","msl"="slp")
  levels <- switch(lower.case(FCST$vname),"p2t"=seq(-4,4,by=0.5),"t2m"=seq(-4,4,by=0.5),
                                          "tp"=seq(-30,30,by=2))
  #print(paste("ERA40_",elem,"_mon.Rdata",sep=""))
  load(paste("ERA40_",elem,"_mon.rda",sep=""))
  ERA40 <- switch(FCST$vname,"p2t"=t2m,"t2m"=t2m,"tp"=prec,"msl"=slp)
  iiyy <- is.element(field$yy,ERA40$yy)
  if ((sum(iiyy))==0) {
    print(range(field$yy))
    print(range(ERA40$yy))
    stop("getecfcst2map: ERROR non-overlapping years!")
  }
  iyy <- max((1:length(field$yy))[iiyy]); lyy <- field$yy[iyy]; lmm <-  field$mm[iyy]
  iyy <- min((1:length(field$yy))[iiyy]); fyy <- field$yy[iyy]; fmm <-  field$mm[iyy]
  print(c(range(field$yy),NA,range(ERA40$yy),NA,lyy,fyy))
  print(paste("1-",cmon[fmm],"-",1961," &     1-",cmon[lmm],"-",1990,sep=""))

# Anomali: x'(1987--2001) = x - clim(1987--2001)
#          x'(1961--1990) = x - clim(1987--2001) + clim(1987--2001) - clim(1961--1990)
#          x'(1961--1990) = x'(1987--2001) + clim(1987--2001) - clim(1961--1990)
#
# adjust = clim(1961--1990) - clim(1987--2001)
#
#          x'(1961--1990) = x'(1987--2001) - adjust
  dev.off()

#  print("ecfcst2map.era40_61-90.png")
  bitmap("ecfcst2map.era40_61-90.png",type="png256",res=res)
  mean.era40.61.90  <- meanField(ERA40,t.rng=c(paste("1-",cmon[fmm],"-",1961,sep=""),
                                               paste("1-",cmon[lmm],"-",1990,sep="")),
                           lon.rng=range(field$lon),lat.rng=range(field$lat),mon=months)
  map(mean.era40.61.90,main="ERA40: 1961-90 climatology",newFig=FALSE)
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  print("ecfcst2map.era40_87-01.png")
  bitmap("ecfcst2map.era40_87-01.png",type="png256",res=res)
  mean.era40.87.01  <- meanField(ERA40,t.rng=c(paste("1-",cmon[fmm],"-",1987,sep=""),
                                               paste("1-",cmon[lmm],"-",2001,sep="")),
                           lon.rng=range(field$lon),lat.rng=range(field$lat),mon=months)
  map(mean.era40.87.01,main="ERA40: 1987-2001 climatology",newFig=FALSE)
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  print("Adjustment")
#           adjust = clim(1961--1990) - clim(1987--2001):

  bitmap("ecfcst2map.adjust.png",type="png256",res=res)
  adjust <- map(mean.era40.61.90,mean.era40.87.01,main="ERA40 (1961--1990) - (1987--2001) adjustment",newFig=FALSE)
  #map(adjust,main="ERA40: (1961--1990) - (1987--2001) adjustment")
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  #print("ecfcst2map.tot_field.png")
  unit.scaling <- 3600*24*30 * 1000  # Needs changing to 28, 29, 30 or 31 days... XXX
  if (lower.case(FCST$vname) == "tp") field$dat <- field$dat * unit.scaling
  bitmap("ecfcst2map.tot_field.png",type="png256",res=res)
  fcst.tot <- mapField(field,what="abs",plot=FALSE)
  print(summary(c(fcst.tot$map)))
  map(fcst.tot,main="ECMWF original seasonal forecast",levels=levels,newFig=FALSE) 
  print("Construct climatology map")
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  print("ecfcst2map.clim.png")
  bitmap("ecfcst2map.clim.png",type="png256",res=res)
  clim <- meanField(ERA40,t.rng=c("1-Jan-1961","31-Dec-1990"),
                    lon.rng=range(field$lon),lat.rng=range(field$lat),mon=months)
  map(clim,main="ERA40: 1961-90 climatology",newFig=FALSE)
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  print("Construct adjusted anomaly map")
  print(paste("kart_",FCST$vname,".png",sep=""))
  bitmap(paste("kart_",FCST$vname,".png",sep=""),type="png256",res=res) 
  print(dev.cur())
  par(cex.main=0.8)
#          x'(1961--1990) = x'(1987--2001) - adjust:
  #print(options()$device); print(levels)
  inv.col <- switch(lower.case(FCST$vname),"tp"=TRUE,"t2m"=TRUE,otherwise=FALSE)
  fcst <- map(fcst.tot,adjust,main="avvik fra ECMWF klimatologi",
              sub=paste(cmon[months[1]],year,"til folgende",cmon[months[3]]),
              levels=levels,inv.col=inv.col,newFig=FALSE)
  print(dev.cur())
  mtext(side=4,paste("(ecfcst2map - getecfcst.R)"))
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  vname <- switch(FCST$vname,"p2t"="t2m","t2m"="t2m","tp"="precip","msl"="slp")

  longname <- switch(FCST$vname,"p2t"="2-meter temperature anomaly from 1961-90 mean (adjusted)",
                                "t2m"="2-meter temperature anomaly from 1961-90 mean (adjusted)",
                                "tp"="total precipitation",
                                "msl"="sea level pressure")
  unit  <- switch(FCST$vname,"p2t"="deg C","t2m"="deg C","tp"="mm","msl"="hPa")

  imon <- is.element(ERA40$mm,months)
  yrs <- as.numeric(rownames(table(ERA40$yy[imon]))); nyrs <- length(yrs)
  ix <- (ERA40$lon >= min(field$lon)) & (ERA40$lon <= max(field$lon)) 
  i1 <- (1:length(ERA40$lon))[ix][1]-1
  iy <- (ERA40$lat >= min(field$lat)) & (ERA40$lat <= max(field$lat))
  j1 <- (1:length(ERA40$lat))[iy][1]-1
  nx.era40 <- length(ERA40$lon[ix]); ny.era40 <- length(ERA40$lat[iy])
  t2m.era40 <- rep(0,nx.era40*ny.era40*nyrs); dim(t2m.era40) <- c(nx.era40,ny.era40,nyrs)
  T2M.era40 <- rep(0,nx*ny*nyrs); dim(T2M.era40) <- c(nx,ny,nyrs)

# Slow, but speed is not essential here and it's important to be sure about dimensions, etc.
  print(paste("Estimate mean for given months for each year: x= ",i1+1," - ",i1+nx.era40,
              ", jj= ", j1+1," - ", j1+ny.era40,sep=""))
  for (it in 1:nyrs) {
    for (i in 1:nx.era40) {
      ii <- i + i1
      for (j in 1:ny.era40) {
        jj <- j + j1
        iyr.imon <- is.element(ERA40$yy,yrs[it]) & imon
        if (sum(iyr.imon)==0) print(paste("No data: ",ERA40$yy[iyr.imon][1]))
        t2m.era40[i,j,it] <- mean(ERA40$dat[iyr.imon,jj,ii],na.rm=TRUE)
      }
    }
  }

  nx <- length(field$lon); ny <- length(field$lat)
  mask <- rep(0,nx*ny); dim(mask) <- c(nx,ny)
  print(c(nx,ny,NA,nx.era40,ny.era40))
  print(c(range(field$lon),NA,range(ERA40$lon[ix])))
  print(c(range(field$lat),NA,range(ERA40$lat[iy])))
  lon.xy <- rep(ERA40$lon[ix],ny.era40)
  lat.xy <- sort(rep(ERA40$lat[iy],nx.era40))
  t2m.era40[is.na(t2m.era40)] <- -999
  for (it in 1:nyrs) {
    T2M.era40[,,it]  <- interp(lon.xy,lat.xy,t2m.era40[,,it],fcst$lon,fcst$lat)$z
  }
  print(summary(c(ERA40$dat[imon,,])))
  print(summary(c(t2m.era40)))
  print(summary(c(T2M.era40)))
  rm(t2m.era40)


  fcst$map[is.na(fcst$map)] <- -999
  adjust$map[is.na(adjust$map)] <- -999
  clim$map[is.na(clim$map)] <- -999
  ensemble[is.na(ensemble)] <- -999

  nx <- length(adjust$lon); ny <- length(adjust$lat)
  lon.xy <- rep(adjust$lon,ny)
  lat.xy <- sort(rep(adjust$lat,nx))
  adj <- interp(lon.xy,lat.xy,adjust$map,fcst$lon,fcst$lat)$z
  adj[abs(adj) > 100] <- -999
  for (i in 1:n.ens) ensemble[,,i] <- ensemble[,,i] - adj
  nx <- length(clim$lon); ny <- length(clim$lat)
  lon.xy <- rep(clim$lon,ny)
  lat.xy <- sort(rep(clim$lat,nx))
  clim.map <- interp(lon.xy,lat.xy,clim$map,fcst$lon,fcst$lat)$z

#x11()
#contour(fcst$lon,fcst$lat,clim.map)
#addland()

  nx <- length(fcst$lon); ny <- length(fcst$lat)
  print("Estimate probabilities of likelihood that n members belong to a given category:")
  probs <- rep(NA,3*nx*ny); dim(probs) <- c(nx,ny,3); counts <- probs; signf <- probs
  tem.cat <- rep(NA,2*nx*ny); dim(tem.cat) <- c(nx,ny,2)
  ii <- 0
  e.count <- round(sum(1:n.ens*dbinom(1:n.ens, n.ens, 0.3333)))
  for (i in 1:nx) {
    for (j in 1:ny) {
      tem.cat[i,j,1] <- round(mean(T2M.era40[i,j,],na.rm=TRUE) - 0.430819 *sd(T2M.era40[i,j,],na.rm=TRUE) -
                               clim.map[i,j],4)
      tem.cat[i,j,2] <- round(mean(T2M.era40[i,j,],na.rm=TRUE) + 0.430819 *sd(T2M.era40[i,j,],na.rm=TRUE) -
                               clim.map[i,j],4)
# FOR testing:
#      tem.cat[i,j,1] <- round(mean(T2M.era40[i,j,],na.rm=TRUE) - 0.430819 *sd(T2M.era40[i,j,],na.rm=TRUE),4)
#      tem.cat[i,j,2] <- round(mean(T2M.era40[i,j,],na.rm=TRUE) + 0.430819 *sd(T2M.era40[i,j,],na.rm=TRUE),4)
      n1 <- sum(ensemble[i,j,] <= tem.cat[i,j,1],na.rm=TRUE)
      n2 <- sum((ensemble[i,j,] > tem.cat[i,j,1]) &
                (ensemble[i,j,] < tem.cat[i,j,2]),na.rm=TRUE)
      n3 <- sum(ensemble[i,j,] >= tem.cat[i,j,2],na.rm=TRUE)
      signf[i,j,1] <- round(binom.test(n1, n.ens, 0.3333)$p.value,3)
      signf[i,j,2] <- round(binom.test(n2, n.ens, 0.3333)$p.value,3)
      signf[i,j,3] <- round(binom.test(n3, n.ens, 0.3333)$p.value,3)
      if (sum(is.finite(ensemble[i,j,])) == n.ens)
         mask[i,j] <-  1 - as.numeric((t.test(c(ensemble[i,j,]))$p.value < 0.05) &
                       ( (signf[i,j,1] < 0.05) | (signf[i,j,3] < 0.05) ))
      else mask[i,j] <- NA
      counts[i,j,1] <- n1;  counts[i,j,2] <- n2;  counts[i,j,3] <- n3 
      probs[i,j,1] <- round(100*n1/n.ens);
      probs[i,j,2] <- round(100*n2/n.ens);probs[i,j,3] <- round(100*n3/n.ens)
      ii <- ii + 1
      if (mod(ii,50)==0) print(paste("Example counts: ",n1,"+", n2,"+",n3,
              "= ",n1+n2+n3," Giving probs: ",
              signf[i,j,1],", ", signf[i,j,2],", and ", signf[i,j,3]," at ",fcst$lon[i],"E,",
              fcst$lat[j],"N. tem.cat 1&2 =", tem.cat[i,j,1]," & ",tem.cat[i,j,2],sep=""))

    }
  }

  #print(paste("options()$device=",options()$device))
  bitmap("ecfcst2map.thresh.png",type="png256",res=res)
  contour(fcst$lon,fcst$lat,tem.cat[,,1],main="Thresholds",col="blue",lwd=2)
  addland()
  contour(fcst$lon,fcst$lat,tem.cat[,,2],add=TRUE,col="red",lwd=2)
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  bitmap("ecfcst2map.prob_low.png",type="png256",res=res)
  contour(fcst$lon,fcst$lat,probs[,,1],main="Probability for cold",col="blue",lwd=2)
  addland()
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  bitmap("ecfcst2map.prob_med.png",type="png256",res=res)
  contour(fcst$lon,fcst$lat,probs[,,2],main="Probability for medium",col="blue",lwd=2)
  addland()
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  bitmap("ecfcst2map.prob_high.png",type="png256",res=res)
  contour(fcst$lon,fcst$lat,probs[,,3],main="Probability for warm",col="blue",lwd=2)
  addland()
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  bitmap("ecfcst2map.count.png",type="png256",res=res)
  contour(fcst$lon,fcst$lat,counts[,,1],main="Counts",col="blue",lwd=2)
  addland()
  contour(fcst$lon,fcst$lat,counts[,,2],add=TRUE,lwd=2)
  contour(fcst$lon,fcst$lat,counts[,,3],add=TRUE,col="red",lwd=2)
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  bitmap("ecfcst2map.sign_low.png",type="png256",res=res)
  contour(fcst$lon,fcst$lat,probs[,,1],main="Category cold p-value",sub="Binomial distr.",col="blue",lwd=2)
  addland()
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  bitmap("ecfcst2map.sign_med.png",type="png256",res=res)
  contour(fcst$lon,fcst$lat,probs[,,2],main="Category medium p-value",sub="Binomial distr.",col="blue",lwd=2)
  addland()
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  bitmap("ecfcst2map.sign_high.png",type="png256",res=res)
  contour(fcst$lon,fcst$lat,probs[,,3],main="Category warm p-value",sub="Binomial distr.",col="blue",lwd=2)
  addland()
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()


  bitmap("ecfcst2map.mask.png",type="png256",res=res)
  image(fcst$lon,fcst$lat,mask,main="Significance mask",col="blue",lwd=2)
  addland()
  if (label) mtext(side=4,date(),cex=0,7,col="grey")
  dev.off()

  signf[is.na(signf)] <- -999
  probs[is.na(probs)] <- -999
  counts[is.na(counts)] <- -999
  Map[is.na(Map)] <- -99

  save(file="ecfcst2map.rda",clim,fcst,adjust,field,fcst.tot,ensemble,probs)

  print("define the nc-dimensions")
  dimlon <- dim.def.ncdf( "lon", "deg E", fcst$lon)
  dimlat <- dim.def.ncdf( "lat", "deg N", fcst$lat)
  dimmem <- dim.def.ncdf( "member", "number", 1:n.ens)
  dimcat <- dim.def.ncdf( "category", "number", 1:3)
  dimlon.a <- dim.def.ncdf( "lon.adj", "deg E", adjust$lon)
  dimlat.a <- dim.def.ncdf( "lat.adj", "deg N", adjust$lat)
  dimlon.c <- dim.def.ncdf( "lon.clim", "deg E", clim$lon)
  dimlat.c <- dim.def.ncdf( "lat.clim", "deg N", clim$lat)
  print("define the nc-variables")
  varfcst <- var.def.ncdf(name=vname,units=unit,dim=list(dimlon,dimlat), missval=-99900, 
        longname=longname, prec="integer")  
  varfcsttot <- var.def.ncdf(name="total",units=unit,dim=list(dimlon,dimlat), missval=-99900, 
        longname=paste("total",vname,"values"), prec="integer")  
  varadjust <- var.def.ncdf(name="adjust",units=unit,dim=list(dimlon.a,dimlat.a),  missval=-99900, 
        longname="Difference in mean ERA40 - ensemble mean forecast", prec="integer")  
  varclim <- var.def.ncdf(name="clim",units=unit,dim=list(dimlon.c,dimlat.c),  missval=-99900, 
        longname="1961-90 climatology", prec="integer")  

  varens <- var.def.ncdf(name='ensemble',units=unit,dim=list(dimlon,dimlat,dimmem),  missval=-99900, 
        longname="Ensemble members", prec="integer")  
  varsignf <- var.def.ncdf(name='Signficance',units="%",dim=list(dimlon,dimlat,dimcat),  missval=-99900, 
        longname="Probability that number of members in category is random", prec="integer")  
  varprob <- var.def.ncdf(name='probabilities',units="%",dim=list(dimlon,dimlat,dimcat),  missval=-99900, 
        longname="Probability for outcome in a given category", prec="integer")  
  varcount <- var.def.ncdf(name='counts',units="%",dim=list(dimlon,dimlat,dimcat),  missval=-99900, 
        longname="Number of members in categories: cold,normal,warm", prec="integer")  
  varmask <- var.def.ncdf(name='mask',units="[0,1]",dim=list(dimlon,dimlat),  missval=-99, 
        longname="5% significance mask: t-test & binomial distrib. (0=signf., 1=not.signf.)", prec="short")  


# Create a netCDF file with this variable
  #print("create the nc-file")
  ncnew <- create.ncdf( "ecfcst.nc", list(varfcst,varfcsttot,varadjust,varclim,
                                          varprob,varsignf,varens,varcount,varmask))
  #print("add forecast")
  put.var.ncdf( ncnew, varfcst,fcst$map*100)
  #print("add total forecast")
  put.var.ncdf( ncnew, varfcsttot,fcst.tot$map*100)
  #print("add adjustment data")
  put.var.ncdf( ncnew, varadjust,adjust$map*100)
  #print("add clim data")
  put.var.ncdf( ncnew, varclim,clim$map*100)

  #print("add probabilities")
  put.var.ncdf( ncnew, varsignf, signf*100)
  put.var.ncdf( ncnew, varprob, probs*1000)
  put.var.ncdf( ncnew, varcount, counts)
  #print("add ensembles")
  #print(dim(ensemble))
  #print(c(length(fcst$lon),length(fcst$lat),n.ens))
  put.var.ncdf( ncnew, varens,ensemble*100)
  put.var.ncdf( ncnew, varmask,mask)

  print("add attribute")
  att.put.ncdf( ncnew, 0, 'dato', paste(cmon[months[1]],year,"til folgende",cmon[months[3]]))
  att.put.ncdf( ncnew, varfcst, 'scale_factor', 100)
  att.put.ncdf( ncnew, varfcsttot, 'scale_factor', 100)
  att.put.ncdf( ncnew, varadjust, 'scale_factor', 100)
  att.put.ncdf( ncnew, varclim, 'scale_factor', 100)
  att.put.ncdf( ncnew, varens, 'scale_factor', 100)
  att.put.ncdf( ncnew, varsignf, 'scale_factor', 100)
  att.put.ncdf( ncnew, varprob, 'scale_factor', 100)
  close.ncdf(ncnew)

  invisible(field)
}


#eccheck <- function(FCST,loc=18700,n.ens=40,label=TRUE,res=300)  {
eccheck <- function(FCST,loc=18700,n.ens=51,label=TRUE,res=300)  {

  cmon <- c("j","f","m","a","m","j","j","a","s","o","n","d")
  if (is.null(FCST)) load("getecfcst.rda")
  fcst.list <- rownames(summary(FCST))
  fcst.list <- fcst.list[grep("fcst.",fcst.list)]
  nf <- length(fcst.list); nx <- FCST$nx; ny <- FCST$ny
  ele <- switch(FCST$vname,"p2t"="TAM","tp"="RR","msl"="POM")
  obs <- KDVH4DS(loc,param=ele)
  obs.ts <- plotStation(obs,l.anom=FALSE,what="t",trend=FALSE,std.lev=FALSE)
  #clim <- colMeans(obs$val[(obs$yy >= 1961) & (obs$yy <= 1990),])

  dat <- rep(NA,nf*nx*ny); dim(dat) <- c(nf,ny,nx)
  tim <- 1:nf; yy <- rep(NA,nf); mm <- yy; dd <- yy
  id.x <- matrix(rep(FCST$vname,ny*nx),ny,nx); id.t <- rep(FCST$vname,nf)
  #print(dim(dat))
  ts <- rep(NA,nf*n.ens*5); dim(ts) <- c(5,nf,n.ens)
  a.mm <- rep(NA,nf*5);  dim(a.mm) <- c(5,nf)
  mon3 <- rep(" ",nf)

  for (lag in 1:5) {
   for (inum in 1:n.ens) {
    for (ifcst in 1:nf) {
      #print(c(ifcst,inum))
      cline <- paste("FCST$",fcst.list[ifcst],"$",FCST$vname,sep="")
      yy[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$yy[",lag,"]",sep="")))
      mm[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$mm[",lag,"]",sep="")))
      dd[ifcst] <- eval(parse(text=paste("FCST$",fcst.list[ifcst],"$dd[",lag,"]",sep="")))
      #print(cline)
      fcst <- eval(parse(text=cline))
      #print(c(length(fcst),NA,dim(fcst)))
      #print(c(inum,lag,lag+2,NA,dim(t(fcst[,,inum,lag])),NA,dim(dat[ifcst,,])))
      dat[ifcst,,] <- t(fcst[,,inum,lag])
    }  
 
    if (FCST$vname=="p2t") dat <- dat - 273
    field  <- list(dat=dat,lon=FCST$lon,lat=FCST$lat,tim=tim,
                   v.name=FCST$vname,id.x=id.x,id.t=id.t,
                   yy=yy,mm=mm,dd=dd,n.fld=1,
                   id.lon=rep(FCST$vname,FCST$nx),id.lat=rep(FCST$vname,FCST$ny),
                   attributes=NULL)
    class(field) <- c("field",FCST$obj.type)

    st <- plotField(field,lon=obs$lon,lat=obs$lat,what="abs",col="grey40",add=TRUE)
    a <- plotStation(st,l.anom=FALSE,what="t",add=TRUE,trend=FALSE,std.lev=FALSE,col="red")
    print(c(length(a$value),length(a$mm),NA,length(a.mm[lag,])))
    ts[lag,,inum] <- a$value
    a.mm[lag,] <- c(a$mm)
   }
  }

  
  vname=switch(FCST$vname,"p2t"="2m temperatur (C)","tp"="Nedbor","msl"="Navniva lufttrykk")
  par(col.axis="white")
  y <- obs.ts$value[is.finite(obs.ts$value)]
  plot(c(1,24), range(c(y,c(ts))),type="n",main=obs$location,xlab="Dato",ylab=vname,
       sub=paste("Interpolert til ",obs$lat,"N/",obs$lon," (monthly means)",sep=""))
  par(col.axis="black")
  axis(1,at=1:24,label=rep(cmon,2))
  axis(2)
  
  grid()
  for (i in 1:n.ens) {
    for (lag in 1:5) points(c(a.mm[lag,],a.mm[lag,]+12),rep(c(ts[lag,,i]),2),pch=20,cex=0.5,col="grey30")
  }
  points(c(obs.ts$mm,obs.ts$mm+12),rep(obs.ts$value,2),pch=5,cex=0.75)

  mtext(side=4,paste("(eccheck - getecfcst.R)"))
  
  #print(range(ts))
  dev2bitmap(paste("feil",loc,".png",sep=""),type="png256",res=res)
}

ec.hitrate <- function(x,plot=TRUE,potent=FALSE,label=TRUE,type="png256"){
  
obs.all <- x$obs
qtiles <- qnorm(seq(0.2,0.8,by=0.2))

levs      <-c( -Inf,
                mean(obs.all, na.rm=TRUE) + qtiles[1]*sd(obs.all, na.rm=TRUE),
                mean(obs.all, na.rm=TRUE) + qtiles[2]*sd(obs.all, na.rm=TRUE),
                mean(obs.all, na.rm=TRUE) + qtiles[3]*sd(obs.all, na.rm=TRUE),
                mean(obs.all, na.rm=TRUE) + qtiles[4]*sd(obs.all, na.rm=TRUE),
                Inf)

n.cat <- length(levs)-1
# Set up matrices:

categ<-c("Very high","High","Normal","Low","Very low")
n.cat<-length(categ)

print("Score matrices")
cases<-matrix(rep(0,n.cat^2),n.cat,n.cat)

print("hit-miss table for regression model")

t.o <- x$obs

# test: t.o<- x$ens[,1]; for (i in 2:40) t.o <- c(t.o,x$ens[,i])

t.p<- x$ens

dims <- dim(t.p); print(dims); print(length(t.o))

for (i in seq(2,length(levs),by=1)) {
      in.class<- (t.o > levs[i-1]) & (t.o <= levs[i])
  for (ii in seq(2,length(levs),by=1)) {
    cases[ii-1,i-1] <- sum(t.p[in.class,] > levs[ii-1] &
                           t.p[in.class,] <= levs[ii],na.rm=T)
    }
  }

cases <- round(100*cases/sum(c(cases)))

# Estimate scores:
print("estimate scores")
p.hit<-round(sum(diag(cases))/sum(cases)*100)/100

 xy.pos<-seq(1,n.cat,by=1)

  if (plot) {
    print(paste("Plotting: ",paste("hitrate_",strip(x$location),x$vname,"_",x$lag,".",substr(type,1,3),sep="")))
    
    bitmap(paste("hitrate_",strip(x$location),x$vname,"_",x$lag,".",substr(type,1,3),sep=""),type=type,res=300) 
    par(col.axis="white")
    plot(c(0,n.cat),c(0,n.cat),type="n",
       main=paste("Hit-Miss table:",x$obs.name),
       sub=paste("[p-hit=",as.character(p.hit),", expected for no skill p-val=",
       as.character(1/n.cat),"]"),
       xlab="Observed",ylab="Predicted")
    par(col.axis="black",ps=10)
 
    for (i in seq(1,n.cat-1,by=1)) {
      lines(c(0,xy.pos[n.cat]),c(xy.pos[i],xy.pos[i]),col="red")
      lines(c(xy.pos[i],xy.pos[i]),c(0,xy.pos[n.cat]),col="red")
    }

    for (i in seq(1,n.cat,by=1)) {
      for (ii in seq(1,n.cat,by=1)) {
        if (i==ii) {
          text(xy.pos[i]-0.5,xy.pos[ii]-0.5,as.character(cases[ii,i]),cex=2)
        } else {
          text(xy.pos[i]-0.5,xy.pos[ii]-0.5,as.character(cases[ii,i]),cex=2,
               col="grey40")
        }
      }
    }
    axis(1,xy.pos-0.5,categ)     
    axis(2,xy.pos-0.5,categ)

    mtext(paste(x$location," (lead=",as.character(x$lag),"): ",
                as.character(x$num)," forecasts based on a sample of ",
                as.character(dims[1]), "(",x$vname,")"),side=4,cex=0.85)
    if (potent)  mtext("ensemble v.s. ensemble mean",side=2,cex=0.85,col="darkblue")
    dev.off()
  }

eval.res <- list(cases=cases,levs=levs,p.hit=p.hit)
invisible(eval.res)
}


met.no.graphics <- function(FCST=NULL,label=TRUE,n.ens=40) {
if (is.null(FCST)) load("getecfcst.rda")
#old.dev <- options()$device
#options(device="none")

ecfcst2map(FCST,label=label,n.ens=n.ens)

ensemble <- ecfcst2eps(FCST,loc=98550,label=label)$ensemble
ecfcst2hist(FCST,loc=98550,label=label,ensemble=ensemble,n.ens=n.ens)
ecfcst2hist(FCST,loc=98550,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=98550,FCST)

ensemble <- ecfcst2eps(FCST,loc=90450,label=label)$ensemble
ecfcst2hist(FCST,loc=90450,label=label,ensemble=ensemble,n.ens=n.ens)
ecfcst2hist(FCST,loc=90450,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=90450,FCST)

ensemble <- ecfcst2eps(FCST,loc=82290,label=label,n.ens=n.ens)$ensemble
ecfcst2hist(FCST,loc=82290,label=label,ensemble=ensemble,n.ens=n.ens)
ecfcst2hist(FCST,loc=82290,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=82290,FCST)

ensemble <- ecfcst2eps(FCST,loc=69100,label=label,n.ens=n.ens)$ensemble
ecfcst2hist(FCST,loc=69100,label=label,ensemble=ensemble,n.ens=n.ens)
ecfcst2hist(FCST,loc=69100,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=69100,FCST)

ensemble <- ecfcst2eps(FCST,loc=50540,label=label,n.ens=n.ens)$ensemble
ecfcst2hist(FCST,loc=50540,label=label,ensemble=ensemble,n.ens=n.ens)
ecfcst2hist(FCST,loc=50540,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=50540,FCST)

ensemble <- ecfcst2eps(FCST,loc=44560,label=label,n.ens=n.ens)$ensemble
ecfcst2hist(FCST,loc=44560,label=label,ensemble=ensemble,n.ens=n.ens)
ecfcst2hist(FCST,loc=44560,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=44560,FCST)

ensemble <- ecfcst2eps(FCST,loc=39040,label=label,n.ens=n.ens)$ensemble
ecfcst2hist(FCST,loc=39040,label=label,ensemble=ensemble,n.ens=n.ens)
ecfcst2hist(FCST,loc=39040,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=39040,FCST)

ensemble <- ecfcst2eps(FCST,loc=18700,label=label,n.ens=n.ens)$ensemble
ecfcst2hist(FCST,loc=18700,label=label,ensemble=ensemble,n.ens=n.ens)
ecfcst2hist(FCST,loc=18700,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=18700,FCST)

ensemble <- ecfcst2eps(FCST,loc=07010,label=label,n.ens=n.ens)$ensemble
ecfcst2hist(FCST,loc=07010,label=label,ensemble=ensembl,n.ens=n.ense)
ecfcst2hist(FCST,loc=07010,label=label,ensemble=ensemble,Prob=TRUE,n.ens=n.ens)
evaluate(loc=07010,FCST)

#eccheck(FCST,loc=98550)
#eccheck(FCST,loc=90450)
#eccheck(FCST,loc=82290)
#eccheck(FCST,loc=69100)
#eccheck(FCST,loc=50540)
#eccheck(FCST,loc=44560)
#eccheck(FCST,loc=39040)
#eccheck(FCST,loc=18700)
#eccheck(FCST,loc=07010)

options(device=old.dev)
}

check.point <- function(FCST=NULL,loc=18700,label=TRUE) {
  options(device="x11")
  print("============= Diagnosing the results ================")
  print("check.point()")

  if (is.null(FCST)) load("getecfcst.rda")
  options(device="none")
  eps.fcst <- ecfcst2eps(FCST,loc=loc)
  unc.fcst <- ecfcst2hist(FCST,loc=loc)
  print(paste(eps.fcst$location,": loc=",loc,sep=""))

  load("ecfcst2map.rda")
  #print(c(nx,ny)); print(fcst$lon); print(fcst$lat)
#   save(file="ecfcst2map.Rdata",clim,fcst,adjust,field,fcst.tot)
#   results <- list(obs=y,ens=ts[ii2,],yy=a$yy,mm=a$mm,num=n.ens,location=obs$location,lag=lag,vname=vname,
#                  klima.61.90= klima, klima.87.01=klima.87.01, lon=obs$lon, lat=obs$lat)
  fcst$map[!is.finite(fcst$map)] <- -999
  #print("fcst:interp")
  nx <- length(fcst$lon); ny <- length(fcst$lat)
  lon.xy <- rep(fcst$lon,ny)
  lat.xy <- sort(rep(fcst$lat,nx))
  map.val <- interp(lon.xy,lat.xy,fcst$map,   eps.fcst$lon,eps.fcst$lat)$z
  tot.val <- interp(lon.xy,lat.xy,fcst.tot$map,   eps.fcst$lon,eps.fcst$lat)$z
  adjust$map[!is.finite(adjust$map)] <- -999
  #print("adjust:interp")
  nx <- length(adjust$lon); ny <- length(adjust$lat)
  lon.xy <- rep(adjust$lon,ny)
  lat.xy <- sort(rep(adjust$lat,nx))
  map.adj <- interp(lon.xy,lat.xy,adjust$map, eps.fcst$lon,eps.fcst$lat)$z
  clim$map[!is.finite(clim$map)] <- -999
  #print("clim:interp")
  nx <- length(clim$lon); ny <- length(clim$lat)
  lon.xy <- rep(clim$lon,ny)
  lat.xy <- sort(rep(clim$lat,nx))
  map.cli <- interp(lon.xy,lat.xy,clim$map,   eps.fcst$lon,eps.fcst$lat)$z

  nt <- length(eps.fcst$yy)
  print(c(dim(eps.fcst$ens),NA,dim(unc.fcst$ens)))
  print("Date & ensemble mean values:")
  print(c(length(rowMeans(eps.fcst$ens)),length(colMeans(eps.fcst$ens))))
  print(rbind(eps.fcst$yy,eps.fcst$mm,rowMeans(eps.fcst$ens)))
  print(paste("nt=",nt))
  print("most recent forecast:")
  print(dim(eps.fcst$ens))
  print(eps.fcst$ens[nt,])
  #print(unc.fcst$ens[nt,])
  print("---------------------------------------------------------------")

  print(paste("ecfcst2eps:  tot= ",round(mean(eps.fcst$ens[,nt]),2),
              " year/month= ",eps.fcst$yy[nt],"/",eps.fcst$mm[nt],
              " klima.61.90 (",-round(eps.fcst$klima.61.90[nt],2),
              ")  + klima.87.01 (",round(eps.fcst$klima.87.01[nt],2),
              ") = adjust =",round(eps.fcst$klima.87.01[nt]-eps.fcst$klima.61.90[nt],2)))

  print(paste("ecfcst2hist: tot= ",round(mean(unc.fcst$ens[,nt]),2),
              " year/month= ",unc.fcst$yy[nt],"/",unc.fcst$mm[nt],
              " klima.61.90=",round(unc.fcst$klima.61.90[nt],2),
              " klima.87.01=",round(unc.fcst$klima.87.01[nt],2),
              " adjust=",round(unc.fcst$klima.87.01[nt]-unc.fcst$klima.61.90[nt],2)))

  print(paste("ecfcst2map: tot=",round(tot.val,2)," + adjust=",-round(map.adj,2)," = fcst.adjusted=",
              round(map.val,2)))

  print(paste("Original ecfcst2eps forecast= ",
        round(mean(eps.fcst$ens[nt,]) + eps.fcst$klima.61.90[nt] - eps.fcst$klima.87.01[nt],2)))

  print("---------------------------------------------------------------")


  bitmap(file="check.point.png",type="png256")
  fcst$map[abs(fcst$map) > 100] <- NA
  contour(fcst$lon,fcst$lat,fcst$map,lwd=2,col="darkblue")
  grid()
  addland()
  contour(fcst$lon,fcst$lat,fcst$map,lwd=2,col="darkblue",add=TRUE)
  points(eps.fcst$lon,eps.fcst$lat,pch=20,col="grey80",cex=1.5)
  text(eps.fcst$lon,eps.fcst$lat,round(map.val,2),font=2,cex=1.5)
  if (label) mtext(side=1,date(),cex=0,7,col="grey")
  dev.off()
}

check.field <- function(FCST=NULL,label=TRUE,n.ens=51) {
  if (is.null(FCST)) load("getecfcst.rda")
  flds <- rownames(summary(FCST)); flds <- flds[grep("fcst.",flds)]
  options(device="x11")
  print("============= Diagnosing the results ================")
  print("check.field()")
  x <- eval(parse(text=paste("FCST$",flds[length(flds)],sep="")))
  bitmap("check.field1.png",type="png256")
  plot(range(x$lon),range(x$lat),type="n",main="check.field1: ensemble members")
  grid()
  colour <- c("black","grey30","red","darkred","blue","darkblue","magenta","green","darkgreen","brown")
  colour <- rep(colour,4)
  addland()
  MAP <- matrix(x$t2m[,,1,2],73,38)*0
  for (i in 1:n.ens) {
    map <- 0.3333*(matrix(x$t2m[,,i,2],73,38)+matrix(x$t2m[,,i,3],73,38)+matrix(x$t2m[,,i,4],73,38))
    contour(x$lon,x$lat,map,add=TRUE,col=colour[i],levels=seq(1,4,by=1))
    contour(x$lon,x$lat,map,add=TRUE,col=colour[i],levels=seq(-4,-1,by=1),lty=2)
    contour(x$lon,x$lat,map,add=TRUE,col=colour[i],levels=c(0,0),lwd=2)
    MAP <- MAP + map
  }
  MAP <- MAP/n.ens
  if (label) mtext(side=1,date(),cex=0,7,col="grey")
  dev.off()

  bitmap("check.field2.png",type="png256")
  plot(range(x$lon),range(x$lat),type="n",main="check.field2: ensemble mean")
  grid()
  addland()
  contour(x$lon,x$lat,MAP,add=TRUE,levels=seq(0.5,4,by=0.5))
  contour(x$lon,x$lat,MAP,add=TRUE,levels=seq(-4,-0.5,by=0.5),lty=2)
  contour(x$lon,x$lat,MAP,add=TRUE,levels=c(0,0),lwd=2)
  if (label) mtext(side=1,date(),cex=0,7,col="grey")
  dev.off()
  results <- list(lon=x$lon,lat=x$lat,map=MAP,tim="1-month-lead",
                  date=flds[length(flds)],description="check.field (getecfcst.R)",
                  attributes=NULL)
  class(results) <- "map"
  attr(results,"long_name") <- "2 meter temperature"
  attr(results,"descr") <- "ECMWF seasonal forecast"
  invisible(results)
}


test.code <- function(stnr = 69100,normal.period=c(1961,1990)) {
  obs <- KDVH4DS(StNr=stnr)
  obs.ts <- plotStation(obs,l.anom=FALSE,what="t",trend=FALSE,std.lev=FALSE)
  obs.ts$value.3m <- seasonal.mean(obs.ts$value)

  klima <- adjust.ts(obs,test=TRUE)$klima
  y <- obs.ts$value.3m - klima; yy <- obs.ts$yy; mm <- obs.ts$mm
  print(summary(y))

  clim <- obs
  clim$val <- obs$val - anomaly.station(obs)$val

  mon.ave <- rep(colMeans(obs$val,na.rm=TRUE),dim(clim$val)[1])
  yy.2 <-  sort(rep(obs$yy,12)); mm.2 <- rep(1:12,length(obs$yy)); yymm.2 <- yy.2+mm.2
  Clim <- plotStation(clim,l.anom=FALSE,what="t",trend=FALSE,std.lev=FALSE,lwd=1,lty=1,col="darkgreen")

  i1 <- is.element(Clim$yymm,yymm.2)
  i2 <- is.element(yymm.2,Clim$yymm)
  ii1 <- is.element(obs.ts$yymm,yymm.2)
  ii2 <- is.element(yymm.2,obs.ts$yymm)
  
  print(c(length(mon.ave),length(Clim$value)))

  plot(Clim$yymm,Clim$value,type="b",pch=19,col="grey35",
       main=obs$location,ylab=obs$obs.name,xlab="time",xlim=c(1990,2010) ,ylim=1.25*range(c(y,Clim$value),na.rm=TRUE))
  lines(obs.ts$yymm,klima,col="blue",lty=2)
  points(obs.ts$yymm,klima,col="blue")
  grid()

  obs.ts <- plotStation(obs,l.anom=TRUE,what="t",trend=FALSE,std.lev=FALSE,normal.period=normal.period)
  obs.ts$value.3m <- seasonal.mean(obs.ts$value)

  plot(obs.ts$yymm,obs.ts$value,pch=19,col="grey35",main=obs$location,ylab=obs$obs.name,xlab="time",xlim=c(1996,2006),
       ylim=1.5*range(y,obs.ts$value,na.rm=TRUE))
  lines(obs.ts$yymm,obs.ts$value.3m,lwd=2,col="red")
  grid()

  points(yy+(mm-0.5)/12, y,pch=20,cex=0.7,col="red")
  points(yy+(mm-0.5)/12, y,pch=21,cex=0.9,col="darkred")
}


#-------------------------------- monthly forecasts!


Monthly <- function(location=18700,vname="T2m",path="/klimadata/rasmusb/monthly/",verbose=FALSE,SAVE=TRUE) {

  obs <- KDVH(location)
  obs$RR[obs$RR<0] <- 0
  if (upper.case(vname)=="T2M") obs$t2m <- ma.filt(obs$TAM,7) else
  if (upper.case(vname)=="TP") obs$t2m <- ma.filt(obs$RR,7)

  N <- length(obs$t2m); w <- 2*pi*seq(1,N,by=1)/365.25
#  x11()
#  plot(obs$Year+(obs$Month-1)/12+(obs$Day-1)/365.25,obs$t2m,col="grey80",pch=20,xlim=c(2004.5,2007.25),cex=0.7,
#       main=obs$Location,sub="ECMWF monthly forecast",xlab="Time",ylab=vname)
#  grid()
  normal.period <- is.element(obs$Year,seq(1961,1990,by=1))
  AC.61.90 <- data.frame(y=obs$t2m[normal.period],x1=cos(w[normal.period]),x2=sin(w[normal.period]),
                   x3=cos(2*w[normal.period]),x4=sin(2*w[normal.period]),x5=cos(3*w[normal.period]),x6=sin(3*w[normal.period]))
  # The anomalies for the monthly forecasts are with respect to 1900-2001:
  ec.period <- is.element(obs$Year,seq(1990,2001,by=1))
  AC.90.01 <- data.frame(y=obs$t2m[ec.period],x1=cos(w[ec.period]),x2=sin(w[ec.period]),
                   x3=cos(2*w[ec.period]),x4=sin(2*w[ec.period]),x5=cos(3*w[ec.period]),x6=sin(3*w[ec.period]))
  AC <- data.frame(y=obs$t2m, x1=cos(w), x2=sin(w), x3=cos(2*w), x4=sin(2*w), x5=cos(3*w), x6=sin(3*w))
  clim.mod <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data=AC.61.90)
  clim.mod.ec <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data=AC.90.01)
  climatology <- predict(clim.mod,newdata=AC)
  clim.adj <- predict(clim.mod,newdata=AC) - predict(clim.mod.ec,newdata=AC)
#  lines(obs$Year+(obs$Month-1)/12+(obs$Day-1)/365.25,climatology,lwd=2,col="red")
#  lines(obs$Year+(obs$Month-1)/12+(obs$Day-1)/365.25,climatology-clim.adj,lty=3,col="red")

  recent <- julday(obs$Month,obs$Day,obs$Year) > julday(6,1,2004)
  N <- sum(recent)
  dat <- rep(NA,N*4*51); dim(dat) <- c(N,51,4); DAT <- dat
  act <- rep(NA,N*4); dim(act) <- c(N,4); day <- act; ano <- act
  
  if (path=="/opdata/sesong/") name.id <- "smma"
  if (substr(path,1,26)=="/klimadata/rasmusb/monthly") name.id <- "MonthlyFcst"
  print(paste("Data path=",path))
  vname2 <- switch(upper.case(vname),"T2M"="2T","MSL"="MSL","TP"="TP",
                    "P2T"="2T","MSL"="MSL","SLP"="MSL","TEMP"="2T","PRE"="TP","PREC"="TP",
                    "PRECIP"="TP","167.171"="2T")
  file.list <- list.files(path=path,pattern=".nc")
  file.list <- file.list[grep("MonthlyFcst",file.list)]
  file.list <- file.list[grep(vname,file.list)]
  print(vname2)
  print(file.list)
  nf <- length(file.list)
  for (ifcst in 1:nf) {
     print(file.list[ifcst])
     dots <- instring(".",file.list[ifcst])
     fcst <- substr(file.list[ifcst],dots[2],dots[3]-1)
     ncid<-open.ncdf(paste(path,file.list[ifcst],sep=""))
     nv <- ncid$nvars
     cdfvars <- rep("-",nv)
     for (i in 1:nv) cdfvars[i] <- ncid$var[[i]]$name
     miss <-  ncid$var[[1]]$missval
     nd <- ncid$var[[1]]$ndims
     cdfdims <- rep("-",nd)
     for (i in 1:nd) cdfdims[i] <- ncid$var[[1]]$dim[[i]]$name
     #print(cdfdims)

     ilon <- grep("lon",lower.case(cdfdims))
     ilat <- grep("lat",lower.case(cdfdims))
     itim <- grep("tim",lower.case(cdfdims))
     inum <- grep("num",lower.case(cdfdims))
     ifcp <- grep("fcperiod",lower.case(cdfdims))
     lon <- get.var.ncdf(ncid,cdfdims[ilon])
     lat <- get.var.ncdf(ncid,cdfdims[ilat])

     #print(lat);print(lon)
     
     if (ifcst==1) { lon0 <- lon; lat0 <- lat }
     num <- get.var.ncdf(ncid,cdfdims[inum])
     tim <- get.var.ncdf(ncid,cdfdims[itim])
     fcperiod <- get.var.ncdf(ncid,cdfdims[ifcp])
     print(fcperiod)
     
     attr(lon,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[ilon],"$units",sep="")))
     attr(lat,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[ilat],"$units",sep="")))
     attr(tim,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[itim],"$units",sep="")))
     attr(num,"unit") <- eval(parse(text=paste("ncid$dim$",cdfdims[inum],"$units",sep="")))
     t.unit <- attr(tim,"unit")

     torg <- substr(t.unit,regexpr("since",t.unit)+6,nchar(t.unit))
     t.unit <- substr(t.unit,1,regexpr("since",t.unit)-2)
  
     if (!is.null(torg)) {
       print(paste("torg=",torg))
       yy0 <- datestr2num(torg)[1]
       mm0 <- datestr2num(torg)[2]
       dd0 <- datestr2num(torg)[3]
     } 

     if (verbose) print(paste("Time unit:",lower.case(t.unit),"  yy0=",yy0," mm0=",mm0," dd0=",dd0))
     mmddyy <- caldat(tim/24 + julday(mm0,dd0,yy0))
     mm <- mmddyy$month
     yy <- mmddyy$year
     dd <- mmddyy$day
     tim <- tim/24
     t.unit <- "day"
     #print(tim)
     print(rbind(yy,mm,dd))
           
     daysayear<- 365.25
     i1 <- sum(lon < obs$lon-2)
     i2 <- sum(lon <= obs$lon+2)-i1
     j1 <- sum(lat >= obs$lat+2)
     j2 <- sum(lat >= obs$lat-2)-j1
     lon <- lon[i1:(i1+i2-1)]; lat <- lat[j1:(j1+j2-1)]
     nx <- length(lon); ny <- length(lat); nt <- length(tim)
     lonxy <- rep(lon,length(lat)); latxy <- sort(rep(lat,length(lon)))
     #print(c(nx,ny,nt))
     #print(c(obs$lon,NA,lon))
     #print(c(obs$lat,NA,lat))
     
     start <- c(i1,j1,1,1,1)
     count <- c(i2,j2,51,4,nt)
     #print("start:")
     #print(start)
     #print("count:")
     #print(count)
     v1 <- cdfvars[1]
     #data <- get.var.ncdf(ncid,v1)
     #print("dim(data):"); print(dim(data))
     data <- get.var.ncdf(ncid,v1,start=start,count=count)
     #print("51 member values:")
     #print(data[1,1,,1,1])
     #print("4 forecast preiods:")
     #print(data[1,1,1,,1])
     print(paste("length(data)=",length(data)))
     arv <- att.get.ncdf(ncid, cdfvars[1], 'scale_factor')
     if( arv$hasatt ) scal <- arv$value else scal <- 1
     arv <- att.get.ncdf(ncid, cdfvars[1], 'add_offset')
     if( arv$hasatt ) offs <- arv$value else offs <- 0
     if (as.numeric(R.Version()$major) < 2) { 
       print(paste("The host is",Sys.info()[4],"running R version",R.Version()$major,"-",
                    R.Version()$minor,": scaling by x",scal,"and adding",offs))
       data <- data * scal
       data <- data + offs
     }
     close.ncdf(ncid)

     print("read the data...")
     dims <- dim(data)
     if (length(dims)==4) {dim(data) <- c(dims[1],dims[2],51,4,1); print(dim(data))}

     # sort the data:
     data <- data[,,,order(fcperiod),]
     if (length(dims)==4) {dim(data) <- c(dims[1],dims[2],51,4,1)}
     
     Clim.adj <- as.numeric(clim.adj[recent])
     Climatology <- climatology[recent]
     Year <- obs$Year[recent]; Month <- obs$Month[recent]; Day <- obs$Day[recent]
     jdays <- julday(Month,Day,Year); jday0 <- julday(1,1,2004)

     #print(dim(data)); print(dim(dat)); print(dim(act))

     if (lower.case(vname)=="tp")  {conv.scal <- 3600*24 * 1000; add.clim <- 0}else
                                   {conv.scal <- 1; ; add.clim <- 1}

     for (ii in 1:length(yy)) {
       for (ix in 1:4) {
         #print(c(mm[ii],dd[ii],yy[ii]))
#         idate <- is.element(jdays,julday(mm[ii],dd[ii],yy[ii])+ix*7+1)
         idate <- is.element(jdays,julday(mm[ii],dd[ii],yy[ii]))
         jdate <- is.element(julday(obs$Month,obs$Day,obs$Year),julday(mm[ii],dd[ii],yy[ii])+ix*7+1)
         #print(c(sum(idate),sum(jdate)))
         if ((sum(idate)==1) & (sum(jdate)==1)) {
           act[idate,ix]    <- obs$t2m[jdate]
           if (lower.case(vname)!="tp") ano[idate,ix]    <- obs$t2m[jdate] - Climatology[idate] else
                                        ano[idate,ix]    <- obs$t2m[jdate]
         }
         day[idate,ix] <- julday(mm[ii],dd[ii],yy[ii]) + 7*ix +1 - jday0
         #print(paste("[idate]=",(1:length(idate))[idate],"  julday(mm[ii],dd[ii],yy[ii])=",julday(mm[ii],dd[ii],yy[ii])))
         if (sum(idate)==1) {
           for (iv in 1:51) {
             map <- data[,,iv,ix,ii]
             if (lower.case(vname)!="tp") {
               dat[idate,iv,ix] <- interp(lonxy,latxy,map,obs$lon,obs$lat)$z*conv.scal +
                                          add.clim*(Climatology[idate]-Clim.adj[idate]) 
               DAT[idate,iv,ix] <- interp(lonxy,latxy,map,obs$lon,obs$lat)$z*conv.scal - add.clim*Clim.adj[idate]
             } else {
               # For precip, the full field is given.
               dat[idate,iv,ix] <- interp(lonxy,latxy,map,obs$lon,obs$lat)$z*conv.scal
               DAT[idate,iv,ix] <- interp(lonxy,latxy,map,obs$lon,obs$lat)$z*conv.scal 
               
             }
#             points(yy[ii]+(mm[ii]-1)/12+(dd[ii]-1)/365.25 + ix*7/365.25,
#                    dat[idate,iv,ix],cex=0.4,col="steelblue")
           }
         }
       }
     }
     rm(data)
   }  # End loop

   # Convert from units = "m s**-1" to units = "mm/day"
   print("Got all the data!")
   #print(dim(dat)); print(length(Year))
#   points(Year+(Month-1)/12+(Day+8)/365.25, act[,1],col="grey30",pch=21,xlim=c(2004.5,2007.25),cex=0.9)
#   points(Year+(Month-1)/12+(Day+15)/365.25,act[,2],col="grey30",pch=21,xlim=c(2004.5,2007.25),cex=0.9)
#   points(Year+(Month-1)/12+(Day+22)/365.25,act[,3],col="grey30",pch=21,xlim=c(2004.5,2007.25),cex=0.9)
#   points(Year+(Month-1)/12+(Day+29)/365.25,act[,4],col="grey30",pch=21,xlim=c(2004.5,2007.25),cex=0.9)
  
#   dev.copy2eps(file=paste("monthly_",strip(obs$Location),".eps"))

   results=list(forecast=dat,julianday=day,jday0=jday0,clim.adj=clim.adj[recent],clim.61.90=Climatology,
                observed=act,time.origin="01.01.2004",Year=Year,Month=Month,Day=Day,
                location=obs$Location,vname=vname,fcperiod=sort(fcperiod),fcst.anom=DAT,obs.anom=ano)
   class(results) <- "monthly.forecasts"
   if (SAVE) save(file=paste("monthly_",strip(obs$Location),".rda"),results)
   invisible(results)

}  # end function


Monthly2ensemble.timeseries <- function(x,lag=1) {
   mid.pt <- round(mean(diff(x$fcperiod))/2)
   mmddyy <- caldat(x$julianday + x$jday0 + x$fcperiod[lag] + mid.pt)
   mm <- mmddyy$month
   yy <- mmddyy$year
   dd <- mmddyy$day
   a <- list(val=x$observed[,lag],yy=yy,mm=mm,dd=dd)
   ensemble.timeseries <- list(ts=x$forecast[,,lag],a=a,
                               months=NA,mon3=NA)
   invisible(ensemble.timeseries)
}

test.reliab <- function(test=1,S2N=1,y0=0,CDF=TRUE,N=1000,M=51,f=1/365.25,ens.spread=0.33,type="png256") {
  # tricking reliab.diag to think it deals with real forecasts:
  FCST <- list(vname="test")
  class(FCST) <- "monthly.forecasts"
  w <- 2*pi*seq(1,N,by=1)*f
  ts <- rep(NA,N*M); dim(ts) <- c(N,M)
  if (test==1) {
    obs <- 0.2*cos(w) + 0.3*sin(2*w) - 0.8*cos(6*w) + 0.03*rnorm(N)
    if (S2N>0) for (i in 1:N) ts[i,] <- rnorm(M,mean=obs[i],sd=1/S2N) else
               for (i in 1:N) ts[i,] <- rnorm(M,mean=rnorm(1)*2,sd=ens.spread)
  } else {
     p <- rep(NA,N); obs <- rep(y0+1,N)
     for (i in 1:N) {
       ts[i,] <- rnorm(M,mean=0.33*cos(4*w[i])*2,sd=ens.spread)
       p[i] <- sum(ts[i,]<y0)/M
     }
     for (i in seq(0.05,1.05,by=0.05)) {
       iii <- (p<i) & (p>=(i-0.05))
       srt <- order(rnorm(sum(iii)))
       ni <- round(i*sum(iii))
       iv <- (1:N)[iii][srt][1:ni]
       obs[iv] <- y0-1
       print(paste("test2: i=",i," ni=",ni," sum(iii)=",sum(iii)))
     }
     print(table(round(p,1),obs))
     
  }
  if (test==1) eval.results <- reliab.diag(obs=obs,ts=ts,y0=y0,res=150,type=type,loc=paste("S2N",S2N,sep=""),
                                           a=NULL,yy=NULL,mm=NULL,dd=NULL,FCST=FCST,CDF=CDF) else
                               
                               reliab.diag(obs=obs,ts=ts,y0=y0,res=150,type=type,loc="test2",
                                           a=NULL,yy=NULL,mm=NULL,dd=NULL,FCST=FCST,CDF=CDF)
#  x11()
#  plot(c(1,N),range(c(obs,ts)),type="n")
#  grid()
#  for (i in 1:N) { points(rep(i,M),ts[i,],col="grey"); points(i,obs[i],pch=19,col="red") }

#  x11()
#  plot(p,obs)
}

test.hitrate <- function(test=1,S2N=1,y0=0,CDF=TRUE,N=1000,M=51,f=1/365.25,type="png256") {
  # tricking reliab.diag to think it deals with real forecasts:
  FCST <- list(vname="test")
  class(FCST) <- "monthly.forecasts"
  w <- 2*pi*seq(1,N,by=1)*f
  obs <- 0.2*cos(w) + 0.3*sin(2*w) - 0.8*cos(6*w) + 0.03*rnorm(N)
  ts <- rep(NA,N*M); dim(ts) <- c(N,M)
  if (S2N>0) for (i in 1:N) ts[i,] <- rnorm(M,mean=obs[i],sd=1/S2N) else
             for (i in 1:N) ts[i,] <- rnorm(M,mean=rnorm(1)*2)*0.33
  x <- list(obs=obs,ens=ts,lag=0,location=paste("S2N",round(S2N),sep=""),num=N,vname="test")
  x$obs.name="test"

  print("Hit-rate:")
  hit.results <- ec.hitrate(x,plot=TRUE,potent=FALSE,label=TRUE,type=type)

  x11()
  plot(c(1,N),range(c(obs,ts)),type="n")
  grid()
  for (i in 1:N) { points(rep(i,M),ts[i,],col="grey"); points(i,obs[i],pch=19,col="red")  }
  
}

eval.Monthly1 <- function(x,lag=1,y0=0,CDF=TRUE,type="png256") {
  mid.pt <- round(mean(diff(x$fcperiod))/2)
  good <- is.finite(x$julianday[,lag]) 
  #print(summary(x$julianday[good,lag] + x$jday0 + x$fcperiod[lag] + mid.pt))
  mmddyy <- caldat(x$julianday[good,lag] + x$jday0 + x$fcperiod[lag] + mid.pt)
  mm <- mmddyy$month
  yy <- mmddyy$year
  dd <- mmddyy$day
  a <- list(val=x$observed[good,lag],yy=yy,mm=mm,dd=dd)
  eval.results <- reliab.diag(obs=x$obs.ano[good,lag],ts=x$fcst.anom[good,,lag],y0=y0,res=150,type=type,loc=strip(x$location),
                              a=a,yy=x$Year[good],mm=x$Month[good],dd=x$Day[good],FCST=x,CDF=CDF,lag=lag)
  eval.results
}

eval.Monthly2 <- function(x,lag=1,type="png256") {
  nt <- length(x$Year)
  dates <- x$Year+(x$Month-1)/12+(x$Day-1)/365.25
  q05 <- rep(NA,nt); q95 <- q05; q90 <- q05; q10 <- q05; q20 <- q05; q30 <- q05; q40 <- q05
  q50 <- q05; q60 <- q05; q70 <- q05; q80 <- q05; mx <- q05;  mn <- q05;  m <- q05 
  for (it in 1:nt) {
    q05[it] <- quantile(x$fcst.anom[it,,lag],0.05,na.rm=TRUE)
    q95[it] <- quantile(x$fcst.anom[it,,lag],0.95,na.rm=TRUE)
    q10[it] <- quantile(x$fcst.anom[it,,lag],0.10,na.rm=TRUE)
    q90[it] <- quantile(x$fcst.anom[it,,lag],0.90,na.rm=TRUE)
    q20[it] <- quantile(x$fcst.anom[it,,lag],0.20,na.rm=TRUE)
    q80[it] <- quantile(x$fcst.anom[it,,lag],0.80,na.rm=TRUE)
    q30[it] <- quantile(x$fcst.anom[it,,lag],0.30,na.rm=TRUE)
    q70[it] <- quantile(x$fcst.anom[it,,lag],0.70,na.rm=TRUE)
    q40[it] <- quantile(x$fcst.anom[it,,lag],0.40,na.rm=TRUE)
    q60[it] <- quantile(x$fcst.anom[it,,lag],0.60,na.rm=TRUE)
    q50[it] <- quantile(x$fcst.anom[it,,lag],0.50,na.rm=TRUE)
    mx[it] <- max(x$fcst.anom[it,,lag],na.rm=TRUE)
    mn[it] <- min(x$fcst.anom[it,,lag],na.rm=TRUE)
    m[it] <- mean(x$fcst.anom[it,,lag],na.rm=TRUE)
  }

  if (upper.case(x$vname)=="TP") x$obs.anom[x$obs.anom<0] <-  0
  good <- is.finite(m) & is.finite(x$obs.anom[,lag])
  r <- round( cor(m[good],x$obs.anom[,lag][good]),2 )
  
  bitmap(paste("eval2_",strip(x$location),x$vname,"_",lag,".",substr(type,1,3),sep=""),type,res=300) 
  plot(range(dates[good]),range(c(x$obs.anom[,lag],x$fcst.anom[,,lag]),na.rm=TRUE),
       type="n",xlab="Time",ylab=x$vname,
       main=paste("Week anomaly fcst:",x$location,x$vname," r=",r))
  grid()
  
  polygon( c(dates[good],reverse(dates[good])),c(mx[good],reverse(mn[good])), border="yellow",lwd=1 )
  polygon( c(dates[good],reverse(dates[good])),c(q95[good],reverse(q05[good])), col="grey90",border="grey90",lwd=1 )
  polygon( c(dates[good],reverse(dates[good])),c(q90[good],reverse(q10[good])), col="grey80",border="grey80",lwd=1 )
  polygon( c(dates[good],reverse(dates[good])),c(q80[good],reverse(q20[good])), col="grey70",border="grey70",lwd=1 )
  polygon( c(dates[good],reverse(dates[good])),c(q70[good],reverse(q30[good])), col="grey60",border="grey60",lwd=1 )
  polygon( c(dates[good],reverse(dates[good])),c(q60[good],reverse(q40[good])), col="grey50",border="grey50",lwd=1 )
  lines(dates[good],q50[good],lwd=3,col="grey40")

  points(dates[good],x$obs.anom[,lag][good],pch=19)
  
  dev.off()

  r
}


eval.Monthly3 <- function(x,lag=1,y0=0,CDF=TRUE,type="png256") {
  x <- list(obs=x$obs.anom[,lag],ens=x$fcst.anom[,,lag],lag=lag,location=x$location,
            num=length(x$fcst.anom[1,,lag]),vname=x$vname)
  x$obs.name="test"

  print("Hit-rate:")
  hit.results <- ec.hitrate(x,plot=TRUE,potent=FALSE,label=TRUE,type=type)
  print(hit.results)
  hit.results
}
 

Monthly4EBL <- function(type="pdfwrite") {
  locations <- c(18700,07010,39040,44560,50540,69100,82290,90450)
  nl <- length(locations)

  table1 <- rep(NA,nl*4*2); dim(table1) <- c(nl,4,2)
  table2 <- table1; table3 <- table1

  b <- rep(NA,4); r <- b; p <- b
  
  i <- 0; ii <- 0
  for (location in locations) {
    i <- i+1
    for (vname in c("T2M","TP")) {
      ii <- ii+1
      x <- Monthly(location=location,vname=vname)
      for (lag in 1:4) b[lag] <- round(eval.Monthly1(x,lag=lag,type=type)$bs,2)
      for (lag in 1:4) r[lag] <- eval.Monthly2(x,lag=lag,type=type)
      for (lag in 1:4) p[lag] <- round(eval.Monthly3(x,lag=lag,type=type)$p.hit,2)
      table2[i,,ii] <-  c(b)
      table1[i,,ii] <-  c(r)
      table3[i,,ii] <-  c(p)
    }

    ii <- 0
  }

  while (dev.cur()>1) dev.off()
  table1a <- table1[,,1]
  rownames(table1a) <- locations
  colnames(table1a) <- c("Day 5-11","Day 12-18","Day 19-25","Day26-32")
  table1b <- table1[,,2]
  rownames(table1b) <- locations
  colnames(table1b) <- c("Day 5-11","Day 12-18","Day 19-25","Day26-32")
  table2a <- table2[,,1]
  rownames(table2a) <- locations
  colnames(table2a) <- c("Day 5-11","Day 12-18","Day 19-25","Day26-32")
  table2b <- table2[,,2]
  rownames(table2b) <- locations
  colnames(table2b) <- c("Day 5-11","Day 12-18","Day 19-25","Day26-32")
  table3a <- table3[,,1]
  rownames(table3a) <- locations
  colnames(table3a) <- c("Day 5-11","Day 12-18","Day 19-25","Day26-32")
  table3b <- table3[,,2]
  rownames(table3b) <- locations
  colnames(table3b) <- c("Day 5-11","Day 12-18","Day 19-25","Day26-32")

  
  print("T2m - correlations:")
  print.table(table1a)
  print("RR - correlations:")
  print.table(table1b)

  print("T2m - Brier scores:")
  print.table(table2a)
  print("RR - Brier Scores:")
  print.table(table2b)

  print("T2m - Hit-ratio:")
  print.table(table3a)
  print("RR - Hit-ratio:")
  print.table(table3b)
  
  test.reliab(type=type)
  test.reliab(S2N=0,type=type)
  test.reliab(S2N=10,type=type)
  test.reliab(test=2,type=type)
  test.hitrate(type=type)
  test.hitrate(S2N=0,type=type)
  test.hitrate(S2N=10,type=type)
  
}

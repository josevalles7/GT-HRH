#--Subrutinas para cómputo de pronósticos por Analogía (pronóstico subestacional a estacional). Versión Beta--#
# Todavía restan corregir algunos errores de código que producen NAs. Por otro lado la aproximación seleccionada fue seleccionar analogías en dominio de Índice Estandarizado de Caudal (ajuste a gamma), con ventanas móviles (por defecto 30 días 6 hacia atrás y 3 hacia adelante para pronos). La métrica de similitud seleccionada es RMSE entre las ventanas móviles precedentes que cubren los mismos días del año que la ventana 'actual' evaluada. Los pesos se construyeros computando (1-R^2)^(-2) de cada ventana móvil análoga (IDW, con k=2), con R^2 siendo el coeficiente de determinación entre el SSI de cada ventana al respecto de la ventana 'actual'. Preferí esto por que el R^2 mide la asociación en términos de fase, mientras un bajo RMSE puede corresponder a una señal próxima pero con tendencia muy distinta (e.g. oscilante si se superpone sobre la actual). Hay que seguir tarbajando, esto lo escribí a la vyuelta de Bogotá y después, vacaciones y comenzó el año con todo. LMG 20230310
#Librerías Requeridas. 
library(fitdistrplus)
library(lubridate)
library(xts)
source("~/01-PROCEDIMIENTOS/02-FUNCIONES/LinearSystems.R")
#función de testeo
TestValue<-function(val1,val2){
  if(val1!=val2){
    return(as.numeric(1))
  }
  else{
    return(as.numeric(0))
  }
}
#Computa la raíz cuadrada del error medio cuadrático (métrica de distancia)
RMSE<-function(sim,obs){
  return(sqrt(mean((obs-sim)^2)))
}
#Computa la longitud de registros para una serie xts y devuelve un vector con TimeStart (date),TimeEnd (date) y duración m años (integer) 
ComputeRecordsLenght<-function(serie,Window=30,N=6,Horizon=3){
  TimeStart=min(index(serie))
  TimeEnd=max(index(serie))
  m=as.numeric(year(TimeEnd))-as.numeric(year(TimeStart))+1
  type=0 #Las observaciones p/ajuste y los pronósticos a elaborar corresponden al mismo año civil
  if((year(TimeEnd+Horizon*Window)>year(TimeEnd))|(year(TimeEnd-N*Window)<year(TimeEnd))){
    type=1 #El pronóstico a elaborar contendrá elementos del próximo año civil o algunas de las observaciones utilizadas corresponden al año civil precedente
  }
  pars=c(TimeStart,TimeEnd,as.numeric(m),as.numeric(type))
  return(pars)
}
#Agrega serie diaria xts a dt=Window [StartPeriod/EndPeriod]. Entrada: serie temporal regularizada a paso diario (objeto xts)
AggDailySerieByWindow<-function(DailySerie,Window=30){
  return(na.approx(rollapply(DailySerie,Window,mean,align = 'right',by.column=T)))
}
#Devuelve valor agregado en ventana de ancho start/end en objeto xts
GetWindowAvgVal<-function(serie,start,end,fun='mean'){
  sub=subset(serie,index(serie)>=start&index(serie)<=end)
  if(fun=='mean'){
    return(mean(as.numeric(sub),na.rm=T))
  } 
  if(fun=='max'){
    return(max(as.numeric(sub),na.rm=T))
  }
  if(fun=='min'){
    return(min(as.numeric(sub),na.rm=T))
  }
}
#Declara xts vacía
CreateForecastXTS<-function(start=Sys.Date(),Window=30,Horizon=3){
  dt=paste0(Window," days")
  t=seq(start-Window*Horizon,start+Window*Horizon,by=dt)
  serie=xts(order.by=t)
  return(serie)
}
#Obtiene estadísticos sobre ventana de agregación temporal con dt=Window y LeadTime en unidades dt. 
GetWindowStats<-function(serie,StartPeriod=as.Date("1991-01-01"),EndPeriod=as.Date("2020-12-31"),LeadTime=1,Window=30){
  Start=max(index(serie))
  Startdoy=yday(Start+Window*(LeadTime-1))
  Enddoy=yday(Start+Window*LeadTime)
  doy=yday(index(serie))
  year=year(index(serie))
  StartRef=year(StartPeriod)
  EndRef=year(EndPeriod)
  stats=c()
  for(y in seq(StartRef,EndRef)){
    if(year(Start+Window*LeadTime)>year(Start)){
      sub=subset(na.omit(serie),doy>=Startdoy&year==y,na.rm=T)
      sub=rbind(na.omit(sub),subset(serie,doy<=Enddoy&year==y+1,na.rm=T))
    }
    else{
      sub=subset(na.omit(serie),doy>=Startdoy&doy<=Enddoy&year==y)
    }
    stats$mean[y-StartRef+1]=mean(as.numeric(na.omit(sub)),na.rm=T)
    stats$max[y-StartRef+1]=max(as.numeric(na.omit(sub)),na.rm=T)
    stats$min[y-StartRef+1]=min(as.numeric(na.omit(sub)),na.rm=T)
  }
  return(stats)
}
#Computa Scores Estandarizados para los últimos N pasos de duración dt=Window, tomando como período de referencia [StartPeriod/EndPeriod], a partir de un objeto xts rol (salida AggDailySerieByWindow)
ComputeScores<-function(rol,StartPeriod=as.Date('1991-01-01'),EndPeriod=as.Date('2020-12-31'),Window=30,N=6){
  pars=ComputeRecordsLenght(rol,Window,N)
  TimeStart=pars[1]+Window+1
  TimeEnd=pars[2]
  m=pars[3]
  years=year(index(to.yearly(rol)))
  anomaly=matrix(nrow=m,ncol=N)
  rownames(anomaly)=years
  colname=c()
  for(i in seq(0,N-1)){
      if(yday(TimeEnd-Window*i)!=366){
        doy=yday(TimeEnd-Window*i)
      }
      else{
        doy=1
      }
      colname[i+1]=doy
      sub=subset(rol,yday(index(rol))==doy&year(index(rol))>=year(StartPeriod)&year(index(rol))<=year(EndPeriod)&!is.na(rol))
      fit=fitdist(as.numeric(sub),distr="gamma",method="mme")
      sub=subset(rol,yday(index(rol))==doy&index(rol)>TimeStart)
      years=year(index(to.yearly(sub)))
      anomaly[paste(years),i+1]=qnorm(pgamma(as.numeric(na.omit(sub)),shape = coefficients(fit)[1],rate = coefficients(fit)[2])) 
  }
  colnames(anomaly)=colname
  checkyears=0
  j=0
  for(i in seq(1,N)){
    if(year(TimeEnd-Window*i)!=year(TimeEnd)){
      checkyears=1
      if(j==0){
        j=i
      }
    }
  } 
  if(checkyears==1){#deben considerarse datos del año precedente (desplazamiento rígido: traslado a fila posterior, en estas columnas. Se elaboraron estas líneas para dar solución a las ejecuciones con tamaño de ventana=30 días. Debe revisarse y generalizar para casos ejecutables (para ser consistente no ha de ser completo!)
        lagMatrix=matrix(ncol=dim(anomaly)[2],nrow=dim(anomaly)[1])
        lagMatrix[2:dim(anomaly)[1],j:N]=anomaly[1:dim(anomaly)[1]-1,j:N]
        if(i>1){
          lagMatrix[1:dim(anomaly)[1],1:j]=anomaly[1:dim(anomaly)[1],1:j]
        }
        rownames(lagMatrix)=rownames(anomaly)
        colnames(lagMatrix)=colnames(anomaly)
        anomaly=lagMatrix
  }
  return(anomaly)
}
#Computa pesos de series análogas. Error es un vector con el error obtenido por cada una de las series análogas (salida de ComputeAnalogyAnalysis<)
ComputeWeights<-function(error,l=2){
  return(error^(-l)/sum(error^(-l)))
}
#Realiza Análisis en Búsqueda de años análogos a partir de una serie diaria xts agregada a dt=Window, tomando un período de referencia, considerando los últimos N intervalos de logitud Window y seleccionando los M años hidrológicos más semejantes como predictores
ComputeAnalogyAnalysis<-function(DailySerie,Window=30,StartPeriod=as.Date('1991-01-01'),EndPeriod=as.Date('2020-12-31'),N=6,M=5,metric='rmse'){
  analysis=list()
  test=0
  rol=AggDailySerieByWindow(DailySerie,Window)
  year=year(index(to.yearly(rol)))
  anomaly=ComputeScores(rol,StartPeriod,EndPeriod,Window,N)
  analysis$pars=ComputeRecordsLenght(rol,Window,N)
  error=c()
  r2=c()
  nash=c()
  bias=c()
  offset=c()
  m=analysis$pars[3]
  for(i in seq(1,m-1)){
    error[i]=RMSE(anomaly[m,1:N],anomaly[i,1:N])
    if(all(is.na(anomaly[i,1:N]))){# se asigna valor 0 para que el peso sea el mínimo y la estimación, en todo caso, sea el valor medio para el mes, si todos los predictores son NA
      bias[i]=0
      offset[i]=0
      r2[i]=0
    }
    else{
      model=summary(lm(anomaly[m,1:N]~anomaly[i,1:N]))
      bias[i]=model$coef[2]
      offset[i]=model$coef[1]  
      r2[i]=model$r.squared
    }
    nash[i]=1-error[i]^2*N/var(anomaly[m,1:N])
  }
  analysis$anomaly=anomaly
  analysis$error=error
  analysis$year=year
  analysis$test=TestValue(length(analysis$year),dim(analysis$anomaly)[1]) #realiza test sobre dimensiones de matriz scores y vector de años (la cantidad de filas debe ser la misma) 
  analysis$best=order(error)[1:M]
  analysis$r2=r2[analysis$best]
  analysis$nash=nash[analysis$best]
  analysis$bias=bias[analysis$best]
  analysis$offset=offset[analysis$best]
  analysis$members=analysis$year[analysis$best]
  if(metric=='rmse'){
    print('Computing weights using rmse')
    analysis$weights=ComputeWeights(analysis$error[analysis$best])
  }
  if(metric=='r2'){
    print('Computing weights using 1-r2')
    analysis$weights=ComputeWeights((1-analysis$r2))
  }
  return(analysis)
}
#Computa la matriz de scores de años semejantes y realiza el pronóstico de score como media ponderada
ComputeForecastScores<-function(rol,Horizon=3,Window=30,N=6,StartPeriod=as.Date('1991-01-01'),EndPeriod=as.Date('2020-12-31')){
  pars=ComputeRecordsLenght(rol,Window,N)
  TimeStart=pars[1]+Window+1
  TimeEnd=pars[2]
  colname=c()
  rowname=c()
  m=as.numeric(pars[3])
  years=year(index(to.yearly(rol)))
  anomaly=matrix(nrow=m,ncol=Horizon)
  rownames(anomaly)=years
  for(i in seq(1,Horizon)){
    colname[i]=yday(TimeEnd+Window*i)
    sub=subset(rol,yday(index(rol))==yday(TimeEnd+Window*i)&year(index(rol))>=year(StartPeriod)&year(index(rol))<=year(EndPeriod)&!is.na(rol))
    fit=fitdist(as.numeric(sub),distr="gamma",method="mme")
    sub=subset(rol,yday(index(rol))==yday(TimeEnd+Window*i)) 
    years=year(index(to.yearly(sub)))
    anomaly[paste(years),i]=qnorm(pgamma(as.numeric(na.approx(sub)),shape = coefficients(fit)[1],rate = coefficients(fit)[2]))
  }
  #anomaly[paste(year(TimeEnd)),1:Horizon]=NA
  colnames(anomaly)=colname
  return(na.approx(anomaly))
}
#realiza análisis de analogía y devuelve una matriz con la serie rol de los años más semejante y el considerado para el pronóstico
ComputePersistenceForecast<-function(DailySerie,Window=30,StartPeriod=as.Date('1991-01-01'),EndPeriod=as.Date('2020-12-31'),N=6,M=5,Horizon=3,adjust='y',metric='r2'){
  print("Analog Records Analysis...")
  analysis=ComputeAnalogyAnalysis(DailySerie,Window,StartPeriod,EndPeriod,N,M,metric)
  print("done")
  best=analysis$best
  w=analysis$weights
  a0=analysis$offset
  a1=analysis$bias
  print("Computing rolling mean...")
  rol=AggDailySerieByWindow(DailySerie,Window)
  print("done")
  print("Computing Forecasts Scores...")
  forward=ComputeForecastScores(rol,Horizon,Window,N,StartPeriod,EndPeriod)
  print("done")
  start=max(index(rol))
  members=matrix(nrow=M+1,ncol=Horizon,data=0)
  forecast=matrix(nrow=1,ncol=Horizon,data=0)
  for(i in seq(1,M)){
    for(step in seq(1,Horizon)){
      if(year(start+Window*step)>year(start)){
        index=best[i]+1
      }
      else{
        index=best[i]
      }
      members[i,step]=forward[index,step]
    }
  }
  for(i in seq(1,Horizon)){
    if(adjust=='y'){
      forecast[1,i]=sum((a1*members[1:M,i]+a0)*w)
    }
    else{
      forecast[1,i]=sum(members[1:M,i]*w)
    }
  }
  members[M+1,1:Horizon]=forecast
  foreCastTitle=paste0('forecast',"_",start,"/",start+Horizon*Window,"_W_",Window)
  if(year(start+Window*Horizon)>year(start)){
    rownames(members)=c(paste0(analysis$year[best],"/",analysis$year[best]+1,sep=""),foreCastTitle)
  }
  else{
    rownames(members)=c(analysis$year[best],foreCastTitle)
  }
  return(members)
}
#Realiza análisis de persistencia y devuelve pronósticos y observaciones pasadas, en xts y matriz de probabilidades
GetForecastTimeSeries<-function(serie,Window=30,Horizon=3,N=6,M=5,StartPeriod=max(index(serie))-365*30,EndPeriod=max(index(serie))-365,metric='r2',K=1.68,bests=M,useFit='y'){
  Start=max(index(serie))
  scores=ComputePersistenceForecast(serie,Window,StartPeriod,EndPeriod,N,M,Horizon,adjust='y',metric) 
  probs=pnorm(scores)
  step=paste0(Window," days")
  t=seq(Start-Window*(N-1),Start+Window*Horizon,by=step)
  obs=c()
  forecasts=xts(order.by=t[c(N+1):c(N+Horizon)])
  fmatrix=matrix(nrow=dim(probs)[1],ncol=dim(probs)[2],data=0)
  for(i in seq(1,N-1)){
    obs[i]=GetWindowAvgVal(serie,start=t[i],end=t[i+1])
  }
  for(LeadTime in seq(1,Horizon)){
    stats=GetWindowStats(serie,StartPeriod,EndPeriod,LeadTime,Window)
    if(useFit=='y'){
      print("Using gamma fit on scores for quantile estimation")
      fit=fitdist(as.numeric(subset(stats$mean,stats$mean>=0)),distr="gamma",method="mme")
      fmatrix[,LeadTime]=qgamma(probs[,LeadTime],shape = coefficients(fit)[1],rate = coefficients(fit)[2])
    }
    else{
      print("Using empirical distribution on scores for quantile estimation")
      fmatrix[,LeadTime]=quantile(stats$mean,probs=c(probs[,LeadTime]),na.rm=T)
    }
  }
  names=c()
  for(m in seq(1:c(M+1))){
    if(m<=M){
      names[m]=paste0("member",m)
      
    }
    else{
      names[m]=paste0("CentralTrend")
    }
    forecasts=cbind(forecasts,fmatrix[m,])
  }
  names[m+1]='max'
  names[m+2]='min'
  f=list()
  f$probs=probs
  f$obs=xts(obs,order.by=t[2:N])
  max=c()
  min=c()
  for(i in seq(1:Horizon)){
    max[i]=forecasts[i,M+1]+K*sd(forecasts[i,1:bests])
    min[i]=forecasts[i,M+1]-K*sd(forecasts[i,1:bests])
  }
  forecasts=cbind(forecasts,as.numeric(max))
  forecasts=cbind(forecasts,as.numeric(min))
  colnames(forecasts)=names
  f$forecasts=as.xts(forecasts)
  return(f)
}
#--Funciones para Evaluación
#Realiza Plot de Pronóstico (Tendencia Central) y límites de confianza. Entrada, objeto f obtenido mediante GetForecastTimeSeries
PlotForecast<-function(f,include_members='y',N=length(f$obs)+1,M=dim(f$probs)[1]-1){
  lims=c(min(f$forecast,f$obs),max(f$forecast,f$obs))
  plot(rbind(f$obs,f$forecasts$CentralTrend),ylim=lims,main='Analogy Forecasts',ylab='[m³/s]')
  lines(rbind(f$obs[N-1],f$forecasts$max),col='red')
  lines(rbind(f$obs[N-1],f$forecasts$min),col='red')
  if(include_members=='y'){
    members=xts(f$obs[N-1])
    for(i in seq(1,M-1)){
      members=cbind(f$obs[N-1],members)
    }
    members=rbind(members,f$forecast[,1:M])
    lines(members,col='blue',lty=2)
  }
}
#Obtiene valores observados agregados mediante serie xts
GetObs<-function(serie,f,Window=30,Horizon=3){
  min=min(index(f$forecast))-Window
  max=max(index(f$forecast))
  v=c()
  for(i in seq(1,Horizon)){
    v[i]=GetWindowAvgVal(serie,start=min+Window*(i-1),end=min+Window*i)
  }
  step=paste0(Window," days")
  v=xts(v,order.by=seq(min+Window,max,by=step))
  colnames(v)=c('obs')
  return(v)
}
#Obtiene valor diferencial para el cómputo de la simetría direccional (consistencia en tendencia)
GetSpedVal<-function(obs0,obs1,sim1){
  if(is.na(obs0)|is.na(obs1)|is.na(sim1)){
    return(NA)
  }
  else
  {
    sObs=(obs1-obs0)
    sSim=(sim1-obs0)
    if(sObs*sSim>0){
      return(1)
    }
    else{
      return(0)
    }
  }
}
#Obtiene el valor diferencial para el cómputo de MAPE
GetdMAPE<-function(obs,sim){
  if(is.na(obs)){
    return(NA)
  }
  else{
    return(abs(sim/obs-1)*100)
  }
}
#Obtiene el valor diferencial para la estimación de probabldiad de intervalo de confianza
GetdInZone<-function(obs,simMax,simMin){
  if(is.na(obs)|is.na(simMax)|is.na(simMin)){
    return(NA)
  }
  else{
    if((obs>=simMin)&(obs<=simMax)){
    return(1)
    }
    else{
      return(0)
    }
  }
}
#Realiza subset de registros disponibles para evaluación de pronóstico (modo hindcast)
GetIntervalRecords<-function(records,forecastDate){
  interval=paste0(min(index(records)),"/",forecastDate)
  return(records[interval])
}
# png("Evals/PRTest_868_20002022.png",width=16,height=10,res=300,units='in')
# plot(forecasts,ylim=c(lims),main='Caudal Medio Móvil en 30 días y Previsión por Analogía a 30 días c/ 15 días  \nEstación Corrientes. 2000 - 2022',ylab='[m³/s]',lty=2,col='red',cex=.75)
# points(forecasts,col='blue',lty=2,pch=16,cex=.75)
# lines(obs,col='black')
# addLegend('topleft',on=1,legend.names=c('Observado','Previsión Tendencia Central'),col=c('black','blue'),lty=c(1,NA),pch=c(NA,16))

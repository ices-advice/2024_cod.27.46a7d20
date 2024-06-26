library(surveyIndex)

load("data/Q34preprocessed.RData")
load("boot/data/commonGrid.RData")

xtabs(NoAtALK~Year+Age,data=dQ34[[1]])
Q34ages=0:7


########################
## ALK
########################

dQ34 = fixAgeGroup(dQ34, 7, n=1, fun="max") # 2004 No 7+
dQ34 = addSpectrum(dQ34,by=1)
dQ34.ysplit = split(dQ34,dQ34$Year)

## Declare settings for ALK model
mf = NULL
ack=TRUE;
useBICs=TRUE;
varCofs=FALSE;
maxKs=50;

dQ34.ALK= mclapply(dQ34.ysplit,fitALK,minAge=min(Q34ages),maxAge=max(Q34ages),autoChooseK=ack,useBIC=useBICs,varCof=varCofs,maxK=maxKs,mc.cores=10)

dQ34.Nage=mclapply(dQ34.ALK,predict,mc.cores=10)
for(i in 1:length(dQ34.ALK)) dQ34.ysplit[[i]]$Nage=dQ34.Nage[[i]];

ddQ34 <- do.call("c",dQ34.ysplit)
dQ34 <- NULL
dQ34.ysplit <- NULL


################
## Model
################

ddQ34$ctime = as.numeric(as.character(ddQ34$Year))
nyears=nlevels(ddQ34$Year)

Q34grid <- gridd


modelP3=rep("Year + Gear + s(lon,lat,bs='ds',k=120,m=c(1,0.5)) + s(lon,lat,bs='ds',m=c(1,0.5),k=9,by=Year,id=1) + s(Depth,bs='ds',m=c(1,0),k=6)+s(TimeShotHour,bs='cc',k=6)+s(Ship,bs='re')+offset(log(HaulDur))",length(Q34ages))
modelZ3=rep("Year + Gear + s(lon,lat,bs='ds',k=80,m=c(1,0.5)) + s(lon,lat,bs='ds',k=7,m=c(1,0.5),by=Year,id=1) + s(Depth,bs='ds',m=c(1,0),k=6)+s(TimeShotHour,bs='cc',k=6)+s(Ship,bs='re')+offset(log(HaulDur))",length(Q34ages))


models = list()


## SELECTED MODEL as per benchmark 2021:

models$SIQ34.3 = getSurveyIdx(ddQ34,Q34ages,myids=NULL,predD=Q34grid,cutOff=0.01,fam="Gamma",mc.cores=10,modelZ=modelZ3,modelP=modelP3,nBoot=1000,gamma=1,control=list(trace=TRUE,maxit=10))


########################
## Subpopulation indices
########################

## Cod areas
Viking = subset(Q34grid,subArea == 1)
Northwest = subset(Q34grid,subArea %in% c(2,3))
South = subset(Q34grid,subArea == 4)

models$Viking = redoSurveyIndex(ddQ34,models[[1]],predD=Viking,predfix=NULL,mc.cores=10)
models$Northwest = redoSurveyIndex(ddQ34,models[[1]],predD=Northwest,predfix=NULL,mc.cores=10)
models$South = redoSurveyIndex(ddQ34,models[[1]],predD=South,predfix=NULL,mc.cores=10)

save(ddQ34, models, file="model/Q34models.Rdata")


#########################
## Retrospective analysis
#########################

SIretro = retro.surveyIdx(models[[1]],ddQ34,predD=Q34grid,npeels=3,control=list(trace=TRUE,maxit=10),mc.cores=10)
save(SIretro, file="model/Q34retro.Rdata")



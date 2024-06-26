library(surveyIndex)

load("data/Q1preprocessed.RData")
load("boot/data/commonGrid.RData")

xtabs(NoAtALK~Year+Age,data=dQ1[[1]])
Q1ages=1:7


########################
## ALK
########################

dQ1 = addSpectrum(dQ1,by=1)
dQ1.ysplit = split(dQ1,dQ1$Year)

## Declare settings for ALK model
mf = NULL
ack=TRUE;
useBICs=TRUE;
varCofs=FALSE;
maxKs=50;

dQ1.ALK= mclapply(dQ1.ysplit,fitALK,minAge=min(Q1ages),maxAge=max(Q1ages),autoChooseK=ack,useBIC=useBICs,varCof=varCofs,maxK=maxKs,mc.cores=10)

dQ1.Nage=mclapply(dQ1.ALK,predict,mc.cores=10)
for(i in 1:length(dQ1.ALK)) dQ1.ysplit[[i]]$Nage=dQ1.Nage[[i]];

ddQ1 <- do.call("c",dQ1.ysplit)
dQ1 <- NULL
dQ1.ysplit <- NULL


################
## Model
################

ddQ1$ctime = as.numeric(as.character(ddQ1$Year))
nyears=nlevels(ddQ1$Year)

# # Run for the benchmark only
# grid <- getGrid(ddQ1, nLon = 55)
# gridd <- subset(ddQ1[[2]], haul.id %in% grid[[3]])
# save(grid, gridd, file="data/commonGrid.Rdata")

Q1grid <- gridd


modelP3=rep("Year + s(lon,lat,bs='ds',k=120,m=c(1,0.5)) + s(lon,lat,bs='ds',m=c(1,0.5),k=9,by=Year,id=1) + s(Depth,bs='ds',m=c(1,0),k=6)+s(TimeShotHour,bs='cc',k=6)+s(Ship,bs='re')+offset(log(HaulDur))",length(Q1ages))
modelZ3=rep("Year + s(lon,lat,bs='ds',k=80,m=c(1,0.5)) + s(lon,lat,bs='ds',k=7,m=c(1,0.5),by=Year,id=1) + s(Depth,bs='ds',m=c(1,0),k=6)+s(TimeShotHour,bs='cc',k=6)+s(Ship,bs='re')+offset(log(HaulDur))",length(Q1ages))


models = list()


## SELECTED MODEL as per benchmark 2021:

system.time(models$SIQ1.3 <- getSurveyIdx(ddQ1,Q1ages,myids=NULL,predD=Q1grid,cutOff=0.01,fam="Gamma",mc.cores=10,modelZ=modelZ3,modelP=modelP3,nBoot=1000,gamma=1,control=list(trace=TRUE,maxit=10)))


########################
## Subpopulation indices
########################

## Cod areas
Viking = subset(Q1grid,subArea == 1)
Northwest = subset(Q1grid,subArea %in% c(2,3))
South = subset(Q1grid,subArea == 4)

pdf("model/subareas.pdf")
plot(Q1grid$lon,Q1grid$lat,type="n")
points(Viking$lon,Viking$lat,pch=16,col='red')
points(Northwest$lon,Northwest$lat,pch=16,col='blue')
points(South$lon,South$lat,pch=16,col='green')
map("worldHires", fill = TRUE, plot = TRUE, add = TRUE, col = grey(0.5))
legend("topright",legend=c("Viking","Northwest","South"),pch=16,col=c("red","blue","green"),bg="white")
dev.off()

models$Viking = redoSurveyIndex(ddQ1,models[[1]],predD=Viking,predfix=NULL,mc.cores=10)
models$Northwest = redoSurveyIndex(ddQ1,models[[1]],predD=Northwest,predfix=NULL,mc.cores=10)
models$South = redoSurveyIndex(ddQ1,models[[1]],predD=South,predfix=NULL,mc.cores=10)

save(ddQ1, models, file="model/Q1models.Rdata")


#########################
## Retrospective analysis
#########################

SIretro = retro.surveyIdx(models[[1]],ddQ1,predD=Q1grid,npeels=3,control=list(trace=TRUE,maxit=10),mc.cores=10)
save(SIretro, file="model/Q1retro.Rdata")


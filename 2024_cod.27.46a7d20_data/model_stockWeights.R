#library(surveyIndex)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(icesTAF)
library(stockassessment)
library(DATRAS)

load("data/Q1preprocessed.RData")
source("utilities.R")

cur.yr <-2024
first.yr <-1983

Q1ages=1:7
plusgp <-max(Q1ages)

Nb <-500 # number of bootstrap replicates

dQ1 = addSpectrum(dQ1,by=1)

ddQ1 <-dQ1
dQ1 <-NULL
ddQ1$ctime = as.numeric(as.character(ddQ1$Year))

########################
## Weight-at-age
########################

ddQ1$subAreaName <- NA
ddQ1$subAreaName[ddQ1$subArea==1] <-"Viking"
ddQ1$subAreaName[ddQ1$subArea %in% c(2,3)] <-"Northwest"
ddQ1$subAreaName[ddQ1$subArea==4] <-"South"

ddQ1$YearsubArea <- factor(paste(ddQ1$Year,ddQ1$subAreaName,sep=":"))

ddQ1=addSpectrum(ddQ1,by=1)

YAsplitQ1= split(ddQ1,ddQ1$YearsubArea)

xtabs(is.na(IndWgt)~Year+Age,data=ddQ1[[1]])
xtabs(!is.na(IndWgt)~Year+Age,data=ddQ1[[1]])

ddQ1[[1]]$SubArea <-do.call("rbind",strsplit(as.character(ddQ1[[2]]$YearsubArea[match(ddQ1[[1]]$haul.id,ddQ1[[2]]$haul.id)]),":"))[,2]

xtabs(is.na(IndWgt)~Year+SubArea,data=ddQ1[[1]])
xtabs(!is.na(IndWgt)~Year+SubArea,data=ddQ1[[1]])

xtabs(is.na(IndWgt)~Year+Survey,data=ddQ1[[1]])
xtabs(!is.na(IndWgt)~Year+Survey,data=ddQ1[[1]])

xtabs(is.na(IndWgt)~Year+Country,data=ddQ1[[1]])
xtabs(!is.na(IndWgt)~Year+Country,data=ddQ1[[1]])

mkdir("model/tables")
# Number of hauls per survey
write.csv(xtabs(!is.na(haul.id)~Year+Survey,data=ddQ1[["HH"]]),file="model/tables/survey.hauls.csv")
# Number of length and weight measurements per survey
write.csv(xtabs(!is.na(IndWgt)~Year+Survey,data=ddQ1[[1]]),file="model/tables/survey.lwt.samples.csv")

#  Number of hauls per sub stock
write.csv(xtabs(!is.na(haul.id)~Year+subAreaName,data=ddQ1[["HH"]]),file="model/tables/subarea.hauls.csv")


##########################################################################
#  Split by area only - given no trend over time
#  Include only 2011-2022 which includes Scottish data
#
##########################################################################
dat <-subset(ddQ1, Year %in% 2011:2022)

dat[[1]] <-dat[[1]][!is.na(dat[[1]]$IndWgt),]
dat[[1]] <-subset(dat[[1]],IndWgt>0)
AsplitQ1 <-split(dat,dat$subAreaName)

# The models by area
wtlenMod <-lapply(AsplitQ1,weight.len.mod)

sum.wtlen <-data.frame(do.call("rbind",lapply(wtlenMod,function(x){summary(x)$coef[,1]})))
colnames(sum.wtlen) <-c("loga","b")

write.csv(sum.wtlen, file=paste0("model/tables/lwt.params.by.area.csv"))

############################################
#  Bootstrapped estimates of weight at age
#############################################

WAA.boot <-lapply(YAsplitQ1,boot.weight.at.age,minAge=min(Q1ages),maxAge=max(Q1ages),WALmods=wtlenMod,Nboot=Nb)

#  Output the substock mean weights by subarea ONLY 1983 onwards

###### Mean and variance of bootstrapped distribution

WAA.mean <-bind_rows(lapply(WAA.boot,function(x){apply(subset(x,rownames(x)!="Value")/1000,2,mean,na.rm=TRUE)}),.id='YearSubArea')
nms <-do.call("rbind",strsplit(WAA.mean$YearSubArea,":"))
WAA.mean$Year <-as.numeric(nms[,1])
WAA.mean$SubArea <-nms[,2]

WAA.sd <-bind_rows(lapply(WAA.boot,function(x){apply(subset(x,rownames(x)!="Value")/1000,2,sd,na.rm=TRUE)}),.id='YearSubArea')
WAA.sd$Year <-as.numeric(nms[,1])
WAA.sd$SubArea <-nms[,2]

names(WAA.mean)[grep("Age",names(WAA.mean))] <-as.character(Q1ages)
names(WAA.sd)[grep("Age",names(WAA.sd))] <-as.character(Q1ages)
WAA.mean <-WAA.mean %>% select(-YearSubArea)
WAA.sd <-WAA.sd %>% select(-YearSubArea)

save(WAA.boot,WAA.mean,WAA.sd, file="model/sw.Rdata")

## Preprocess data, write TAF data tables

## Before: Assessment input files
## After: data.R

library(icesTAF)
library(stockassessment)

mkdir("data")

source("utilities.R")

maxAge <- 7
maxYear <- Inf
minYear <- 1983

path <- "boot/data"

################################################################################
#### Shared catch data
################################################################################

## CN
cn_WS <- read.ices(file.path(path,"cn_WS.dat"))
cn_NS <- read.ices(file.path(path,"cn_NS.dat"))

keyFleetStock <- matrix(0,0,3, dimnames=list(NULL,c("Northwest","South","Viking")))

cn <- cutSum(cn_NS) + cutSum(cn_WS)
keyFleetStock <- rbind(keyFleetStock,"Catch"=c(1,1,1))

## CW
cw_WS <- read.ices(file.path(path,"cw_WS.dat"))
cw_NS <- read.ices(file.path(path,"cw_NS.dat"))
cw <- (na2zero(cutMean(cw_NS,cn_NS) * cutSum(cn_NS)) + na2zero(cutMean(cw_WS,cn_WS) * cutSum(cn_WS))) / cn 

## LN, first as landing in numbers
ln_WS <- read.ices(file.path(path,"lf_WS.dat")) * cn_WS
ln_NS <- read.ices(file.path(path,"lf_NS.dat"))
ln <- cutSum(ln_NS) + cutSum(ln_WS)

## LF, Convert to fraction
lf <- (cutSum(ln_WS) + cutSum(ln_NS)) / cn

## LW
lw_WS <- read.ices(file.path(path,"lw_WS.dat"))
lw_NS <- read.ices(file.path(path,"lw_NS.dat"))
lw <- (na2zero(cutMean(lw_NS,ln_NS) * cutSum(ln_NS)) + na2zero(cutMean(lw_WS,ln_WS) * cutSum(ln_WS))) / ln 

## DN (discard numbers)
dn_WS <- cn_WS - ln_WS
dn_NS <- cn_NS - ln_NS 
dn <- cutSum(dn_NS) + cutSum(dn_WS)

## DW
dw_WS <- read.ices(file.path(path,"dw_WS.dat"))
dw_NS <- read.ices(file.path(path,"dw_NS.dat"))
dw <- (na2zero(cutMean(dw_NS,dn_NS) * cutSum(dn_NS)) + na2zero(cutMean(dw_WS,dn_WS) * cutSum(dn_WS))) / dn
dw[dn==0] <- 0


################################################################################
#### "Stock-wise" biological input
################################################################################

mo <- list(Northwest = cutDrop(read.ices(file.path(path,"mo_northwest.dat"))),
           South = cutDrop(read.ices(file.path(path,"mo_south.dat"))),
           Viking = cutDrop(read.ices(file.path(path,"mo_viking.dat"))))

nm0 <- cutDrop(read.ices(file.path(path,"nm_2023.dat")))
# nm0 <- read.ices(file.path(path,"nm_2020.dat"))
# nm0[(nrow(nm0)-4):nrow(nm0),] <- matrix(nm0["2019",],nrow=5,ncol=15, byrow=TRUE)
# nm0 <- cutDrop(nm0)

## No stock-wise information
nm <- list(Northwest = nm0,
           South = nm0,
           Viking = nm0)

## PF (all 0, handled by setup SAM data)
pf <- NULL
## PM (all 0, handled by setup SAM data)
pm <- NULL

## SW
sw <- list(Northwest = cutDrop(read.ices(file.path(path,"sw_northwest.dat"))),
           South = cutDrop(read.ices(file.path(path,"sw_south.dat"))),
           Viking = cutDrop(read.ices(file.path(path,"sw_viking.dat"))))

## Surveys
#### Write to valid file
sf <- tempfile("survey",fileext=".dat")
cat("survey.dat auto written",sep="\n",file = sf)
cat("One more line that isn't used",sep="\n",file = sf, append=TRUE)
survfold <- ""
sf2 <- ""

cat(readLines(sprintf(file.path(path,"%s/q34index.txt"),survfold)),sep="\n",file=sf,append=TRUE)
keyFleetStock <- rbind(keyFleetStock,"Survey_Q34"=c(1,1,1))
## Q1 indices
cat(readLines(sprintf(file.path(path,"%s/NorthwestQ1.txt"),paste0(survfold,sf2))),sep="\n",file=sf,append=TRUE)
keyFleetStock <- rbind(keyFleetStock,"Survey_Q1_NW"=c(1,0,0))
cat(readLines(sprintf(file.path(path,"%s/SouthQ1.txt"),paste0(survfold,sf2))),sep="\n",file=sf,append=TRUE)
keyFleetStock <- rbind(keyFleetStock,"Survey_Q1_SO"=c(0,1,0))
cat(readLines(sprintf(file.path(path,"%s/VikingQ1.txt"),paste0(survfold,sf2))),sep="\n",file=sf,append=TRUE)
keyFleetStock <- rbind(keyFleetStock,"Survey_Q1_VI"=c(0,0,1))

cat(readLines(sprintf(file.path(path,"%s/Northwest34recruits.txt"),survfold)),sep="\n",file=sf,append=TRUE)
keyFleetStock <- rbind(keyFleetStock,"Survey_Rec_NW"=c(1,0,0))
cat(readLines(sprintf(file.path(path,"%s/South34recruits.txt"),survfold)),sep="\n",file=sf,append=TRUE)
keyFleetStock <- rbind(keyFleetStock,"Survey_Rec_SO"=c(0,0.75,0))
cat(readLines(sprintf(file.path(path,"%s/Viking34recruits.txt"),survfold)),sep="\n",file=sf,append=TRUE)        
keyFleetStock <- rbind(keyFleetStock,"Survey_Rec_VI"=c(0,0.25,1))

surveys <- stockassessment:::read.surveys(sf) #Note: this does not give a proper plus group index
names(surveys) <- gsub(" ","_",names(surveys))
surveys <- lapply(surveys, cutDrop)

## Years are already forward shifted - ensure times is c(0,0) and age is 1
colnames(surveys[[5]]) <- 1
attr(surveys[[5]],"times") <- c(0,0)
colnames(surveys[[6]]) <- 1
attr(surveys[[6]],"times") <- c(0,0)
colnames(surveys[[7]]) <- 1
attr(surveys[[7]],"times") <- c(0,0)

## Shared Q34 index
sdQ34 <- as.matrix(read.table(sprintf(file.path(path,"%s/q34sds.txt"),survfold)))
colnames(sdQ34) <- gsub("X","",colnames(sdQ34))
attr(surveys[[1]],"weight") <- cutDrop(1/sdQ34^2, as.numeric(maxAge))
## NW Q1
sdQ1_NW <- as.matrix(read.table(sprintf(file.path(path,"%s/NorthwestQ1sds.txt"),paste0(survfold,sf2))))
colnames(sdQ1_NW) <- gsub("X","",colnames(sdQ1_NW))
attr(surveys[[2]],"weight") <- cutDrop(1/sdQ1_NW^2, as.numeric(maxAge))
## SO Q1
sdQ1_SO <- as.matrix(read.table(sprintf(file.path(path,"%s/SouthQ1sds.txt"),paste0(survfold,sf2))))
colnames(sdQ1_SO) <- gsub("X","",colnames(sdQ1_SO))
attr(surveys[[3]],"weight") <- cutDrop(1/sdQ1_SO^2, as.numeric(maxAge))
## VI Q1
sdQ1_VI <- as.matrix(read.table(sprintf(file.path(path,"%s/VikingQ1sds.txt"),paste0(survfold,sf2))))
colnames(sdQ1_VI) <- gsub("X","",colnames(sdQ1_VI))
attr(surveys[[4]],"weight") <- cutDrop(1/sdQ1_VI^2, as.numeric(maxAge))

## NW Rec
sdQRec_NW <- as.matrix(read.table(sprintf(file.path(path,"%s/Northwest34Rsds.txt"),survfold)))
colnames(sdQRec_NW) <- gsub("X","",colnames(sdQRec_NW))
colnames(sdQRec_NW) <- 1
attr(surveys[[5]],"weight") <- cutDrop(1/sdQRec_NW^2)
## SO Rec
sdQRec_SO <- as.matrix(read.table(sprintf(file.path(path,"%s/South34Rsds.txt"),survfold)))
colnames(sdQRec_SO) <- gsub("X","",colnames(sdQRec_SO))
colnames(sdQRec_SO) <- 1
attr(surveys[[6]],"weight") <- cutDrop(1/sdQRec_SO^2)
## VI Rec
sdQRec_VI <- as.matrix(read.table(sprintf(file.path(path,"%s/Viking34Rsds.txt"),survfold)))
colnames(sdQRec_VI) <- gsub("X","",colnames(sdQRec_VI))
colnames(sdQRec_VI) <- 1
attr(surveys[[7]],"weight") <- cutDrop(1/sdQRec_VI^2)

auxData <- NULL

e <- new.env()
load(file.path(path,"WGNSSK2024_quartYearProp.RData"),e)
## Remove quarter 4 (no stochasticity given the first three quarters)
e$quartYearProp <- e$quartYearProp[e$quartYearProp$year <= max(as.numeric(rownames(ln)))   &
                                     as.numeric(e$quartYearProp$year) >= minYear &
                                     as.numeric(e$quartYearProp$year) <= maxYear,]
qyp0 <- do.call("rbind",lapply(split(e$quartYearProp,e$quartYearProp$year),function(x){ x$prop <- log(x$prop/x$prop[as.numeric(x$quarter)==4]);x}))
qyp <- qyp0[as.numeric(qyp0$quarter) < 4,]
qtimes <- c(0,0.25,0.5,0.75,1)
auxData <- list(as.matrix(data.frame(fleetType = 80,
                                     year = as.numeric(qyp$year),
                                     age = -1, ## Combined for all ages
                                     logobs = qyp$prop,
                                     catchFleet = 1,
                                     seasonStart = qtimes[as.numeric(qyp$quarter)],
                                     seasonEnd = qtimes[as.numeric(qyp$quarter)+1],
                                     catchType = 3, #Landed weight
                                     seasonNumber = as.numeric(qyp$quarter),
                                     numberOfSeasons = max(as.numeric(qyp$quarter))+1)))
keyFleetStock <- rbind(keyFleetStock,"AUX_SeasonProp"=c(1,1,1))

e <- new.env()
load(file.path(path,"WGNSSK2024_SubstockProportions.RData"),e)
## Remove Viking (no stochasticity given the two others) and years without landing information
## Transform proportions to a useful scale for missing obs
apd0 <- do.call("rbind",lapply(split(e$substockProps,e$substockProps$year),function(x){ x$prop <- log(x$prop/x$prop[x$areaPopulation=="Viking"]);x}))
apd <- apd0[apd0$areaPopulation != "Viking" &
              as.numeric(apd0$year) <= max(as.numeric(rownames(ln_NS))) &
              as.numeric(apd0$year) >= min(as.numeric(rownames(ln_NS)))  &
              as.numeric(apd0$year) >= minYear &
              as.numeric(apd0$year) <= maxYear,]
auxData <- append(auxData,list(as.matrix(data.frame(fleetType = 90,
                                                    year = as.numeric(apd$year),
                                                    age = -1, ## Combined for all ages
                                                    logobs = apd$prop,
                                                    catchFleet = 0,
                                                    seasonStart = 0,
                                                    seasonEnd = 0.25,
                                                    catchType = 3, #Landed weight
                                                    stockNumber = match(gsub("(^[A-Z])(.+)","\\1",apd$areaPopulation),c("N","S","V"))))))
keyFleetStock <- rbind(keyFleetStock,"AUX_SubNew"=c(1,1,1))


Stock_data <- lapply(as.list(c("Northwest","South","Viking")), function(nam){
  cat(nam,"\n")
  tmp <-  setup.sam.data(surveys=surveys,
                         residual.fleet=cn,
                         aux.fleets=auxData,
                         prop.mature=mo[[nam]], 
                         stock.mean.weight=sw[[nam]], 
                         catch.mean.weight=cw, 
                         dis.mean.weight=dw, 
                         land.mean.weight=lw,
                         ## prop.f=pf, 
                         ## prop.m=pm, 
                         natural.mortality=nm[[nam]], 
                         land.frac=lf)
  ##tmp$logobs <- NA (observations are needed for defconf)
  tmp
})
names(Stock_data) <- c("Northwest","South","Viking")

sharedDatIn <- setup.sam.data(surveys = surveys,
                              residual.fleet = cn,
                              aux.fleets=auxData,                               
                              prop.mature=mo[[1]],
                              stock.mean.weight=sw[[1]],
                              catch.mean.weight=cw,
                              dis.mean.weight=dw,
                              land.mean.weight=lw,
                              prop.f=pf[[1]],
                              prop.m=pm[[1]],
                              natural.mortality=nm[[1]],
                              land.frac=lf)

## Reduce to needed information,
## add keyFleetStock (row is fleet, col is stock, value is proportion of stock covered by fleet. Use NA to estimate the valye)
sharedDat <- c(list(hasSharedObs = 1L, maxAgePlusGroup = 1),
               sharedDatIn[c("aux","auxData","fleetTypes","logobs","idx1","idx2","idxCor","noYears","noFleets","weight","corList")],
               list(keyFleetStock = keyFleetStock,
                    ## Method to combine sub-stock variances
                    ## 0: Delta method A; 1: Delta method B; 2: Average (as babysam)
                    covCombine = 2L
               ))
## Add fleetNames
attr(sharedDat, "fleetNames") <- attr(sharedDatIn,"fleetNames")

dat <- list(Stock_data = Stock_data,
            Shared_data = sharedDat)   

save(dat, file="data/data.Rdata")
  
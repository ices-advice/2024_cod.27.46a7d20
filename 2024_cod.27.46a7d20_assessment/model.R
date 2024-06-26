## Run analysis, write model results

## Before: data.Rdata, conf.Rdata
## After: model.Rdata

library(icesTAF)
library(stockassessment)
library(multiStockassessment)

mkdir("model")

source("utilities.R")

# 1 Fit multistock SAM
load('data/data.Rdata')
load('config/conf.Rdata')
system.time(fit <- runFit(conf, run = TRUE))
# AIC: -2681.798
save(fit, file="model/model.Rdata")

## 2 Residuals
RES <- residuals(fit)
save(RES, file="model/residuals.Rdata")

## 3 Retrospective runs
system.time(RETRO <- retro(fit, 5, 1, initializePars = TRUE, initializeRandomEffects = TRUE))
stockassessment::mohn(RETRO, addTotal = TRUE)
stockassessment::mohn(RETRO, addTotal=TRUE, lag = 1)
save(RETRO, file="model/retro.Rdata")

## 4 Reference points
mkdir("model/refpts")
# Extract sd of ln(SSB)
v1 <- summary(attr(fit,"m_sdrep"))
v2 <- v1[grepl(".+_logssb",rownames(v1)),]
sigmas <- lapply(split(v2[,2],rownames(v2)),tail,n=1)
save(sigmas, file="model/refpts/sigmas.Rdata")

source.taf("model_northwesternRPs.R")
source.taf("model_southernRPs.R")
source.taf("model_vikingRPs.R")

## 5 Forecasts
load("model/model.Rdata")
load("model/refpts/northwestern/NW_1997_type2_refPoints.Rdata")
load("model/refpts/southern/S_1997_type2_5yr_refPoints.Rdata")
load("model/refpts/viking/V_1997_typeM_refPoints.Rdata")
NW <- NW_1997_type2_refPoints
SO <- S_1997_type2_5yr_refPoints
VI <- V_1997_typeM_refPoints

cores <- 1
DF <- FALSE        # deterministicF 
UNLC <- !DF        # useNonLinearityCorrection 
UFH <- TRUE        # useFHessian
PNF <- FALSE       # processNoiseF
FFD <- FALSE       # fixedFDeviation
NOSIM <- 1000       # Needs changing to 1000 

## Average biological input over last three years
aveYears <- lapply(fit,function(x) tail(head(x$data$years,-1),3))
## Use recruitment from 1998 onwards instead of log-random walk
recYears <- lapply(fit,function(x) 1998:max(x$data$years))

AdviceYear <- 2025

## Reference points
Fmsy <- c("Northwest" = as.numeric(round(NW['Fmsy'],3)), "South" = as.numeric(round(SO['Fmsy'],3)), "Viking" = as.numeric(round(VI['Fmsy'],3)))
Flower <- c("Northwest" = as.numeric(round(NW['Fmsylower'],3)), "South" = as.numeric(round(SO['Fmsylower'],3)), "Viking" = as.numeric(round(VI['Fmsylower'],3)))
Fupper <- c("Northwest" = as.numeric(round(NW['Fmsyupper'],3)), "South" = as.numeric(round(SO['Fmsyupper'],3)), "Viking" = as.numeric(round(VI['Fmsyupper'],3)))
MSYBtrigger <- c("Northwest" = as.numeric(round(NW['MSYbtrigger'])), "South" = as.numeric(round(SO['MSYbtrigger'])), "Viking" = as.numeric(round(VI['MSYbtrigger'])))
Blim <- c("Northwest" = as.numeric(round(NW['Blim'])), "South" = as.numeric(round(SO['Blim'])), "Viking" = as.numeric(round(VI['Blim'])))
Fpa <- c("Northwest" = as.numeric(round(NW['Fpa'],3)), "South" = as.numeric(round(SO['Fpa'],3)), "Viking" = as.numeric(round(VI['Fpa'],3)))
Flim <- c("Northwest" = as.numeric(round(NW['Flim'],3)), "South" = as.numeric(round(SO['Flim'],3)), "Viking" = as.numeric(round(VI['Flim'],3)))


## List to store results
FC <- list()

## F = FMSY
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",
                 rep(sprintf("F=%f",Fmsy[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "F = F~MSY~",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

## F = FmsyLower
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",
                 rep(sprintf("F=%f",Flower[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "F = F~MSY lower~",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

## F = FmsyUpper
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",
                 rep(sprintf("F=%f",Fupper[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "F = F~MSY upper~",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

## F = 0
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",rep("F=0.000001",2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "F = 0",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

## F = Fpa
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",rep(sprintf("F=%f",Fpa[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "F = F~pa~",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

## F = Flim
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",rep(sprintf("F=%f",Flim[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "F = F~lim~",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))



## SSB = Blim
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",
                 rep(sprintf("SSB=%f",Blim[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = sprintf("SSB (%d) = B~lim~",AdviceYear+1),
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

## SSB = MSY Btrigger
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",
                 rep(sprintf("SSB=%f",MSYBtrigger[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = sprintf("SSB (%d) = MSY B~trigger~",AdviceYear+1),
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

## F = F2024
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",rep("F=1*",2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "F = F~2024~",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

## SSB(2026) = SSB(2025)
SSB <- unlist(lapply(FC[[1]]$result,function(x)attr(x,"shorttab")[3,3]))
cstr <- vector("list",3)
for(i in 1:3)
  cstr[[i]] <- c("F=1*",
                 rep(sprintf("SSB=%f",SSB[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "SSB(2026)=SSB(2025)",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))

#save(FC, file="model/forecast.Rdata")

## ----------------------------------------------

# HCRs (if applicable)

SSB <- unlist(lapply(FC[[1]]$result,function(x)attr(x,"shorttab")[3,3]))

if(any(SSB<MSYBtrigger)){
  
  ## ICES Advice rule, FMSY
  fval <- pmin(Fmsy*SSB/MSYBtrigger, Fmsy)
  cstr <- vector("list",3)
  for(i in 1:3)
    cstr[[i]] <- c("F=1*",
                   rep(sprintf("F=%f",fval[i]),2))
  set.seed(12345)
  FC[[length(FC)+1]] <- list(name = "MSY approach: F~MSY~ x SSB (2025) / MSY B~trigger~",
                             constraints = cstr,
                             result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                    deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                    nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))
  
  ## ICES Advice rule, Flower
  fval <- pmin(Flower*SSB/MSYBtrigger, Flower)
  cstr <- vector("list",3)
  for(i in 1:3)
    cstr[[i]] <- c("F=1*",
                   rep(sprintf("F=%f",fval[i]),2))
  set.seed(12345)
  FC[[length(FC)+1]] <- list(name = "F = F~MSY lower~ x SSB (2025) / MSY B~trigger~",
                             constraints = cstr,
                             result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                    deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                    nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))
}

save(FC, file="model/forecast.Rdata")


## Common F reduction
SSB <- unlist(lapply(FC[[1]]$result,function(x)attr(x,"shorttab")[3,3]))
fval <- pmin(Fmsy*SSB/MSYBtrigger, Fmsy)
Fmult <- fval[[2]] / attr(FC[[1]]$result[[2]], "shorttab")[1,2]
cstr <- vector("list",3)
for(i in c(1,3))
  cstr[[i]] <- c("F=1*",
                 rep(paste0(sprintf("F=%f",Fmult),'*'),2))
for(i in 2)
  cstr[[i]] <- c("F=1*",
                 rep(sprintf("F=%f",fval[i]),2))
set.seed(12345)
FC[[length(FC)+1]] <- list(name = "F = 0.39F~2024~",
                           constraints = cstr,
                           result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                  deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                  nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))
save(FC, file="model/forecast.Rdata")


## Common F reduction range
FCR1 <- list()
SSB <- unlist(lapply(FC[[1]]$result,function(x)attr(x,"shorttab")[3,3]))
fval <- pmin(Fmsy*SSB/MSYBtrigger, Fmsy)
# Which multipliers give an F < Fmsy for the southern substock?
mults <- seq(1,2,by=0.1)[which(fval[[2]]*seq(1,2,by=0.1) < Fmsy[[2]])]
for (m in mults){
  Fmult <- m * fval[[2]] / attr(FC[[1]]$result[[2]], "shorttab")[1,2]
  cstr <- vector("list",3)
  for(i in c(1,3))
    cstr[[i]] <- c("F=1*",
                   rep(paste0(sprintf("F=%f",Fmult),'*'),2))
  for(i in 2)
    cstr[[i]] <- c("F=1*",
                   rep(sprintf("F=%f",m*fval[i]),2))
  set.seed(12345)
  FCR1[[length(FCR1)+1]] <- list(name = "F = mF~2024~",
                             constraints = cstr,
                             result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                    deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                    nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))
}

save(FC, FCR1, file="model/forecast.Rdata")


# Forecasts in F increments
FCR2 <- list()

for (fval in c(0.000001, seq(0.01,round(max(Fupper), digits=2),0.01))){
  
  cstr <- vector("list",3)
  for(i in 1:3)
    cstr[[i]] <- c("F=1*",
                   rep(sprintf("F=%f",fval),2))
  set.seed(12345)
  FCR2[[length(FCR2)+1]] <- list(name = "Frange",
                                 constraints = cstr,
                                 result = modelforecast(fit, cstr, ave.years=aveYears, rec.years=recYears,year.base=2023, useModelLastN=TRUE, resampleFirst = TRUE, splitLD = TRUE,
                                                        deterministicF = DF, useNonLinearityCorrection= UNLC, processNoiseF = PNF, fixedFDeviation = FFD,
                                                        nosim=NOSIM,useFHessian=UFH, ncores = cores, progress = TRUE))
}

save(FC, FCR1, FCR2, file="model/forecast.Rdata")

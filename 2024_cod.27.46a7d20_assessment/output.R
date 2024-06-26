## Extract results of interest, write TAF output tables

## Before:
## After:

library(icesTAF)
library(FLCore)
library(stockassessment)
library(multiStockassessment)
library(icesAdvice)
library(xlsx)
library(reshape2)

mkdir("output")

## 1 FLR stock objects
load("model/model.Rdata")
source(system.file("flr_convenience.R",package="multiStockassessment"))

# Model inputs (excluding shared catches)
stks_dat <- as.FLStock.msam(fit, biopar=FALSE)

ay <- range(stks_dat[[1]])["maxyear"]
ny <- length(1983:ay)

for (i in 1:3){
  # Fbar range
  stks_dat[[i]]@range[c("minfbar","maxfbar")] <- c(2,4)
  
  # Get rid of assessment estimates
  stks_dat[[i]]@stock.n <- FLQuant(matrix(0, nrow = 7, ncol=ny), dimnames=list(age=1:7,year=1983:range(stks_dat[[i]])["maxyear"]))
  stks_dat[[i]]@stock <- computeStock(stks_dat[[i]])
  stks_dat[[i]]@harvest <- FLQuant(matrix(0, nrow = 7, ncol=ny), dimnames=list(age=1:7,year=1983:range(stks_dat[[i]])["maxyear"]))
  
  # Name & description
  substock <- name(stks_dat[[i]])
  name(stks_dat[[i]]) <- paste("Northern Shelf cod (cod.27.46a7d20).", substock, "sub-stock (WGNSSK 2024)")
  desc(stks_dat[[i]]) <- paste("cod.27.46a7d20", substock, "- FLStock containing observational data")
  
  # Units
  nmes <- names(units(stks_dat[[i]]))
  unit.lst <- as.list(c(rep(c("tonnes","thousands","kg"),4),rep("NA",2),"f",rep("NA",2)))
  names(unit.lst) <- nmes
  units(stks_dat[[i]]) <- unit.lst
  
  print(plot(stks_dat[[i]]))
}

# Model outputs
stks <- as.FLStock.msam(fit, predicted=TRUE)

for (i in 1:3){
  # Fbar range
  stks[[i]]@range[c("minfbar","maxfbar")] <- c(2,4)
  
  # Catches
  stks[[i]]@catch <- computeCatch(stks[[i]])
  stks[[i]]@landings <- computeLandings(stks[[i]])
  stks[[i]]@discards <- computeDiscards(stks[[i]])
  
  # Name & description
  substock <- name(stks[[i]])
  name(stks[[i]]) <- paste("Northern Shelf cod (cod.27.46a7d20).", substock, "sub-stock (WGNSSK 2024)")
  desc(stks[[i]]) <- paste("cod.27.46a7d20", substock, "- FLStock containing stock assessment model outputs")
  
  # Units
  nmes <- names(units(stks[[i]]))
  unit.lst <- as.list(c(rep(c("tonnes","thousands","kg"),4),rep("NA",2),"f",rep("NA",2)))
  names(unit.lst) <- nmes
  units(stks[[i]]) <- unit.lst
  
  print(plot(stks[[i]]))
}

save(stks_dat, file="output/cod46a7d20_FLstockObjects_initial_2024.Rdata")
save(stks, file="output/cod46a7d20_FLstockObjects_assessment_2024.Rdata")

## 2 Advice
load('model/model.Rdata')
summ <- summary(fit, returnList=TRUE)
for(i in 1:length(fit)){
  summ[[i]] <- as.data.frame(cbind(summ[[i]], rbind(catchtable(fit, returnList=TRUE)[[i]], rep(NA,3))))
  colnames(summ[[i]])[c(10)] <- c("Catches")
  summ[[i]][,10:12] <- round(summ[[i]][,10:12])
  summ[[i]] <- summ[[i]][,c(2,1,3,5,4,6,11,10,12,8,7,9)]
  summ[[i]][nrow(summ[[i]]), 10:12] <- NA
  for (col in 10:12){
    summ[[i]][,col] <- icesRound(as.numeric(summ[[i]][,col]))
  }
  write.xlsx(summ[[i]], file="output/summary.xlsx", sheetName = names(summ)[[i]], append=(i!=1))
}

## 3 Change in advice
source('utilities_multistockPlots.R')

# load last year's assessment and forecast
load("boot/data/WGNSSK2023_model.Rdata") # last year's assessment
prev_ass <- fit
load("boot/data/WGNSSK2023_forecast.Rdata")  # last year's forecast
prev_FC1 <- FC[[13]]

# load this year's assessment and forecast
load("model/model.RData") # this year's assessment
load("model/forecast.RData") # this year's forecast

substocks <- c("Northwest","South","Viking")
ay <- 2024

# Forecast_assumptions
assumptions <- list()
for (s in 1:length(fit)){
  Last <- data.frame(Year=names(getFCTab(prev_FC1, F='rec', s)[,1]),
                     Recruitment=getFCTab(prev_FC1, F='rec', s)[,1],
                     SSB=getFCTab(prev_FC1, F='ssb', s)[,1][1:3],
                     Fbar=getFCTab(prev_FC1, F='fbar', s)[,1],
                     'Total catch'=getFCTab(prev_FC1, F='catch', s)[,1],
                     Type=c("Data year", "Intermediate year", "Advice year"),
                     Source=paste("WGNSSK", ay-1))
  
  This <- data.frame(Year=names(getFCTab(FC[[13]], F='rec', s)[,1]),
                     Recruitment=getFCTab(FC[[13]], F='rec', s)[,1],
                     SSB=getFCTab(FC[[13]], F='ssb', s)[,1][1:3],
                     Fbar=getFCTab(FC[[13]], F='fbar', s)[,1],
                     'Total catch'=getFCTab(FC[[13]], F='catch', s)[,1],
                     Type=c("Data year", "Intermediate year", "Advice year"),
                     Source=paste("WGNSSK", ay))
  
  assumptions[[s]] <- rbind(Last,This)
  assumptions[[s]] <- melt(assumptions[[s]], id.vars = c("Source","Year","Type"))
}

write.xlsx(assumptions[[1]],file="output/Forecast_assumptions.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(assumptions[[2]],file="output/Forecast_assumptions.xlsx", sheetName=substocks[2], 
           row.names=F,append=TRUE)
write.xlsx(assumptions[[3]],file="output/Forecast_assumptions.xlsx", sheetName=substocks[3], 
           row.names=F,append=TRUE)

# Forecast_stockwts
stockwts <- list()
sway1 <- multiStockassessment:::splitParameter(exp(attr(prev_ass,"m_pl")$logSW))
sway2 <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logSW))
for (s in 1:length(fit)){
  Last <- data.frame(Age=1:7,
                     'Data year'=sway1[[s]][(nrow(sway1[[s]])-11),],
                     'Intermediate year'=sway1[[s]][(nrow(sway1[[s]])-10),],
                     'Advice year'=sway1[[s]][(nrow(sway1[[s]])-9),],
                     Source=paste("WGNSSK", ay-1))
  Last <- melt(Last, id.vars = c("Age","Source"))
  Last$Year <- rep(c(ay-2,ay-1,ay), each=7)
  
  This <- data.frame(Age=1:7,
                     'Data year'=sway2[[s]][(nrow(sway2[[s]])-11),],
                     'Intermediate year'=sway2[[s]][(nrow(sway2[[s]])-10),],
                     'Advice year'=sway2[[s]][(nrow(sway2[[s]])-9),],
                     Source=paste("WGNSSK", ay))
  This <- melt(This, id.vars = c("Age","Source"))
  This$Year <- rep(c(ay-1,ay,ay+1), each=7)
  
  stockwts[[s]] <- rbind(Last,This)
  colnames(stockwts[[s]])[c(3,4)] <- c("Type","Weight")
}

write.xlsx(stockwts[[1]],file="output/Forecast_stockwts.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(stockwts[[2]],file="output/Forecast_stockwts.xlsx", sheetName=substocks[2], 
           row.names=F,append=TRUE)
write.xlsx(stockwts[[3]],file="output/Forecast_stockwts.xlsx", sheetName=substocks[3], 
           row.names=F,append=TRUE)

# Forecast_selectivity
Fsel <- list()
for (s in 1:length(fit)){
  sel1 <- do.call(rbind, lapply(prev_FC1$result[[s]], function(x) apply(exp(x$sim[,8:14]), 2, median)))[1:3,]
  sel1 <- t(apply(sel1, 1, function(x) x/max(x)))
  rownames(sel1) <- rownames(attr(prev_FC1$result[[s]],"tab"))[1:3]
  
  sel2 <- do.call(rbind, lapply(FC[[13]]$result[[s]], function(x) apply(exp(x$sim[,8:14]), 2, median)))[1:3,]
  sel2 <- t(apply(sel2, 1, function(x) x/max(x)))
  rownames(sel2) <- rownames(attr(FC[[13]]$result[[s]],"tab"))[1:3]
  
  Last <- data.frame(Age=1:7,
                     'Data year'=sel1[1,],
                     'Intermediate year'=sel1[2,],
                     'Advice year'=sel1[3,],
                     Source=paste("WGNSSK", ay-1))
  Last <- melt(Last, id.vars = c("Age","Source"))
  Last$Year <- rep(c(ay-2,ay-1,ay), each=7)
  
  This <- data.frame(Age=1:7,
                     'Data year'=sel2[1,],
                     'Intermediate year'=sel2[2,],
                     'Advice year'=sel2[3,],
                     Source=paste("WGNSSK", ay))
  This <- melt(This, id.vars = c("Age","Source"))
  This$Year <- rep(c(ay-1,ay,ay+1), each=7)
  
  Fsel[[s]] <- rbind(Last,This)
  colnames(Fsel[[s]])[4] <- "Selectivity"
}

write.xlsx(Fsel[[1]],file="output/Forecast_selectivity.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(Fsel[[2]],file="output/Forecast_selectivity.xlsx", sheetName=substocks[2], 
           row.names=F,append=TRUE)
write.xlsx(Fsel[[3]],file="output/Forecast_selectivity.xlsx", sheetName=substocks[3], 
           row.names=F,append=TRUE)

# N at age compared to previous forecast
Natage <- list()
for (s in 1:length(fit)){
  naa1 <- do.call(rbind, lapply(prev_FC1$result[[s]], function(x) apply(exp(x$sim[,1:7]), 2, median)))[1:3,]
  rownames(naa1) <- rownames(attr(prev_FC1$result[[s]],"tab"))[1:3]
  
  naa2 <- do.call(rbind, lapply(FC[[13]]$result[[s]], function(x) apply(exp(x$sim[,1:7]), 2, median)))[1:3,]
  rownames(naa2) <- rownames(attr(FC[[13]]$result[[s]],"tab"))[1:3]
  
  Last <- data.frame(Age=1:7,
                     'Data year'=naa1[1,],
                     'Intermediate year'=naa1[2,],
                     'Advice year'=naa1[3,],
                     Source=paste("WGNSSK", ay-1))
  Last <- melt(Last, id.vars = c("Age","Source"))
  Last$Year <- rep(c(ay-2,ay-1,ay), each=7)
  
  This <- data.frame(Age=1:7,
                     'Data year'=naa2[1,],
                     'Intermediate year'=naa2[2,],
                     'Advice year'=naa2[3,],
                     Source=paste("WGNSSK", ay))
  This <- melt(This, id.vars = c("Age","Source"))
  This$Year <- rep(c(ay-1,ay,ay+1), each=7)
  
  Natage[[s]] <- rbind(Last,This)
  colnames(Natage[[s]])[c(3,4)] <- c("Type","N")
}

write.xlsx(Natage[[1]],file="output/Forecast_N_at_age.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(Natage[[2]],file="output/Forecast_N_at_age.xlsx", sheetName=substocks[2], 
           row.names=F,append=TRUE)
write.xlsx(Natage[[3]],file="output/Forecast_N_at_age.xlsx", sheetName=substocks[3], 
           row.names=F,append=TRUE)

# compare N at age table 
n_now <- list()
n_prev <- list()
for (s in 1:length(fit)){
  # this year's assessment and forecast
  n_now[[s]] <- as.data.frame(cbind(c(1983:(ay)),ntable(fit, returnList=TRUE)[[s]]))
  colnames(n_now[[s]]) <- c("Year",1:7)
  
  # last year's assessment and forecast
  n_prev[[s]] <- as.data.frame(cbind(c(1983:(ay-1)),ntable(prev_ass, returnList=TRUE)[[s]]))
  colnames(n_prev[[s]]) <- c("Year",1:7)
}

write.xlsx(n_now[[1]],file="output/Now_assessment_N_at_age.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(n_now[[2]],file="output/Now_assessment_N_at_age.xlsx", sheetName=substocks[2],
           row.names=F,append=TRUE)
write.xlsx(n_now[[3]],file="output/Now_assessment_N_at_age.xlsx", sheetName=substocks[3],
           row.names=F,append=TRUE)

write.xlsx(n_prev[[1]],file="output/Prev_assessment_N_at_age.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(n_prev[[2]],file="output/Prev_assessment_N_at_age.xlsx", sheetName=substocks[2], 
           row.names=F,append=TRUE)
write.xlsx(n_prev[[3]],file="output/Prev_assessment_N_at_age.xlsx", sheetName=substocks[3], 
           row.names=F,append=TRUE)

# biomass at age compared to previous forecast 
Batage <- list()
for (s in 1:length(fit)){
  Batage[[s]] <- merge(Natage[[s]],stockwts[[s]])
  Batage[[s]]$B <- Batage[[s]]$N * Batage[[s]]$Weight
}

write.xlsx(Batage[[1]][,c("Age","Source","Type","Year","B")],file="output/Forecast_B_at_age.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(Batage[[2]][,c("Age","Source","Type","Year","B")],file="output/Forecast_B_at_age.xlsx", sheetName=substocks[2], 
           row.names=F,append=TRUE)
write.xlsx(Batage[[3]][,c("Age","Source","Type","Year","B")],file="output/Forecast_B_at_age.xlsx", sheetName=substocks[3], 
           row.names=F,append=TRUE)

# compare B at age table
b_now <- list()
b_prev <- list()
for (s in 1:length(fit)){
  # this year's assessment
  b_now[[s]] <- n_now[[s]][,ac(1:7)]*sway2[[s]][1:length(1983:ay),]
  b_now[[s]]$Year <- n_now[[s]]$Year
  
  # last year's assessment
  b_prev[[s]] <- n_prev[[s]][,ac(1:7)]*sway1[[s]][1:length(1983:(ay-1)),]
  b_prev[[s]]$Year <- n_prev[[s]]$Year
}

write.xlsx(b_now[[1]],file="output/Now_assessment_B_at_age.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(b_now[[2]],file="output/Now_assessment_B_at_age.xlsx", sheetName=substocks[2], 
           row.names=F, append=TRUE)
write.xlsx(b_now[[3]],file="output/Now_assessment_B_at_age.xlsx", sheetName=substocks[3], 
           row.names=F, append=TRUE)

write.xlsx(b_prev[[1]],file="output/Prev_assessment_B_at_age.xlsx", sheetName=substocks[1], 
           row.names=F)
write.xlsx(b_prev[[2]],file="output/Prev_assessment_B_at_age.xlsx", sheetName=substocks[2], 
           row.names=F, append=TRUE)
write.xlsx(b_prev[[3]],file="output/Prev_assessment_B_at_age.xlsx", sheetName=substocks[3], 
           row.names=F, append=TRUE)


# 4 Additional forecasts
load('model/forecast.Rdata')
load("model/refpts/northwestern/NW_1997_type2_refPoints.Rdata")
load("model/refpts/southern/S_1997_type2_5yr_refPoints.Rdata")
load("model/refpts/viking/V_1997_typeM_refPoints.Rdata")
NW <- NW_1997_type2_refPoints
SO <- S_1997_type2_5yr_refPoints
VI <- V_1997_typeM_refPoints

ssbAdviceYear <- lapply(FC[[1]]$result,function(x)attr(x,"shorttab")[3,3])

# checks
# for(i in 1:length(FCR1))print(sapply(FCR1[[i]]$result,function(x)attr(x,"shorttab")[1,2]))
# for(i in 1:length(FCR1))print(sapply(FCR1[[i]]$result,function(x)attr(x,"shorttab")[3,3]))
# for(i in 1:length(FCR1))print(sapply(FCR1[[i]]$result,function(x)attr(x,"shorttab")[2,2]))
# for(i in 1:length(FCR1))print(sapply(FCR1[[i]]$result,function(x)attr(x,"shorttab")[2,3]))
# FCR1[[1]]$result[[1]][[1]]$fbar

Fupper <- c("Northwest" = as.numeric(round(NW['Fmsyupper'],3)), "South" = as.numeric(round(SO['Fmsyupper'],3)), "Viking" = as.numeric(round(VI['Fmsyupper'],3)))
Blim <- c("Northwest" = as.numeric(round(NW['Blim'])), "South" = as.numeric(round(SO['Blim'])), "Viking" = as.numeric(round(VI['Blim'])))
advice <- c("Northwest"=13529, "South"=3922, "Viking"=5240)

AdviceYear = 2025

# Intermediate reductions
scenarios <- as.character(seq(1,1.3,by=0.1))
ns <- length(FCR1)

frange1 <- list()

for (s in 1:3){
  frange1[[s]] <- data.frame(Basis=scenarios,
                             catchAdvice = rep(NA,ns),
                             Landings = rep(NA,ns),
                             Discards = rep(NA,ns),
                             Ftotal = rep(NA,ns),
                             Fland = rep(NA,ns),
                             Fdiscard = rep(NA,ns),
                             ssbAdviceplus1 = rep(NA,ns),
                             "%SSB" = rep(NA,ns),
                             "%advice" = rep(NA,ns),
                             Risk = rep(NA,ns))
  
  names(frange1[[s]])[names(frange1[[s]])=="catchAdvice"] <- paste0("Catch",AdviceYear)
  names(frange1[[s]])[names(frange1[[s]])=="ssbAdviceplus1"] <- paste0("SSB",AdviceYear+1)
  
  for (i in c(1:ns)){
    frange1[[s]][i,2] <- attr(FCR1[[i]]$result[[s]],"tab")[3,10]
    frange1[[s]][i,3] <- attr(FCR1[[i]]$result[[s]],"tab")[3,19]
    frange1[[s]][i,4] <- attr(FCR1[[i]]$result[[s]],"tab")[3,22]
    frange1[[s]][i,5] <- attr(FCR1[[i]]$result[[s]],"tab")[3,1]
    frange1[[s]][i,6] <- attr(FCR1[[i]]$result[[s]],"tab")[3,13]
    frange1[[s]][i,7] <- attr(FCR1[[i]]$result[[s]],"tab")[3,16]
    frange1[[s]][i,8] <- attr(FCR1[[i]]$result[[s]],"tab")[4,7]
    frange1[[s]][i,11] <- lapply(FCR1[[i]]$result[[s]], function(xx) length(which(xx$ssb < Blim[[s]]))/10)[[4]]
  }
  
  frange1[[s]][,9] <- (frange1[[s]][,8]-ssbAdviceYear[[s]])/ssbAdviceYear[[s]]*100
  frange1[[s]][,10] <- (frange1[[s]][,2]-advice[s])/advice[s]*100  
  
  # ICES round F (landings and discards) & %s
  for (col in c(6:7,9:11)){
    frange1[[s]][,col] <- icesRound(frange1[[s]][,col])
  }
}

names(frange1) <- names(FCR1[[1]]$result)

for (s in 1:3){
  write.xlsx(frange1[[s]], file = "output/frange1.xlsx",
             sheetName = names(frange1)[[s]], row.names = FALSE, append = s!=1)
}

# F increments

# checks
# for(i in 1:length(FCR2))print(sapply(FCR2[[i]]$result,function(x)attr(x,"shorttab")[1,2]))
# for(i in 1:length(FCR2))print(sapply(FCR2[[i]]$result,function(x)attr(x,"shorttab")[3,3]))
# for(i in 1:length(FCR2))print(sapply(FCR2[[i]]$result,function(x)attr(x,"shorttab")[2,2]))
# for(i in 1:length(FCR2))print(sapply(FCR2[[i]]$result,function(x)attr(x,"shorttab")[2,3]))
# FCR2[[1]]$result[[1]][[1]]$fbar

scenarios <- as.character(c(0.000001, seq(0.01,round(max(Fupper), digits=2),0.01)))
ns <- length(FCR2)

frange2 <- list()

for (s in 1:3){
  frange2[[s]] <- data.frame(Basis=scenarios,
                             catchAdvice = rep(NA,ns),
                             Landings = rep(NA,ns),
                             Discards = rep(NA,ns),
                             Ftotal = rep(NA,ns),
                             Fland = rep(NA,ns),
                             Fdiscard = rep(NA,ns),
                             ssbAdviceplus1 = rep(NA,ns),
                             "%SSB" = rep(NA,ns),
                             "%advice" = rep(NA,ns),
                             Risk = rep(NA,ns))
  
  names(frange2[[s]])[names(frange2[[s]])=="catchAdvice"] <- paste0("Catch",AdviceYear)
  names(frange2[[s]])[names(frange2[[s]])=="ssbAdviceplus1"] <- paste0("SSB",AdviceYear+1)
  
  for (i in c(1:ns)){
    frange2[[s]][i,2] <- attr(FCR2[[i]]$result[[s]],"tab")[3,10]
    frange2[[s]][i,3] <- attr(FCR2[[i]]$result[[s]],"tab")[3,19]
    frange2[[s]][i,4] <- attr(FCR2[[i]]$result[[s]],"tab")[3,22]
    frange2[[s]][i,5] <- attr(FCR2[[i]]$result[[s]],"tab")[3,1]
    frange2[[s]][i,6] <- attr(FCR2[[i]]$result[[s]],"tab")[3,13]
    frange2[[s]][i,7] <- attr(FCR2[[i]]$result[[s]],"tab")[3,16]
    frange2[[s]][i,8] <- attr(FCR2[[i]]$result[[s]],"tab")[4,7]
    frange2[[s]][i,11] <- lapply(FCR2[[i]]$result[[s]], function(xx) length(which(xx$ssb < Blim[[s]]))/10)[[4]]
  }
  
  frange2[[s]][,9] <- (frange2[[s]][,8]-ssbAdviceYear[[s]])/ssbAdviceYear[[s]]*100
  frange2[[s]][,10] <- (frange2[[s]][,2]-advice[s])/advice[s]*100  
  
  # ICES round F & %s
  for (col in c(5:7,9:11)){
    frange2[[s]][,col] <- icesRound(frange2[[s]][,col])
  }
}

names(frange2) <- names(FCR2[[1]]$result)

for (s in 1:3){
  write.xlsx(frange2[[s]], file = "output/frange.xlsx",
             sheetName = names(frange2)[[s]], row.names = FALSE, append = s!=1)
}

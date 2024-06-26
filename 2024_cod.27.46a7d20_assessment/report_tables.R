library(icesTAF)
library(stockassessment)
library(multiStockassessment)
library(icesAdvice)
library(xlsx)

mkdir("report/tables")

## 1 Assessment
load('data/data.Rdata')
load('model/model.Rdata')

# Configuration
saveConf(fit[[1]]$conf, file="report/tables/conf.cfg")

# Vulnerability matrix
write.table(dat$Shared_data$keyFleetStock,file="report/tables/vulnerability.txt",row.names=TRUE,col.names=TRUE)

# Model fit
write.xlsx(modeltable(fit), file="report/tables/Tab 4.7d fit.xlsx", sheetName="fit")

# F-at-age
faa <- faytable(fit, returnList=TRUE)
for(i in 1:length(fit)){
  faa[[i]] <- as.data.frame(cbind(faa[[i]], fbartable(fit,returnList=TRUE)[[i]][,1])[-nrow(faa[[i]]),])
  colnames(faa[[i]]) <- c(as.character(1:6), "7+", "Fbar 2-4")
  faa[[i]] <- round(faa[[i]], 3)
  write.xlsx(faa[[i]], file="report/tables/Tab 4.8 Faa.xlsx", sheetName = names(faa)[[i]], append=(i!=1))
}

# N-at-age
naa <- ntable(fit, returnList=TRUE)
for(i in 1:length(fit)){
  naa[[i]] <- as.data.frame(naa[[i]])
  naa[[i]]$Total <- rowSums(naa[[i]])
  colnames(naa[[i]]) <- c(as.character(1:6), "7+", "Total")
  naa[[i]] <- round(naa[[i]])
  write.xlsx(naa[[i]], file="report/tables/Tab 4.9 Naa.xlsx", sheetName = names(naa)[[i]], append=(i!=1))
}

# Catch-at-age
rp <- attr(fit, "m_rep")
caa <- list()
for(i in 1:length(fit)){
  caa[[i]] <- as.data.frame(t(apply(exp(rp$logCatchByFleetAge[[i]][,,]),1:2,sum)))
  caa[[i]] <- caa[[i]][-nrow(caa[[i]]),]
  colnames(caa[[i]]) <- c(as.character(1:6),"7+")
  rownames(caa[[i]]) <- 1983:(1983+nrow(caa[[i]])-1)
  caa[[i]] <- round(caa[[i]])
  write.xlsx(caa[[i]], file="report/tables/Tab 4.10 Caa.xlsx", sheetName = c("Northwest","South","Viking")[[i]], append=(i!=1))
}

# Stock weights-at-age
sw <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logSW))
for(i in 1:length(fit)){
  sw[[i]] <- as.data.frame(sw[[i]][1:(dim(sw[[i]])[1]-10),])
  colnames(sw[[i]]) <- c(as.character(1:6),"7+")
  rownames(sw[[i]]) <- 1983:(1983+nrow(sw[[i]])-1)
  sw[[i]] <- round(sw[[i]],3)
  write.xlsx(sw[[i]], file="report/tables/Tab 4.11 SW.xlsx", sheetName = c("Northwest","South","Viking")[[i]], append=(i!=1))
}

# Maturity-at-age
mo <- multiStockassessment:::splitParameter(plogis(attr(fit,"m_pl")$logitMO))
for(i in 1:length(fit)){
  mo[[i]] <- as.data.frame(mo[[i]][1:(dim(mo[[i]])[1]-10),])
  colnames(mo[[i]]) <- c(as.character(1:6),"7+")
  rownames(mo[[i]]) <- 1983:(1983+nrow(mo[[i]])-1)
  mo[[i]] <- round(mo[[i]],3)
  write.xlsx(mo[[i]], file="report/tables/Tab 4.12 MO.xlsx", sheetName = c("Northwest","South","Viking")[[i]], append=(i!=1))
}

# Natural mortality-at-age
nm <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logNM))
for(i in 1:length(fit)){
  nm[[i]] <- as.data.frame(nm[[i]][1:(dim(nm[[i]])[1]-10),])
  colnames(nm[[i]]) <- c(as.character(1:6),"7+")
  rownames(nm[[i]]) <- 1983:(1983+nrow(nm[[i]])-1)
  nm[[i]] <- round(nm[[i]],3)
  write.xlsx(nm[[i]], file="report/tables/Tab 4.13 NM.xlsx", sheetName = c("Northwest","South","Viking")[[i]], append=(i!=1))
}

# Summary table
summ <- summary(fit, returnList=TRUE)
for(i in 1:length(fit)){
  summ[[i]] <- as.data.frame(cbind(summ[[i]],tsbtable(fit, returnList=TRUE)[[i]], rbind(catchtable(fit, returnList=TRUE)[[i]], rep(NA,3))))
  colnames(summ[[i]])[c(10,13)] <- c("TSB","Catches")
  summ[[i]][,10:15] <- round(summ[[i]][,10:15])
  summ[[i]] <- summ[[i]][,c(1:3,10:12,4:6,13:15,7:9)]
  summ[[i]][nrow(summ[[i]]), 13:15] <- NA
  write.xlsx(summ[[i]], file="report/tables/Tab 4.14 Summary.xlsx", sheetName = names(summ)[[i]], append=(i!=1))
}

## 2 Forecasts
load('model/forecast.Rdata')
load("model/refpts/northwestern/NW_1997_type2_refPoints.Rdata")
load("model/refpts/southern/S_1997_type2_5yr_refPoints.Rdata")
load("model/refpts/viking/V_1997_typeM_refPoints.Rdata")
NW <- NW_1997_type2_refPoints
SO <- S_1997_type2_5yr_refPoints
VI <- V_1997_typeM_refPoints

ssbAdviceYear <- lapply(FC[[1]]$result,function(x)attr(x,"shorttab")[3,3])

# checks
# for(i in 1:length(FC))print(sapply(FC[[i]]$result,function(x)attr(x,"shorttab")[1,2]))
# for(i in 1:length(FC))print(sapply(FC[[i]]$result,function(x)attr(x,"shorttab")[3,3]))
# for(i in 1:length(FC))print(sapply(FC[[i]]$result,function(x)attr(x,"shorttab")[2,2]))
# for(i in 1:length(FC))print(sapply(FC[[i]]$result,function(x)attr(x,"shorttab")[2,3]))
# FC[[1]]$result[[1]][[1]]$fbar

Blim <- c("Northwest" = as.numeric(round(NW['Blim'])), "South" = as.numeric(round(SO['Blim'])), "Viking" = as.numeric(round(VI['Blim'])))
advice <- c("Northwest"=13529, "South"=3922, "Viking"=5240)  

AdviceYear = 2025

scenarios <- sapply(FC,function(x)x$name)
ns <- length(FC)

cs <- list()
ca <- list()

for (s in 1:3){
  cs[[s]] <- data.frame(Basis=scenarios,
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
  
  names(cs[[s]])[names(cs[[s]])=="catchAdvice"] <- paste0("Catch",AdviceYear)
  names(cs[[s]])[names(cs[[s]])=="ssbAdviceplus1"] <- paste0("SSB",AdviceYear+1)
  
  for (i in c(1:ns)){
    cs[[s]][i,2] <- attr(FC[[i]]$result[[s]],"tab")[3,10]
    cs[[s]][i,3] <- attr(FC[[i]]$result[[s]],"tab")[3,19]
    cs[[s]][i,4] <- attr(FC[[i]]$result[[s]],"tab")[3,22]
    cs[[s]][i,5] <- attr(FC[[i]]$result[[s]],"tab")[3,1]
    cs[[s]][i,6] <- attr(FC[[i]]$result[[s]],"tab")[3,13]
    cs[[s]][i,7] <- attr(FC[[i]]$result[[s]],"tab")[3,16]
    cs[[s]][i,8] <- attr(FC[[i]]$result[[s]],"tab")[4,7]
    cs[[s]][i,11] <- lapply(FC[[i]]$result[[s]], function(xx) length(which(xx$ssb < Blim[[s]]))/10)[[4]]
  }
  
  cs[[s]][,9] <- (cs[[s]][,8]-ssbAdviceYear[[s]])/ssbAdviceYear[[s]]*100
  cs[[s]][,10] <- (cs[[s]][,2]-advice[s])/advice[s]*100
  
  assumptions <- c(
    paste0("Fbar(", AdviceYear-1,")"),
    paste0("SSB(", AdviceYear,")"),
    paste0("R(", AdviceYear-1,")"),
    paste0("R(", AdviceYear,")"),
    paste0("Catch(", AdviceYear-1,")"),
    paste0("Landings(", AdviceYear-1,")"),
    paste0("Discards(", AdviceYear-1,")"),
    #paste0("TAC(", AdviceYear-1,")"),
    paste0("Advice for ", AdviceYear-1))
  
  ca[[s]] <- data.frame(variable = assumptions, value = rep(NA,8))
  ca[[s]][1,2] <- attr(FC[[1]]$result[[s]],"tab")[2,1]
  ca[[s]][2,2] <- attr(FC[[1]]$result[[s]],"tab")[3,7]
  ca[[s]][3,2] <- attr(FC[[1]]$result[[s]],"tab")[2,4]
  ca[[s]][4,2] <- attr(FC[[1]]$result[[s]],"tab")[3,4]
  ca[[s]][5,2] <- attr(FC[[1]]$result[[s]],"tab")[2,10]
  ca[[s]][6,2] <- attr(FC[[1]]$result[[s]],"tab")[2,19]
  ca[[s]][7,2] <- attr(FC[[1]]$result[[s]],"tab")[2,22]
  #ca[[s]][8,2] <- TAC
  ca[[s]][8,2] <- advice[s]
  
  # ICES round F (landings and discards) & %s
  for (col in c(6:7,9:11)){
    cs[[s]][,col] <- icesRound(cs[[s]][,col])
  }
  
  # ICES round non-RP Fs.
  # Check for zeros as the last decimal!
  cs[[s]][c(4,7:10,13),5] <- as.numeric(icesRound(cs[[s]][c(4,7:10,13),5]))
  
}

names(ca) <- names(cs) <- names(FC[[1]]$result)

for (s in 1:3){
  write.xlsx(ca[[s]], file = "report/tables/Tab 4.15 forecasts.xlsx",
             sheetName = paste0(names(ca)[[s]],"_assumptions"), row.names = FALSE, append = s!=1)
  write.xlsx(cs[[s]], file = "report/tables/Tab 4.15 forecasts.xlsx",
             sheetName = paste0(names(ca)[[s]],"_scenarios"), row.names = FALSE, append = TRUE)
}

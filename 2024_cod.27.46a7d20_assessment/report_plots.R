## Prepare plots for report

## Before: Assessment input files, data.Rdata, model.Rdata, residuals.Rdata, retro.Rdata, forecast.Rdata
## After: Report plots

library(icesTAF)
library(ggplot2)
library(reshape2)
library(stockassessment)
library(multiStockassessment)

source("utilities.r")
source("utilities_multiStockPlots.R")

mkdir("report/plots")


## 1 Assessment
load('model/model.Rdata')
load('model/residuals.Rdata')
load('model/retro.Rdata')
load('model/forecast.Rdata')

# SSB
png("report/plots/Fig 4.9a SSB.png", width=952, height=609)
makePlot(fit, "ssb", cols=c("blue","green","red","black"))
dev.off()

# Fbar
png("report/plots/Fig 4.9b Fbar.png", width=952, height=609)
makePlot(fit, "fbar", cols=c("blue","green","red","black"))
dev.off()

# Recruitment
png("report/plots/Fig 4.9c Rec.png", width=952, height=609)
makePlot(fit, "rec", cols=c("blue","green","red","black"))
dev.off()

# Catch
png("report/plots/Fig 4.9d Catch.png", width=952, height=609)
makePlot(fit, "catch", cols=c("blue","green","red","black"))
v <- observedCatch(fit)
points(as.numeric(names(v)),v, col = "black",pch="x",cex=1.5)
dev.off()

# TSB
png("report/plots/Fig 4.9e TSB.png", width=952, height=609)
makePlot(fit, "tsb", cols=c("blue","green","red","black"))
dev.off()

# Difference between observed and predicted catches
# Not a plot but the best place to extract this information
cdif <- data.frame(Observed=v,
                   Estimated=catchtable(fit, addTotal=TRUE,returnList = TRUE)$Total[-42,1])
cdif$diff <- cdif$Observed - cdif$Estimated
cdif$perc <- (cdif$Observed / cdif$Estimated -1) * 100

# F-at-age
fay <- faytable(fit,returnList = TRUE)
png("report/plots/Fig 4.10a Faa.png", width=952, height=609)
par(mfrow=c(1,3))
matplot(1983:(1983+nrow(fay[[1]])-2), fay[[1]][-nrow(fay[[1]]),], type="l", lty="solid", col=1:ncol(fay[[1]]), lwd=2, xlab="Year", ylab="F-at-age", main="Northwest", ylim=c(0,1.5))
matplot(1983:(1983+nrow(fay[[2]])-2), fay[[2]][-nrow(fay[[2]]),], type="l", lty="solid", col=1:ncol(fay[[2]]), lwd=2, xlab="Year", ylab="F-at-age", main="South", ylim=c(0,1.5))
matplot(1983:(1983+nrow(fay[[3]])-2), fay[[3]][-nrow(fay[[3]]),], type="l", lty="solid", col=1:ncol(fay[[3]]), lwd=2, xlab="Year", ylab="F-at-age", main="Viking", ylim=c(0,1.5))
legend('topright', col=1:ncol(fay[[3]]), legend=c(colnames(fay[[3]])[1:6],"Age 7+"), bty='n', lwd=2)
dev.off()

# Selectivity in F
png("report/plots/Fig 4.10b Fsel.png", width=952, height=609)
par(mfrow=c(1,3))
for (s in 1:3){
  # From the forecasts
  sel <- do.call(rbind, lapply(FC[[13]]$result[[s]], function(x) apply(exp(x$sim[,8:14]), 2, median)))
  sel <- t(apply(sel, 1, function(x) x/max(x)))
  rownames(sel) <- rownames(attr(FC[[1]]$result[[s]],"tab"))

  # And from the assessment
  fay[[s]] <- fay[[s]][-c((nrow(fay[[s]])-1),nrow(fay[[s]])),]
  fay[[s]] <- rbind(fay[[s]], sel)

  fmat <- fay[[s]]
  barplot(t(fmat/rowSums(fmat)), border = NA, space = c(0),
          xlab = "Year", main = "Selectivity in F")
  text(1, cumsum(t(fmat/rowSums(fmat))[, 1]) - 0.5 * t(fmat/rowSums(fmat))[, 1], label = as.character(1:ncol(fmat)),
       adj = c(0, 0.2))
  abline(v=nrow(fmat)-3, lty=2)
}
dev.off()

# N-at-age
nay <- ntable(fit,returnList = TRUE)
png("report/plots/Fig 4.11 Naa.png", width=952, height=609)
par(mfrow=c(1,3))
matplot(1983:(1983+nrow(nay[[1]])-1), log(nay[[1]]), type="l", lty="solid", col=1:ncol(nay[[1]]), lwd=2, xlab="Year", ylab="log N-at-age", main="Northwest", ylim=c(2,14))
matplot(1983:(1983+nrow(nay[[2]])-1), log(nay[[2]]), type="l", lty="solid", col=1:ncol(nay[[2]]), lwd=2, xlab="Year", ylab="log N-at-age", main="South", ylim=c(2,14))
matplot(1983:(1983+nrow(nay[[3]])-1), log(nay[[3]]), type="l", lty="solid", col=1:ncol(nay[[3]]), lwd=2, xlab="Year", ylab="log N-at-age", main="Viking", ylim=c(2,14))
legend('topright', col=1:ncol(nay[[3]]), legend=c(colnames(nay[[3]])[1:6],"Age 7+"), bty='n', lwd=2)
dev.off()

# Catch-at-age
rp <- attr(fit, "m_rep")
cay <- list()
for (i in 1:length(fit)){
  cay[[i]] <- t(apply(exp(rp$logCatchByFleetAge[[i]][,,]),1:2,sum))
}
png("report/plots/Fig 4.12 Caa.png", width=952, height=609)
par(mfrow=c(1,3))
matplot(1983:(1983+nrow(cay[[1]])-2), log(cay[[1]][-nrow(cay[[1]]),]), type="l", lty="solid", col=1:ncol(cay[[1]]), lwd=2, xlab="Year", ylab="log Catch-at-age", main="Northwest", ylim=c(1,12))
matplot(1983:(1983+nrow(cay[[2]])-2), log(cay[[2]][-nrow(cay[[1]]),]), type="l", lty="solid", col=1:ncol(cay[[2]]), lwd=2, xlab="Year", ylab="log Catch-at-age", main="South", ylim=c(1,12))
matplot(1983:(1983+nrow(cay[[3]])-2), log(cay[[3]][-nrow(cay[[1]]),]), type="l", lty="solid", col=1:ncol(cay[[3]]), lwd=2, xlab="Year", ylab="log Catch-at-age", main="Viking", ylim=c(1,12))
legend('topright', col=1:ncol(cay[[3]]), legend=c(paste("Age",1:6),"Age 7+"), bty='n', lwd=2)
dev.off()

# Stock weights
sway <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logSW))
sway_lo <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logSW-2*attr(fit,"m_plsd")$logSW))
sway_up <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logSW+2*attr(fit,"m_plsd")$logSW))
png("report/plots/Fig 4.13 SW.png", width=952, height=609)
par(mfrow=c(1,3))
matplot(1983:(1983+nrow(sway[[1]])-9), sway[[1]][1:(nrow(sway[[1]])-8),], type="l", lty="solid", col=1:ncol(sway[[1]]), lwd=2, xlab="Year", ylab="Stock weights", main="Northwest", ylim=c(0,16))
matplot(1983:(1983+nrow(fit[[1]]$data$stockMeanWeight)-1), fit[[1]]$data$stockMeanWeight, col=1:ncol(sway[[1]]), add=TRUE)
matplot(1983:(1983+nrow(sway_lo[[1]])-9), sway_lo[[1]][1:(nrow(sway_lo[[1]])-8),], type="l", lty="dotted", col=1:ncol(sway[[1]]), add=TRUE)
matplot(1983:(1983+nrow(sway_up[[1]])-9), sway_up[[1]][1:(nrow(sway_up[[1]])-8),], type="l", lty="dotted", col=1:ncol(sway[[1]]), add=TRUE)
matplot(1983:(1983+nrow(sway[[2]])-9), sway[[2]][1:(nrow(sway[[1]])-8),], type="l", lty="solid", col=1:ncol(sway[[2]]), lwd=2, xlab="Year", ylab="Stock weights", main="South", ylim=c(0,16))
matplot(1983:(1983+nrow(fit[[2]]$data$stockMeanWeight)-1), fit[[2]]$data$stockMeanWeight, col=1:ncol(sway[[2]]), add=TRUE)
matplot(1983:(1983+nrow(sway_lo[[2]])-9), sway_lo[[2]][1:(nrow(sway_lo[[2]])-8),], type="l", lty="dotted", col=1:ncol(sway[[2]]), add=TRUE)
matplot(1983:(1983+nrow(sway_up[[2]])-9), sway_up[[2]][1:(nrow(sway_up[[2]])-8),], type="l", lty="dotted", col=1:ncol(sway[[2]]), add=TRUE)
matplot(1983:(1983+nrow(sway[[3]])-9), sway[[3]][1:(nrow(sway[[1]])-8),], type="l", lty="solid", col=1:ncol(sway[[3]]), lwd=2, xlab="Year", ylab="Stock weights", main="Viking", ylim=c(0,16))
matplot(1983:(1983+nrow(fit[[3]]$data$stockMeanWeight)-1), fit[[3]]$data$stockMeanWeight, col=1:ncol(sway[[3]]), add=TRUE)
matplot(1983:(1983+nrow(sway_lo[[3]])-9), sway_lo[[3]][1:(nrow(sway_lo[[3]])-8),], type="l", lty="dotted", col=1:ncol(sway[[3]]), add=TRUE)
matplot(1983:(1983+nrow(sway_up[[3]])-9), sway_up[[3]][1:(nrow(sway_up[[3]])-8),], type="l", lty="dotted", col=1:ncol(sway[[3]]), add=TRUE)
legend('topright', col=1:ncol(sway[[3]]), legend=paste("Age",c(1:6,"7+")), bty='n', lwd=2)
dev.off()

## Maturity
mo <- multiStockassessment:::splitParameter(plogis(attr(fit,"m_pl")$logitMO))
mo_lo <- multiStockassessment:::splitParameter(plogis(attr(fit,"m_pl")$logitMO-2*attr(fit,"m_plsd")$logitMO))
mo_up <- multiStockassessment:::splitParameter(plogis(attr(fit,"m_pl")$logitMO+2*attr(fit,"m_plsd")$logitMO))
png("report/plots/Fig 4.14 MO.png", width=952, height=609)
par(mfrow=c(1,3))
matplot(1983:(1983+nrow(mo[[1]])-9), mo[[1]][1:(nrow(mo[[1]])-8),], type="l", lty="solid", col=1:ncol(mo[[1]]), lwd=2, xlab="Year", ylab="Maturity", main="Northwest", ylim=c(0,1.2))
matplot(1983:(1983+nrow(fit[[1]]$data$propMat)-1), fit[[1]]$data$propMat, col=1:ncol(mo[[1]]), add=TRUE)
matplot(1983:(1983+nrow(mo_lo[[1]])-9), mo_lo[[1]][1:(nrow(mo_lo[[1]])-8),], type="l", lty="dotted", col=1:ncol(mo[[1]]), add=TRUE)
matplot(1983:(1983+nrow(mo_up[[1]])-9), mo_up[[1]][1:(nrow(mo_up[[1]])-8),], type="l", lty="dotted", col=1:ncol(mo[[1]]), add=TRUE)
matplot(1983:(1983+nrow(mo[[2]])-9), mo[[2]][1:(nrow(mo[[1]])-8),], type="l", lty="solid", col=1:ncol(mo[[2]]), lwd=2, xlab="Year", ylab="Maturity", main="South", ylim=c(0,1.2))
matplot(1983:(1983+nrow(fit[[2]]$data$propMat)-1), fit[[2]]$data$propMat, col=1:ncol(mo[[2]]), add=TRUE)
matplot(1983:(1983+nrow(mo_lo[[2]])-9), mo_lo[[2]][1:(nrow(mo_lo[[2]])-8),], type="l", lty="dotted", col=1:ncol(mo[[2]]), add=TRUE)
matplot(1983:(1983+nrow(mo_up[[2]])-9), mo_up[[2]][1:(nrow(mo_up[[2]])-8),], type="l", lty="dotted", col=1:ncol(mo[[2]]), add=TRUE)
matplot(1983:(1983+nrow(mo[[3]])-9), mo[[3]][1:(nrow(mo[[1]])-8),], type="l", lty="solid", col=1:ncol(mo[[3]]), lwd=2, xlab="Year", ylab="Maturity", main="Viking", ylim=c(0,1.2))
matplot(1983:(1983+nrow(fit[[3]]$data$propMat)-1), fit[[3]]$data$propMat, col=1:ncol(mo[[3]]), add=TRUE)
matplot(1983:(1983+nrow(mo_lo[[3]])-9), mo_lo[[3]][1:(nrow(mo_lo[[3]])-8),], type="l", lty="dotted", col=1:ncol(mo[[3]]), add=TRUE)
matplot(1983:(1983+nrow(mo_up[[3]])-9), mo_up[[3]][1:(nrow(mo_up[[3]])-8),], type="l", lty="dotted", col=1:ncol(mo[[3]]), add=TRUE)
legend('topright', col=1:ncol(mo[[3]]), legend=paste("Age",c(1:6,"7+")), bty='n', lwd=2)
dev.off()

# Natural mortality
nm <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logNM))
nm_lo <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logNM-2*attr(fit,"m_plsd")$logNM))
nm_up <- multiStockassessment:::splitParameter(exp(attr(fit,"m_pl")$logNM+2*attr(fit,"m_plsd")$logNM))
png("report/plots/Fig 4.15 NM.png", width=952, height=609)
par(mfrow=c(1,3))
matplot(1983:(1983+nrow(nm[[1]])-9), nm[[1]][1:(nrow(nm[[1]])-8),1:4], type="l", lty="solid", col=1:ncol(nm[[1]]), lwd=2, xlab="Year", ylab="Natural mortality", main="Northwest", ylim=c(0,1.6))
matplot(1983:(1983+nrow(fit[[1]]$data$natMor)-1), fit[[1]]$data$natMor[,1:4], col=1:ncol(nm[[1]]), add=TRUE)
matplot(1983:(1983+nrow(nm_lo[[1]])-9), nm_lo[[1]][1:(nrow(nm_lo[[1]])-8),1:4], type="l", lty="dotted", col=1:ncol(nm[[1]]), add=TRUE)
matplot(1983:(1983+nrow(nm_up[[1]])-9), nm_up[[1]][1:(nrow(nm_up[[1]])-8),1:4], type="l", lty="dotted", col=1:ncol(nm[[1]]), add=TRUE)
matplot(1983:(1983+nrow(nm[[2]])-9), nm[[2]][1:(nrow(nm[[1]])-8),1:4], type="l", lty="solid", col=1:ncol(nm[[2]]), lwd=2, xlab="Year", ylab="Natural mortality", main="South", ylim=c(0,1.6))
matplot(1983:(1983+nrow(fit[[2]]$data$natMor)-1), fit[[2]]$data$natMor[,1:4], col=1:ncol(nm[[2]]), add=TRUE)
matplot(1983:(1983+nrow(nm_lo[[2]])-9), nm_lo[[2]][1:(nrow(nm_lo[[2]])-8),1:4], type="l", lty="dotted", col=1:ncol(nm[[2]]), add=TRUE)
matplot(1983:(1983+nrow(nm_up[[2]])-9), nm_up[[2]][1:(nrow(nm_up[[2]])-8),1:4], type="l", lty="dotted", col=1:ncol(nm[[2]]), add=TRUE)
matplot(1983:(1983+nrow(nm[[3]])-9), nm[[3]][1:(nrow(nm[[1]])-8),1:4], type="l", lty="solid", col=1:ncol(nm[[3]]), lwd=2, xlab="Year", ylab="Natural mortality", main="Viking", ylim=c(0,1.6))
matplot(1983:(1983+nrow(fit[[3]]$data$natMor)-1), fit[[3]]$data$natMor[,1:4], col=1:ncol(nm[[3]]), add=TRUE)
matplot(1983:(1983+nrow(nm_lo[[3]])-9), nm_lo[[3]][1:(nrow(nm_lo[[3]])-8),1:4], type="l", lty="dotted", col=1:ncol(nm[[3]]), add=TRUE)
matplot(1983:(1983+nrow(nm_up[[3]])-9), nm_up[[3]][1:(nrow(nm_up[[3]])-8),1:4], type="l", lty="dotted", col=1:ncol(nm[[3]]), add=TRUE)
legend('topright', col=1:4, legend=paste("Age",c(1:4)), bty='n', lwd=2)
dev.off()

# Residuals
png("report/plots/Fig 4.16 RES.png", width=1228, height=787)
plot(RES)
dev.off()

# Retros
makeRetroPlot(RETRO, "ssb",1:4, cols = c("blue","green","red","black"))
png("report/plots/Fig 4.17a SSBretro.png", width=721,height=725)
par(mfrow=c(3,1))
makeRetroPlot(RETRO, "ssb",1, cols = c("blue"))
makeRetroPlot(RETRO, "ssb",2, cols = c("green"))
makeRetroPlot(RETRO, "ssb",3, cols = c("red"))
dev.off()

makeRetroPlot(RETRO, "fbar",1:4, cols = c("blue","green","red","black"))
png("report/plots/Fig 4.17b Fretro.png", width=721,height=725)
par(mfrow=c(3,1))
makeRetroPlot(RETRO, "fbar",1, cols = c("blue"))
makeRetroPlot(RETRO, "fbar",2, cols = c("green"))
makeRetroPlot(RETRO, "fbar",3, cols = c("red"))
dev.off()

makeRetroPlot(RETRO, "rec",1:4, cols = c("blue","green","red","black"))
png("report/plots/Fig 4.17c Rretro.png", width=721,height=725)
par(mfrow=c(3,1))
makeRetroPlot(RETRO, "rec",1, cols = c("blue"))
makeRetroPlot(RETRO, "rec",2, cols = c("green"))
makeRetroPlot(RETRO, "rec",3, cols = c("red"))
dev.off()

makeRetroPlot(RETRO, "tsb",1:4, cols = c("blue","green","red","black"))
png("report/plots/Fig 4.17d TSBretro.png", width=721,height=725)
par(mfrow=c(3,1))
makeRetroPlot(RETRO, "tsb",1, cols = c("blue"))
makeRetroPlot(RETRO, "tsb",2, cols = c("green"))
makeRetroPlot(RETRO, "tsb",3, cols = c("red"))
dev.off()

# Assessment comparisons
WG <- fit
load("../../WGNSSK 2023/TAF2023/model/model.Rdata")
WG23 <- fit
# load("model/model20.Rdata")
# old_M <- fit
# load("model/model20fix.Rdata")
# old_Mfix <- fit

# Northwest
ssb <- data.frame(Year=as.numeric(row.names(ssbtable(WG,returnList=TRUE)[[1]])),
                  WGNSSK2023=c(ssbtable(WG23,returnList=TRUE)[[1]][,1],NA),
                  WGNSSK2024=ssbtable(WG,returnList=TRUE)[[1]][,1])#,
                  # key20=ssbtable(old_M,returnList=TRUE)[[1]][,1],
                  # key20fix=ssbtable(old_Mfix,returnList=TRUE)[[1]][,1])
ssb <- melt(ssb, id = "Year")
ssb$plot <- "SSB"
fbar <- data.frame(Year=as.numeric(row.names(fbartable(WG,returnList=TRUE)[[1]])),
                  WGNSSK2023=c(fbartable(WG23,returnList=TRUE)[[1]][,1],NA),
                  WGNSSK2024=fbartable(WG,returnList=TRUE)[[1]][,1])#,
                  # key20=fbartable(old_M,returnList=TRUE)[[1]][,1],
                  # key20fix=fbartable(old_Mfix,returnList=TRUE)[[1]][,1])
rownames(fbar) <- fbar$Year
fbar["2023","WGNSSK2023"] <- fbar["2023","WGNSSK2024"] <- NA #<- fbar["2023","key20"] <- fbar["2023","key20fix"] <- NA 
fbar <- melt(fbar, id = "Year")
fbar$plot <- "Fbar"
rec <- data.frame(Year=as.numeric(row.names(rectable(WG,returnList=TRUE)[[1]])),
                  WGNSSK2023=c(rectable(WG23,returnList=TRUE)[[1]][,1],NA),
                  WGNSSK2024=rectable(WG,returnList=TRUE)[[1]][,1])#,
                  # key20=rectable(old_M,returnList=TRUE)[[1]][,1],
                  # key20fix=rectable(old_Mfix,returnList=TRUE)[[1]][,1])
rec <- melt(rec, id = "Year")
rec$plot <- "Recruitment"
plot <- rbind(ssb,fbar,rec)
plot$plot=factor(plot$plot, levels=c("SSB","Fbar","Recruitment"))
ggplot(plot, aes(x=Year,y=value)) + geom_line(aes(col=variable),lwd=1) +
  facet_grid(plot~., scales="free_y") +
  theme_bw() +
  ylab("")
ggsave("report/plots/Fig 4.18a NWcomp.png")

# South
ssb <- data.frame(Year=as.numeric(row.names(ssbtable(WG,returnList=TRUE)[[2]])),
                  WGNSSK2023=c(ssbtable(WG23,returnList=TRUE)[[2]][,1],NA),
                  WGNSSK2024=ssbtable(WG,returnList=TRUE)[[2]][,1])#,
                  # key20=ssbtable(old_M,returnList=TRUE)[[2]][,1],
                  # key20fix=ssbtable(old_Mfix,returnList=TRUE)[[2]][,1])
ssb <- melt(ssb, id = "Year")
ssb$plot <- "SSB"
fbar <- data.frame(Year=as.numeric(row.names(fbartable(WG,returnList=TRUE)[[2]])),
                   WGNSSK2023=c(fbartable(WG23,returnList=TRUE)[[2]][,1],NA),
                   WGNSSK2024=fbartable(WG,returnList=TRUE)[[2]][,1])#,
                   # key20=fbartable(old_M,returnList=TRUE)[[2]][,1],
                   # key20fix=fbartable(old_Mfix,returnList=TRUE)[[2]][,1])
rownames(fbar) <- fbar$Year
fbar["2023","WGNSSK2023"] <- fbar["2023","WGNSSK2024"] <- NA # <- fbar["2023","key20"] <- fbar["2023","key20fix"] <- NA 
fbar <- melt(fbar, id = "Year")
fbar$plot <- "Fbar"
rec <- data.frame(Year=as.numeric(row.names(rectable(WG,returnList=TRUE)[[2]])),
                  WGNSSK2023=c(rectable(WG23,returnList=TRUE)[[2]][,1],NA),
                  WGNSSK2024=rectable(WG,returnList=TRUE)[[2]][,1])#,
                  # key20=rectable(old_M,returnList=TRUE)[[2]][,1],
                  # key20fix=rectable(old_Mfix,returnList=TRUE)[[2]][,1])
rec <- melt(rec, id = "Year")
rec$plot <- "Recruitment"
plot <- rbind(ssb,fbar,rec)
plot$plot=factor(plot$plot, levels=c("SSB","Fbar","Recruitment"))
ggplot(plot, aes(x=Year,y=value)) + geom_line(aes(col=variable),lwd=1) +
  facet_grid(plot~., scales="free_y") +
  theme_bw() +
  ylab("")
ggsave("report/plots/Fig 4.18b Scomp.png")

# Viking
ssb <- data.frame(Year=as.numeric(row.names(ssbtable(WG,returnList=TRUE)[[3]])),
                  WGNSSK2023=c(ssbtable(WG23,returnList=TRUE)[[3]][,1],NA),
                  WGNSSK2024=ssbtable(WG,returnList=TRUE)[[3]][,1])#,
                  # key20=ssbtable(old_M,returnList=TRUE)[[3]][,1],
                  # key20fix=ssbtable(old_Mfix,returnList=TRUE)[[3]][,1])
ssb <- melt(ssb, id = "Year")
ssb$plot <- "SSB"
fbar <- data.frame(Year=as.numeric(row.names(fbartable(WG,returnList=TRUE)[[3]])),
                   WGNSSK2023=c(fbartable(WG23,returnList=TRUE)[[3]][,1],NA),
                   WGNSSK2024=fbartable(WG,returnList=TRUE)[[3]][,1])#,
                   # key20=fbartable(old_M,returnList=TRUE)[[3]][,1],
                   # key20fix=fbartable(old_Mfix,returnList=TRUE)[[3]][,1])
rownames(fbar) <- fbar$Year
fbar["2023","WGNSSK2023"] <- fbar["2023","WGNSSK2024"] <- NA # <- fbar["2023","key20"] <- fbar["2023","key20fix"] <- NA 
fbar <- melt(fbar, id = "Year")
fbar$plot <- "Fbar"
rec <- data.frame(Year=as.numeric(row.names(rectable(WG,returnList=TRUE)[[3]])),
                  WGNSSK2023=c(rectable(WG23,returnList=TRUE)[[3]][,1],NA),
                  WGNSSK2024=rectable(WG,returnList=TRUE)[[3]][,1]) #,
                  # key20=rectable(old_M,returnList=TRUE)[[3]][,1],
                  # key20fix=rectable(old_Mfix,returnList=TRUE)[[3]][,1])
rec <- melt(rec, id = "Year")
rec$plot <- "Recruitment"
plot <- rbind(ssb,fbar,rec)
plot$plot=factor(plot$plot, levels=c("SSB","Fbar","Recruitment"))
ggplot(plot, aes(x=Year,y=value)) + geom_line(aes(col=variable),lwd=1) +
  facet_grid(plot~., scales="free_y") +
  theme_bw() +
  ylab("")
ggsave("report/plots/Fig 4.18c Vcomp.png")


## 2 Forecasts
load("model/model.Rdata")
load("model/forecast.Rdata")
load("model/refpts/northwestern/NW_1997_type2_refPoints.Rdata")
load("model/refpts/southern/S_1997_type2_5yr_refPoints.Rdata")
load("model/refpts/viking/V_1997_typeM_refPoints.Rdata")
NW <- NW_1997_type2_refPoints
SO <- S_1997_type2_5yr_refPoints
VI <- V_1997_typeM_refPoints

MSYBtrigger <- c("Northwest" = as.numeric(round(NW['MSYbtrigger'])), "South" = as.numeric(round(SO['MSYbtrigger'])), "Viking" = as.numeric(round(VI['MSYbtrigger'])))
Blim <- c("Northwest" = as.numeric(round(NW['Blim'])), "South" = as.numeric(round(SO['Blim'])), "Viking" = as.numeric(round(VI['Blim'])))
Fmsy <- c("Northwest" = as.numeric(round(NW['Fmsy'],3)), "South" = as.numeric(round(SO['Fmsy'],3)), "Viking" = as.numeric(round(VI['Fmsy'],3)))

png("report/plots/Fig 4.19a FCssb.png", width=955,height=615)
makeFCPanel(fit, list(FC[[13]]), "ssb")
dev.off()

png("report/plots/Fig 4.19b FCfbar.png", width=955,height=615)
makeFCPanel(fit, list(FC[[13]]), "fbar")
dev.off()

png("report/plots/Fig 4.19c FCrec.png", width=955,height=615)
makeFCPanel(fit, list(FC[[13]]), "rec")
dev.off()

png("report/plots/Fig 4.19d FCcatch.png", width=955,height=615)
makeFCPanel(fit, list(FC[[13]]), "catch")
dev.off()

## 3 Additional

## Landings weight proportions
# Q1
rp <- attr(fit,"m_rep")
shD <- attr(fit,"m_data")$sharedObs

auxL <- split(as.data.frame(shD$aux),shD$aux[,2])
auxDL <- split(as.data.frame(shD$auxData),shD$aux[,2])
predL <- split(as.data.frame(rp$predPerStock),shD$aux[,2])

e <- new.env()
load(file.path("boot/data/WGNSSK2024_SubstockProportions.RData"),e)

PredNew <- exp(predL[[10]])
ObsNew <- xtabs(prop ~ year + areaPopulation, e$substockProps)

yy2 <- as.numeric(rownames(ObsNew))

png("report/plots/LPfitQ1.png", width=762,height=571)
matplot(yy2,ObsNew,col=c('blue','green','red'),type="p",pch=16,cex=2)
matplot(yy2,PredNew[seq_len(nrow(PredNew))%%2==1,],type="l",add=TRUE, col=c('blue','green','red'),lty=1,lwd=3)
legend("top",c("Northwest","South","Viking"),col=c('blue',"green","red"),lwd=3,bty="n",ncol=3)
dev.off()

## Whole year
maxAge <- 7
maxYear <- Inf
minYear <- 1983
ln_WS <- read.ices("boot/data/lf_WS.dat") * read.ices("boot/data/cn_WS.dat")
ln_NS <- read.ices("boot/data/lf_NS.dat")
ln <- cutSum(ln_NS) + cutSum(ln_WS)
lw_WS <- read.ices("boot/data/lw_WS.dat")
lw_NS <- read.ices("boot/data/lw_NS.dat")
lw <- (na2zero(cutMean(lw_NS,ln_NS) * cutSum(ln_NS)) + na2zero(cutMean(lw_NS,ln_NS) * cutSum(ln_WS))) / ln
cn_WS <- read.ices("boot/data/cn_WS.dat")
cn_NS <- read.ices("boot/data/cn_nS.dat")
cn <- cutSum(cn_NS) + cutSum(cn_WS)

lf <- ln / cn

ly <- data.frame(Northwestern=rowSums(cay[[1]][-nrow(cay[[1]]),] * lf * lw),
                 Southern=rowSums(cay[[2]][-nrow(cay[[2]]),] * lf * lw),
                 Viking=rowSums(cay[[3]][-nrow(cay[[3]]),] * lf * lw))

e2 <- new.env()
load(file.path("boot/data/FY_SubstockProportions_WGNSSK2024.RData"),e2)

PredFY <- proportions(as.matrix(ly), margin=1)[-(1:12),]
ObsFY <- xtabs(prop ~ year + areaPopulation, e2$substockProps)

png("report/plots/LPfit.png", width=762,height=571)
matplot(yy2,ObsFY,col=c('blue','green','red'),type="p",pch=16,cex=2, xlab="Year", ylab="Proportion")
matplot(yy2,PredFY,type="l",add=TRUE, col=c('blue','green','red'),lty=1,lwd=3)
legend("top",c("Northwest","South","Viking"),col=c('blue',"green","red"),lwd=3,bty="n",ncol=3)
dev.off()

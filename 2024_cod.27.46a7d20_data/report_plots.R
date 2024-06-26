## Prepare plots for report

## Before: 
## After: Report plots

library(icesTAF)
library(ggplot2)
library(reshape2)
library(surveyIndex)
library(FLCore)
library(stockassessment)
library(MASS)
library(boot)
library(grid)

source("utilities.r")

mkdir("report/plots")


## 1 Catch plots 
maxAge <- 7
maxYear <- Inf
minYear <- 1983

# Landings and discards
cn_WS <- read.ices("boot/data/cn_WS.dat")
cn_NS <- read.ices("boot/data/cn_nS.dat")
cn <- cutSum(cn_NS) + cutSum(cn_WS)
cw_WS <- read.ices("boot/data/cw_WS.dat")
cw_NS <- read.ices("boot/data/cw_NS.dat")
cw <- (na2zero(cutMean(cw_NS,cn_NS) * cutSum(cn_NS)) + na2zero(cutMean(cw_WS,cn_WS) * cutSum(cn_WS))) / cn 

ln_WS <- read.ices("boot/data/lf_WS.dat") * cn_WS
ln_NS <- read.ices("boot/data/lf_NS.dat")
ln <- cutSum(ln_NS) + cutSum(ln_WS)
lw_WS <- read.ices("boot/data/lw_WS.dat")
lw_NS <- read.ices("boot/data/lw_NS.dat")
lw <- (na2zero(cutMean(lw_NS,ln_NS) * cutSum(ln_NS)) + na2zero(cutMean(lw_WS,ln_WS) * cutSum(ln_WS))) / ln 

dn_WS <- cn_WS - ln_WS
dn_NS <- cn_NS - ln_NS
dn <- cutSum(dn_NS) + cutSum(dn_WS)
dw_WS <- read.ices("boot/data/dw_WS.dat")
dw_NS <- read.ices("boot/data/dw_NS.dat")
dw <- (na2zero(cutMean(dw_NS,dn_NS) * cutSum(dn_NS)) + na2zero(cutMean(dw_WS,dn_WS) * cutSum(dn_WS))) / dn
dw[dn==0] <- 0

catch <- data.frame(year = as.numeric(row.names(ln)),
                    landings = rowSums(ln * lw),
                    discards = rowSums(dn * dw))
catch <- melt(catch, id="year")
catch$variable <- factor(catch$variable, levels=c("discards","landings"))

ggplot(catch, aes(x=year, y=value)) + geom_area(aes(fill=variable), col="black") +
  scale_fill_manual(values=c("dark grey","light grey")) +
  xlim(1983,2023) +
  theme_bw() +
  theme(legend.position = c(0.8,0.8))
ggsave("report/plots/Fig 4.1 Landings and discards.png")


# Proportion discarded at age
disno <- as.data.frame(dn / cn)[, 1:4]
disno$Year <- rownames(disno)
disno <- melt(disno, id.vars = "Year", value.name = "Proportion")
colnames(disno)[2] <- "Age"
disno$Year <- as.numeric(disno$Year)

ggplot(disno) + geom_line(aes(x = Year, y = Proportion, linetype = Age)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Proportion discarded at age (by number)")
ggsave("report/plots/Fig 4.1 Numbers-at-age discarded.png")

# Proportion of numbers discarded
disnotot <- as.data.frame(rowSums(dn) / rowSums(cn))
disnotot$Year <- rownames(disnotot)
colnames(disnotot)[1] <- "Proportion"
disnotot$Year <- as.numeric(disnotot$Year)

ggplot(disnotot) + geom_line(aes(x = Year, y = Proportion)) +
  theme_bw() + 
  ggtitle("Proportion of all cod discarded (by number)")
ggsave("report/plots/Fig 4.1 Numbers discarded.png")

# Proportion discarded by weight
diswt <- as.data.frame(rowSums(dn * dw) / rowSums(cn * cw))
diswt$Year <- rownames(diswt)
colnames(diswt)[1] <- "Proportion"
diswt$Year <- as.numeric(diswt$Year)

ggplot(diswt) + geom_line(aes(x = Year, y = Proportion)) + 
  theme_bw() +
  ggtitle("Proportion of all cod discarded (by weight)")
ggsave("report/plots/Fig 4.1 Weight discarded.png")

# Catch weights
lwt <- as.data.frame(lw)
lwt$Year <- rownames(lwt)
lwt$Type <- "Landings"
dwt <- as.data.frame(dw)
dwt$Year <- rownames(dwt)
dwt$Type <- "Discards"
cwt <- as.data.frame(cw)
cwt$Year <- rownames(cwt)
cwt$Type <- "Catch"
mnwt <- rbind(lwt,dwt,cwt)
mnwt <- melt(mnwt, id.vars = c("Year","Type"), value.name = "Weight")
colnames(mnwt)[3] <- "Age"
mnwt$Year <- as.numeric(mnwt$Year)
mnwt$Type <- factor(mnwt$Type, levels = c("Landings","Discards","Catch"))

ggplot(mnwt) + geom_line(aes(x = Year, y = Weight, linetype = Age)) + #, color = Age)) +
  facet_wrap(Type~.) +
  theme_bw()
ggsave("report/plots/Fig 4.2 Catch weights.png")


## 2 Biological data

# Stock weights
swNW <- read.ices("output/sw_Northwest.dat")
swS <- read.ices("output/sw_South.dat")
swV <- read.ices("output/sw_Viking.dat")

png("report/plots/Fig 4.2 SW inputs.png", width=853, height=500)
par(mfrow=c(1,3))
matplot(1983:(1983+nrow(swV)-1), swV, xlab="Year", ylab="Viking stock weights")
matplot(1983:(1983+nrow(swV)-1), swV, type="l", add=TRUE, lty="solid", lwd=2)
matplot(1983:(1983+nrow(swNW)-1), swNW, xlab="Year", ylab="Northwest stock weights")
matplot(1983:(1983+nrow(swNW)-1), swNW, type="l", add=TRUE, lty="solid", lwd=2)
matplot(1983:(1983+nrow(swS)-1), swS, xlab="Year", ylab="South stock weights")
matplot(1983:(1983+nrow(swS)-1), swS, type="l", add=TRUE, lty="solid", lwd=2)
dev.off()

# Maturity
moNW <- read.ices("output/mo_northwest.dat")
moS <- read.ices("output/mo_south.dat")
moV <- read.ices("output/mo_viking.dat")

png("report/plots/Fig 4.2 MO inputs.png", width=853, height=500)
par(mfrow=c(1,3))
matplot(1983:(1983+nrow(moNW)-1), moNW, xlab="Year", ylab="Northwest maturity")
matplot(1983:(1983+nrow(moNW)-1), moNW, type="l", add=TRUE, lty="solid", lwd=2)
matplot(1983:(1983+nrow(moS)-1), moS, xlab="Year", ylab="South stock maturity")
matplot(1983:(1983+nrow(moS)-1), moS, type="l", add=TRUE, lty="solid", lwd=2)
matplot(1983:(1983+nrow(moV)-1), moV, xlab="Year", ylab="Viking stock maturity")
matplot(1983:(1983+nrow(moV)-1), moV, type="l", add=TRUE, lty="solid", lwd=2)
dev.off()

# Natural mortality
png("report/plots/Fig 4.2 NM input.png", width=600, height=500)
nm <- read.ices("output/nm_2023.dat")
matplot(1983:(1983+nrow(nm)-1), nm[, 1:4], xlab="Year", ylab="Natural mortality")
matplot(1983:(1983+nrow(nm)-1), nm[, 1:4], type="l", add=TRUE, lty="solid", lwd=2)
dev.off()


## 3 Survey plots
load("boot/data/commonGrid.RData")
Q1grid <- Q34grid <- gridd

# Q1 plots
load("model/Q1models.Rdata")
load("model/Q1retro.Rdata")
# Maps
png("report/plots/Fig 4.3 Map Q1.png", width=750, height=538)
surveyIdxPlots(models[[1]],ddQ1,predD=Q1grid,select="absolutemap",year=as.character(2020:2024),colors=rev(heat.colors(8)),
               par=list(mfrow=c(7,5),mar=c(0,0,2,0),oma=c(0,0,2,0)),map.cex=1.3,legend=FALSE,cols=1:7, plotByAge = FALSE)
dev.off()
# Retro
pdf("report/plots/Fig 4.3 retroQ1.pdf",width=12,height=8)
plot(SIretro,base=models[[1]])
SImohn = mohn.surveyIdx(SIretro,base=models[[1]])
title(paste0("Retrospective analysis Q1, avg mohn:",round(mean(SImohn),2)),outer=TRUE)
dev.off()
# Sub-stock indices
png("report/plots/Fig 4.3 Viking index Q1.png", width=952, height=797)
surveyIdxPlots(models[[2]],ddQ1,cols=1:7,par=list(mfrow=c(4,2)),legend=FALSE,select="index",plotByAge=FALSE)
dev.off()
png("report/plots/Fig 4.3 Northwest index Q1.png", width=952, height=797)
surveyIdxPlots(models[[3]],ddQ1,cols=1:7,par=list(mfrow=c(4,2)),legend=FALSE,select="index",plotByAge=FALSE)
dev.off()
png("report/plots/Fig 4.3 South index Q1.png", width=952, height=797)
surveyIdxPlots(models[[4]],ddQ1,cols=1:7,par=list(mfrow=c(4,2)),legend=FALSE,select="index",plotByAge=FALSE) 
dev.off()

indices <- list()
for (i in 1:3){
  indices[[i]] <- FLIndex(index=FLQuant(t(models[[i+1]]$idx), dimnames=list(age=1:7, year=1983:2024)),
                          effort=FLQuant(rep(1,length(1983:2024)), dimnames=list(age='all', year=1983:2024)),
                          name=paste0(names(models)[i+1], " cod indices Q1"))
  indices[[i]]@catch.n <- indices[[i]]@index
  indices[[i]]@range[6] <- 0.25/2
  indices[[i]]@range[7] <- 0.25/2
}
rm(models,ddQ1, Q1grid, SIretro)

load("model/Q34models.Rdata")
load("model/Q34retro.Rdata")
# Maps
png("report/plots/Fig 4.3 Map Q34.png", width=750, height=538)
surveyIdxPlots(models[[1]],ddQ34,predD=Q34grid,select="absolutemap",year=as.character(2020:2023),colors=rev(heat.colors(8)),
               par=list(mfrow=c(7,4),mar=c(0,0,2,0),oma=c(0,0,2,0)),map.cex=1.3,legend=FALSE,cols=1:7, plotByAge = FALSE)
dev.off()
# Retro
pdf("report/plots/Fig 4.3 retroQ34.pdf",width=12,height=8)
plot(SIretro,base=models[[1]])
SImohn = mohn.surveyIdx(SIretro,base=models[[1]])
title(paste0("Retrospective analysis Q3, avg mohn:",round(mean(SImohn),2)),outer=TRUE)
dev.off()
# Indices
png("report/plots/Fig 4.3 Index Q34.png", width=952, height=797)
surveyIdxPlots(models[[1]],ddQ34,cols=2:8,par=list(mfrow=c(4,2)),legend=FALSE,select="index",plotByAge=FALSE)
dev.off()
png("report/plots/Fig 4.3 Viking Rindex.png", width=952, height=450)
surveyIdxPlots(models[[2]],ddQ34,cols=1,par=list(mfrow=c(1,1)),legend=FALSE,select="index",plotByAge=FALSE, main="Viking")
dev.off()
png("report/plots/Fig 4.3 Northwest Rindex.png", width=952, height=450)
surveyIdxPlots(models[[3]],ddQ34,cols=1,par=list(mfrow=c(1,1)),legend=FALSE,select="index",plotByAge=FALSE, main="Northwest")
dev.off()
png("report/plots/Fig 4.3 South Rindex.png", width=952, height=450)
surveyIdxPlots(models[[4]],ddQ34,cols=1,par=list(mfrow=c(1,1)),legend=FALSE,select="index",plotByAge=FALSE, main="South")
dev.off()

indices[[4]] <- FLIndex(index=FLQuant(t(models[[1]]$idx[,-1]), dimnames=list(age=1:7, year=1992:2023)),
                        effort=FLQuant(rep(1,length(1992:2023)), dimnames=list(age='all', year=1992:2023)),
                        name="NS cod indices Q3 & Q4")
indices[[4]]@catch.n <- indices[[4]]@index
indices[[4]]@range[6] <- 1.5/2
indices[[4]]@range[7] <- 1.5/2

for (i in 5:7) indices[[i]] <- FLIndex(index=FLQuant(t(models[[i-3]]$idx[,1]), dimnames=list(age=1, year=1993:2024)),
                                       name=paste0(names(models)[i-3], " recruitment index"))
rm(models,ddQ34, Q34grid, SIretro)


## 4 Data analyses

# Survey-based
cod.tun <- FLIndices(indices)

# Catch curve analyses
png("report/plots/Fig 4.4 CC Viking Q1.png", width=952, height=554)
par(mfrow=c(2,2))
single_cpue_plot_year(FLIndices(cod.tun[[1]]))
single_cpue_plot_cohort(FLIndices(cod.tun[[1]]))
single_catchcurve_index(FLIndices(cod.tun[[1]]))
single_catchcurvegrad_index(FLIndices(cod.tun[[1]]),list(c(2,4)))
dev.off()

png("report/plots/Fig 4.4 CC Northwest Q1.png", width=952, height=554)
par(mfrow=c(2,2))
single_cpue_plot_year(FLIndices(cod.tun[[2]]))
single_cpue_plot_cohort(FLIndices(cod.tun[[2]]))
single_catchcurve_index(FLIndices(cod.tun[[2]]))
single_catchcurvegrad_index(FLIndices(cod.tun[[2]]),list(c(2,4)))
dev.off()

png("report/plots/Fig 4.4 CC South Q1.png", width=952, height=554)
par(mfrow=c(2,2))
single_cpue_plot_year(FLIndices(cod.tun[[3]]))
single_cpue_plot_cohort(FLIndices(cod.tun[[3]]))
single_catchcurve_index(FLIndices(cod.tun[[3]]))
single_catchcurvegrad_index(FLIndices(cod.tun[[3]]),list(c(2,4)))
dev.off()

png("report/plots/Fig 4.4 CC Q34.png", width=952, height=554)
par(mfrow=c(2,2))
single_cpue_plot_year(FLIndices(cod.tun[[4]]))
single_cpue_plot_cohort(FLIndices(cod.tun[[4]]))
single_catchcurve_index(FLIndices(cod.tun[[4]]))
single_catchcurvegrad_index(FLIndices(cod.tun[[4]]),list(c(2,4)))
dev.off()

# Internal consistency
png("report/plots/Fig 4.5 IC Viking.png", width=952, height=600)
corplotswithin_index(cod.tun,1)
dev.off()
png("report/plots/Fig 4.5 IC Northwest.png", width=952, height=600)
corplotswithin_index(cod.tun,2)
dev.off()
png("report/plots/Fig 4.5 IC South.png", width=952, height=600)
corplotswithin_index(cod.tun,3)
dev.off()
png("report/plots/Fig 4.5 IC Q34.png", width=952, height=600)
corplotswithin_index(cod.tun,4)
dev.off()

# External consistency - recruitment
png("report/plots/Fig 4.5 EC Viking.png", width=554, height=554)
corplotsbetween_index1(cod.tun,5,1)
dev.off()
png("report/plots/Fig 4.5 EC Northwest.png", width=554, height=554)
corplotsbetween_index1(cod.tun,6,2)
dev.off()
png("report/plots/Fig 4.5 EC South.png", width=554, height=554)
corplotsbetween_index1(cod.tun,7,3)
dev.off()


# Catch-based
# Catch numbers
catch <- as.data.frame(cn)
catch$Year <- rownames(catch)
catch <- melt(catch, id.vars = "Year", value.name = "Catch")
colnames(catch)[2] <- "Age"
catch$Year <- as.numeric(catch$Year)

ggplot(catch) + geom_point(aes(x = Year, y = Age, size = Catch), fill = 'gray', pch = 21) +
  scale_size(range = c(0.1, 10)) +
  theme_bw() + 
  ggtitle("Catch numbers-at-age") +
  guides(size="none")
ggsave("report/plots/Fig 4.6a_cn.png")

# Catch numbers (mean standardised)
catch <- as.data.frame(cn)
catch <- sweep(catch, 1, rowMeans(catch), `/`)
catch$Year <- rownames(catch)
catch <- melt(catch, id.vars = "Year", value.name = "Catch")
colnames(catch)[2] <- "Age"
catch$Year <- as.numeric(catch$Year)

ggplot(catch) + geom_point(aes(x = Year, y = Age, size = Catch), fill = 'gray', pch = 21) +
  scale_size(range = c(0.1, 10)) +
  theme_bw() +
  ggtitle("Catch numbers-at-age (mean standardised)")

# Standardised proportions-at-age
catch <- as.data.frame(cn)
catch <- sweep(catch, 1, rowMeans(catch), `/`)
mns <- colMeans(catch)
sds <- apply(catch, 2, sd)
catch <- sweep(catch, 2, mns, `-`)
catch <- sweep(catch, 2, sds, `/`)
catch$Year <- rownames(catch)
catch <- melt(catch, id.vars = "Year", value.name = "Catch")
colnames(catch)[2] <- "Age"
catch$Year <- as.numeric(catch$Year)
catch$Sign <- ifelse(catch$Catch < 0, "-", "+")
catch$Catch <- abs(catch$Catch) # Otherwise all negative residuals end up as small bubbles

ggplot(catch) + geom_point(aes(x = Year, y = Age, size = Catch, fill = Sign), pch = 21) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_size(range = c(0.1, 15)) +
  theme_bw() +
  ggtitle("Standardised proportions-at-age") +
  guides(size="none", fill="none")
ggsave("report/plots/Fig 4.6b_stcn.png")

# Catch curves, gradients & correlations
cn.flq <- FLQuant(t(cn), dimnames=list(age=1:7, year = 1983:2023))
cod <- FLStock(catch.n=cn.flq)
png("report/plots/Fig 4.7 CC Catch.png", width=952, height=554)
catchcurve(cod)
dev.off()
png("report/plots/Fig 4.7 CC grad Catch.png", width=952, height=554)
catchcurvegrad(cod, c(2,4))
dev.off()
png("report/plots/Fig 4.8 Catch corr.png", width=952, height=554)
plot.index.corr(list(FLIndex(catch.n=cod@catch.n, name="commercial catch")), wndows = FALSE)
dev.off()

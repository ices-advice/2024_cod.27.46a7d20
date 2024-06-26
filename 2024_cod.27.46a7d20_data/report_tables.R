## Prepare tables for report

## Before: 
## After: Report tables

library(icesTAF)
library(stockassessment)
library(icesAdvice)
library(xlsx)

source("utilities.R")

mkdir("report/tables")


## 1 Catch tables
maxAge <- 7
maxYear <- Inf
minYear <- 1983

# Landings and discards
cn_WS <- read.ices("boot/data/cn_WS.dat")
cn_NS <- read.ices("boot/data/cn_NS.dat")
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

write.xlsx(ln, file = "report/tables/Tab 4.2-4.xlsx", sheetName = "ln")
write.xlsx(lw, file = "report/tables/Tab 4.2-4.xlsx", sheetName = "lw", append = TRUE)
write.xlsx(dn, file = "report/tables/Tab 4.2-4.xlsx", sheetName = "dn", append = TRUE)
write.xlsx(dw, file = "report/tables/Tab 4.2-4.xlsx", sheetName = "dw", append = TRUE)
write.xlsx(cn, file = "report/tables/Tab 4.2-4.xlsx", sheetName = "cn", append = TRUE)
write.xlsx(cw, file = "report/tables/Tab 4.2-4.xlsx", sheetName = "cw", append = TRUE)


## 2 Biological data

# Stock weights
swNW <- read.ices("output/sw_Northwest.dat")
swS <- read.ices("output/sw_South.dat")
swV <- read.ices("output/sw_Viking.dat")

write.xlsx(swNW, file = "report/tables/Tab_sw_input.xlsx", sheetName = "Northwest")
write.xlsx(swS, file = "report/tables/Tab_sw_input.xlsx", sheetName = "South", append=TRUE)
write.xlsx(swV, file = "report/tables/Tab_sw_input.xlsx", sheetName = "Viking", append=TRUE)

# Maturity
moNW <- read.ices("output/mo_northwest.dat")
moS <- read.ices("output/mo_south.dat")
moV <- read.ices("output/mo_viking.dat")

write.xlsx(moNW, file = "report/tables/Tab_mo_input.xlsx", sheetName = "Northwest")
write.xlsx(moS, file = "report/tables/Tab_mo_input.xlsx", sheetName = "South", append=TRUE)
write.xlsx(moV, file = "report/tables/Tab_mo_input.xlsx", sheetName = "Viking", append=TRUE)

# Natural mortality
nm <- read.ices("output/nm_2023.dat")
write.xlsx(nm, file = "report/tables/Tab_nm_input.xlsx", sheetName = "Natural mortality")


## 3 Indices

# Q1
load("model/Q1models.Rdata")
indices = lapply(models,function(x) x$idx )
los = lapply(models,function(x) x$lo )
ups = lapply(models,function(x) x$up )

write.xlsx(cbind(indices[[3]],(log(ups[[3]]) - log(los[[3]]))/4), file = "report/tables/Tab_indices.xlsx", sheetName = "Northwest")
write.xlsx(cbind(indices[[4]],(log(ups[[4]]) - log(los[[4]]))/4), file = "report/tables/Tab_indices.xlsx", sheetName = "South", append = TRUE)
write.xlsx(cbind(indices[[2]],(log(ups[[2]]) - log(los[[2]]))/4), file = "report/tables/Tab_indices.xlsx", sheetName = "Viking", append = TRUE)
rm(ddQ1,models,indices,los,ups)

# Q34
load("model/Q34models.Rdata")
indices = lapply(models,function(x) x$idx )
los = lapply(models,function(x) x$lo )
ups = lapply(models,function(x) x$up )
write.xlsx(cbind(indices[[1]][,-1],(log(ups[[1]][,-1]) - log(los[[1]][,-1]))/4), file = "report/tables/Tab_indices.xlsx", sheetName = "Q34", append = TRUE)

Rindices <- data.frame("Northwest"=indices[[3]][,1],
                       "South"=indices[[4]][,1],
                       "Viking"=indices[[2]][,1],
                       "Northwest"=(log(ups[[3]][,1]) - log(los[[3]][,1]))/4,
                       "South"=(log(ups[[4]][,1]) - log(los[[4]][,1]))/4,
                       "Viking"=(log(ups[[2]][,1]) - log(los[[2]][,1]))/4)
row.names(Rindices) <- as.numeric(row.names(Rindices))+1
write.xlsx(Rindices, file = "report/tables/Tab_indices.xlsx", sheetName = "Recruit", append = TRUE)

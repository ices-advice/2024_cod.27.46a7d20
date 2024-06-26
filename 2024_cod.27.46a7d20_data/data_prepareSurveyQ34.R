library(surveyIndex)
library(mapplots)
library(XLConnect)

load("boot/data/codData100424.Rdata")
ICESarea <- read.csv("boot/data/ICESrect.csv")
wb=XLConnect::loadWorkbook("boot/data/Area_NSCOD3_table_6a_inshore_offshore.xlsx")

dQ34 <- c(dQ3,dQ4)

dQ1 = NULL
dQ3 = NULL
dQ4 = NULL

dQ34 = subset(dQ34, !is.na(lon))

years=1990:2023
dQ34 = subset(dQ34, Year %in% years)

mkdir("data/plots/Q34")
for (yr in years){
  png(paste0("data/plots/Q34/",yr,".png"), width = 15, height = 15, units = "cm", res=600)
  plot(subset(dQ34, Year==yr), col=subset(dQ34, Year==yr)$Survey)
  dev.off()
}

# Check surveys and years
xtabs(~Year+Survey,data=dQ34[[2]])

# Fill in missing StatRecs
dQ34$StatRec <- as.character(dQ34$StatRec)
sel <- which(is.na(dQ34$StatRec)) 
dQ34$StatRec[sel] <- ices.rect2(dQ34$lon[sel], dQ34$lat[sel])
dQ34$StatRec <- factor(dQ34$StatRec)

# Add ICES areas
ICESarea$ICESNAME <- ices.rect2(ices.rect(ICESarea$ICESNAME)$lon, ices.rect(ICESarea$ICESNAME)$lat)
dQ34$StatRec <- ices.rect2(ices.rect(dQ34$StatRec)$lon, ices.rect(dQ34$StatRec)$lat)
dQ34$ICESarea <- ICESarea$Area_27[match(dQ34$StatRec,ICESarea$ICESNAME)]

# Subset to stock definition
dQ34 <- subset(dQ34, ICESarea %in% c("4.a","4.b","4.c","6.a","3.a.20"))
plot(dQ34, col=dQ34$Survey)

# Add cod subareas
areadef=XLConnect::readWorksheet(wb,sheet=1)
areadef$SubArea <- NA
areadef$SubArea[areadef$areaPopulation=="Viking"] <- 1
areadef$SubArea[is.element(areadef$areaPopulation, c("Northwestern offshore"))] <- 2
areadef$SubArea[is.element(areadef$areaPopulation, c("Northwestern inshore"))] <- 3
areadef$SubArea[areadef$areaPopulation=="Southern"] <- 4
#areadef$SubArea[is.element(areadef$Statistical.Rectangle, c("40E4","40E5","39E4","39E5"))] <- "C"

dQ34$subArea <- areadef$SubArea[match(dQ34$StatRec, areadef$Statistical.Rectangle)]
plot(dQ34, col=dQ34$subArea)
points(subset(dQ34, is.na(subArea)), col="red")


# Take care of Scottish vessels - make it Ship+Survey
unique(subset(dQ34, Country=="GB-SCT")[[2]]$Ship)
xtabs(~Year+Ship, data=subset(dQ34, Ship %in% c("749S","748S"))[[2]])

dQ34$Ship = as.character(dQ34$Ship)
sel = c(grep("749S",dQ34$Ship), grep("748S",dQ34$Ship))
dQ34$Ship[sel] = paste(dQ34$Ship[sel],as.numeric(dQ34$Survey[sel]),sep=":")
dQ34$Ship = factor(dQ34$Ship)
xtabs(~Ship+Survey,dQ34[[2]])


# Check gears
xtabs(~Year+Gear, dQ34[[2]]) # Mix of gears in 1991
dQ34 = subset(dQ34, !Year %in% 1990:1991) 


# Check for enough age data
xtabs(NoAtALK~Year+Age,data=dQ34[[1]])


## impute missing depths
dmodelQ34=gam(log(Depth) ~ s(lon,lat,k=200),data=dQ34[[2]])
selQ34=subset(dQ34,is.na(Depth))
selQ34$Depth=0; ## Guard against NA-error
dQ34$Depth[is.na(dQ34$Depth)]=exp(predict(dmodelQ34,newdata=selQ34[[2]]))

save(dQ34,file="data/Q34preprocessed.RData")

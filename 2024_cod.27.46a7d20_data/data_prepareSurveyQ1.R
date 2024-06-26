library(surveyIndex)
library(mapplots)
library(XLConnect)

load("boot/data/codData100424.Rdata")
ICESarea <- read.csv("boot/data/ICESrect.csv")
wb=XLConnect::loadWorkbook("boot/data/Area_NSCOD3_table_6a_inshore_offshore.xlsx")

dQ3 = NULL
dQ4 = NULL

dQ1 = subset(dQ1, !is.na(lon)) # Only one haul in 2024 and no cod caught

years=1983:2024
dQ1 = subset(dQ1, Year %in% years)

mkdir("data/plots/Q1")
for (yr in years){
  png(paste0("data/plots/Q1/",yr,".png"), width = 15, height = 15, units = "cm", res=600)
  plot(subset(dQ1, Year==yr), col=subset(dQ1, Year==yr)$Survey)
  dev.off()
}

# Check surveys and years
xtabs(~Year+Survey,data=dQ1[[2]])

# Fill in missing StatRecs
sel <- which(is.na(dQ1$StatRec))
dQ1$StatRec[sel] <- ices.rect2(dQ1$lon[sel], dQ1$lat[sel])

# Add ICES areas
ICESarea$ICESNAME <- ices.rect2(ices.rect(ICESarea$ICESNAME)$lon, ices.rect(ICESarea$ICESNAME)$lat)
dQ1$StatRec <- ices.rect2(ices.rect(dQ1$StatRec)$lon, ices.rect(dQ1$StatRec)$lat) # Guard against exponents
dQ1$ICESarea <- ICESarea$Area_27[match(dQ1$StatRec,ICESarea$ICESNAME)]

# Subset to stock definition
dQ1 <- subset(dQ1, ICESarea %in% c("4.a","4.b","4.c","6.a","3.a.20"))
plot(dQ1, col=dQ1$Survey)

# Add cod subareas
areadef=XLConnect::readWorksheet(wb,sheet=1)
areadef$SubArea <- NA
areadef$SubArea[areadef$areaPopulation=="Viking"] <- 1
areadef$SubArea[is.element(areadef$areaPopulation, c("Northwestern offshore"))] <- 2
areadef$SubArea[is.element(areadef$areaPopulation, c("Northwestern inshore"))] <- 3
areadef$SubArea[areadef$areaPopulation=="Southern"] <- 4
#areadef$SubArea[is.element(areadef$Statistical.Rectangle, c("40E4","40E5","39E4","39E5"))] <- "C"

dQ1$subArea <- areadef$SubArea[match(dQ1$StatRec, areadef$Statistical.Rectangle)]
plot(dQ1, col=dQ1$subArea)
points(subset(dQ1, is.na(subArea)), col="red")


# Take care of Scottish vessels - make it Ship+Survey
unique(subset(dQ1, Country=="GB-SCT")[[2]]$Ship)
xtabs(~Year+Ship, data=subset(dQ1, Ship %in% c("74EX","749S","748S"))[[2]])
unique(subset(dQ1, Ship=="74EX")$Survey) # NS-IBTS only

dQ1$Ship = as.character(dQ1$Ship)
sel = c(grep("749S",dQ1$Ship), grep("748S",dQ1$Ship))
dQ1$Ship[sel] = paste(dQ1$Ship[sel],as.numeric(dQ1$Survey[sel]),sep=":")
dQ1$Ship = factor(dQ1$Ship)
xtabs(~Ship+Survey,dQ1[[2]])


# Check gears
xtabs(~Year+Gear, dQ1[[2]]) # 63 H18 hauls excluded
dQ1 = subset(dQ1, Gear=="GOV")


# Check for enough age data
xtabs(NoAtALK~Year+Age,data=dQ1[[1]])
dQ1[[1]]=subset(dQ1[[1]],Age>0)


## impute missing depths
dmodelQ1=gam(log(Depth) ~ s(lon,lat,k=200),data=dQ1[[2]])
selQ1=subset(dQ1,is.na(Depth))
selQ1$Depth=0; ## Guard against NA-error
dQ1$Depth[is.na(dQ1$Depth)]=exp(predict(dmodelQ1,newdata=selQ1[[2]]))

save(dQ1,file="data/Q1preprocessed.RData")

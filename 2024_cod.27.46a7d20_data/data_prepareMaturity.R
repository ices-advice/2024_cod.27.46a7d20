set.seed(12345)

library(DATRAS)
library(plyr)
library(dplyr)
library(mapplots)
library(ggplot2)
library(xlsx)

load("boot/data/codData100424.Rdata")
fil <- paste0("boot/data/Area_NSCOD3_table_6a_inshore_offshore.xlsx")

working_year <- 2024
years <- 1983:working_year
ages <- c(1:7)

# Species specific parameters
genus = "Gadus"
bfamily = "morhua"

# Subareas
wb=XLConnect::loadWorkbook(fil)
areadef=XLConnect::readWorksheet(wb,sheet=1)
areadef$SubArea <- NA
areadef$SubArea[areadef$areaPopulation=="Viking"] <- "V"
areadef$SubArea[is.element(areadef$areaPopulation, c("Northwestern offshore","Northwestern inshore"))] <- "NW"
areadef$SubArea[areadef$areaPopulation=="Southern"] <- "S"
#areadef$SubArea[is.element(areadef$Statistical.Rectangle, c("40E4","40E5","39E4","39E5"))] <- "C" # Inshore
any(is.na(areadef$SubArea))

# Load DATRAS data
dQ3=NULL; dQ4=NULL
cod <- subset(dQ1, !is.na(Maturity))

# Recode Maturity 0 == Immature, 1 == Mature
cod[["CA"]]$Maturity = as.character(cod[["CA"]]$Maturity)
cod[["CA"]]$Maturity = revalue(cod[["CA"]]$Maturity, c("61"="0", "62"="1", "63"="1", "64"="1", "65"="1", "1"="0", "2"="1", "3"="1", "4"="1", "5"="1", "I"="0", "M"="1", "A"="0", "B"="1", "Ba"="0", "Bb"="1", "C"="1", "Ca"="1", "Cb"="1", "D"="1", "Da"="1", "Db"="1", "E"="1"))
cod[["CA"]] = filter(cod[["CA"]], Maturity!=66 & Maturity!=6 & Maturity!="F") # Remove abnormal
table(cod[["CA"]]$Maturity) # Check for new scales!
cod[["CA"]]$Maturity = as.numeric(cod[["CA"]]$Maturity)

# Create +group
cod[["CA"]]$Age[cod[["CA"]]$Age >= max(ages)] <- max(ages)

# Add raw number of observed individuals per length group to HH data
cod = addSpectrum(cod,by=1)

# Standardise for an effort of 60 mins
cod[["HH"]]$N <- cod[["HH"]]$N / cod[["HH"]]$HaulDur * 60

# Subset
cod <- subset(cod, Year %in% years) 
cod <- subset(cod, Age %in% ages) 

# Total numbers per haul
cod[["HH"]]$nbL <- rowSums(cod[["HH"]]$N)

# Fill in missing StatRecs
cod$StatRec <- as.character(cod$StatRec)
sel <- which(is.na(cod$StatRec))
cod$StatRec[sel] <- ices.rect2(cod$lon[sel], cod$lat[sel])
any(is.na(cod$StatRec))
cod$StatRec <- ices.rect2(ices.rect(cod$StatRec)$lon, ices.rect(cod$StatRec)$lat) # Some rects not matching due to exponential notation
cod$StatRec <- factor(cod$StatRec)

# Assign cod hauls to population subareas
cod$SubArea <- areadef$SubArea[match(cod$StatRec,areadef$Statistical.Rectangle)]

# Add SubArea column to CA data and set to character
cod[["HH"]] = transform(cod[["HH"]], SubArea=as.character(SubArea), haul.id=as.character(haul.id))
SUBA = as.data.frame(cbind(haul.id = cod[["HH"]]$haul.id, SubArea = cod[["HH"]]$SubArea))
cod[["CA"]] = merge(cod[["CA"]], SUBA, by=c("haul.id"))

plot(subset(cod, !is.na(lon), is.na(SubArea))) # All outside stock definition
cod[["HH"]]$haul.id <- factor(cod[["HH"]]$haul.id) # Otherwise haul.id disappears from CA
cod <- subset(cod, !is.na(SubArea))

mkdir("data/tables")
write.xlsx(ddply(cod[[1]], .(Year, Age, SubArea), summarise, N.CA=sum(NoAtALK), Mature=sum(Maturity)), file="data/tables/sampling.xlsx", sheetName="Maturity", row.names=FALSE)


##############################
## Region specific catch rates
##############################

subHH <- droplevels(cod[["HH"]])

subHH$SubArea <- as.factor(subHH$SubArea)

# Mean catch within StatRecs and then across StatRecs within SubAreas
rsw <- ddply(subHH,. (Year,StatRec), summarise, nbL=mean(nbL), SubArea=unique(SubArea))
rsw <- ddply(rsw,. (Year,SubArea), summarise, nbL=mean(nbL))

# Number of rectangles in SubArea / Total number of rectangles
props <- count(areadef,SubArea)
props$q <- props$n / sum(props$n)

# Area-weighted mean catches / sum of area-weighted mean catches
rsw <- merge(rsw, props[,c(1,3)], by="SubArea", all.x=TRUE)
rsw$wnbL <- rsw$nbL * rsw$q
rsw <- merge(rsw, ddply(rsw,.(Year), summarise, wtSum=sum(wnbL)), by="Year", all.x = TRUE)
rsw$wtC <- rsw$wnbL / rsw$wtSum

# test <- ddply(rsw,.(Year), summarise, swtC=sum(wtC)) # Should sum to 1


######################
# Statistical weights 
######################

# Numbers at age- and length- from the CA data (by haul)
smalk <- cod[["CA"]]
smalk$SubArea <- as.factor(smalk$SubArea)
smalk$length <- floor(smalk$LngtCm)
smalk <- ddply(smalk,.(haul.id,Age,length), summarise, NoAtAlk=sum(NoAtALK), NoMature=sum(Maturity), Year=unique(Year), SubArea=unique(SubArea))

# Numbers at length from the HH data (by haul)
wL <- as.data.frame(cod[["HH"]]$N)
wL$length <- as.numeric(sub('.([^,]+),.*','\\1', wL$sizeGroup))

# Calculate raising factors as no. fish sampled / no. fish measured
smalk <- merge(smalk, wL[,c(1,3:4)], by=c("haul.id", "length"), all.x=TRUE)
table(smalk$NoAtAlk<=smalk$Freq) # Hmmm needs some explanation!
smalk$nbL <- pmax(smalk$NoAtAlk,smalk$Freq)
smalk$rg <- smalk$nbL / smalk$NoAtAlk

# Sum raising factors (Ra=sum(rg)) and the number of sampled fish (ma) across ages...
# And calculate statistical weights as ma x rg / Ra
smalk <- merge(smalk, ddply(smalk,.(haul.id,Age), summarise, Ra=sum(rg), ma=sum(NoAtAlk)), by=c("haul.id","Age"), all.x=TRUE)
smalk$wi <- smalk$ma * smalk$rg / smalk$Ra

# # From WKMOG: "For each age group the sum of statistical weights will equal the number of fish with biological data"
# test <- ddply(smalk,.(haul.id,Age), summarise, NoAtAlk=sum(NoAtAlk), wi=sum(wi))
# all(test$NoAtAlk==test$wi)
# test[test$NoAtAlk!=test$wi,] # They look the same to me!

mat_data <- merge(smalk, rsw[,-c(3:6)], by=c("Year","SubArea"), all.x=TRUE)
save(mat_data, file="data/mat_data.Rdata")

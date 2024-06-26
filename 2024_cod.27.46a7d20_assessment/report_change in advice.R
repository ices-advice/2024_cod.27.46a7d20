# Based on script developed at WGNSSK 2023 (Harriet Cole, Iago Mosqueria, Yves Reecht)

# Inputs:
# Compare_forecast_assumptions.csv - table of recruitment, SSB, Fbar and total catch from this year's WG and last year's WG for the data year, intermediate year and advice year.
# Compare_forecast_stockwts.csv - table of stock weights at age used in the forecast from this year's WG and last year's WG for the data year, intermediate year and advice year.
# Compare_forecast_selectivity.csv - table of selectivity at age used in the forecast from this year's WG and last year's WG for the data year, intermediate year and advice year.
# Compare_forecast_N_at_age.csv - table of stock numbers at age used in the forecast from this year's WG and last year's WG for the data year, intermediate year and advice year.
# Compare_forecast_B_at_age.csv - table of stock biomass at age used in the forecast from this year's WG and last year's WG for the data year, intermediate year and advice year.
# Now_assessment_N_at_age.csv - table of stock numbers at age from THIS year's assessment. Usually a standard output for the WG report.
# Now_assessment_B_at_age.csv - table of stock biomass at age from THIS year's assessment. Multiply "now_assessment_N_at_age.csv" by stock weights at age
# Prev_assessment_N_at_age.csv - table of stock numbers at age from LAST year's assessment. Usually a standard output for the WG report.
# Prev_assessment_B_at_age.csv - table of stock biomass at age from LAST year's assessment. Multiply "now_assessment_N_at_age.csv" by stock weights at age
  
# setup ------------------------------------####
rm(list=ls())
graphics.off()

# libraries
library(icesTAF)
library(reshape2)
library(ggplot2)
library(xlsx)

# settings
stock.name <- paste("cod.27.46a7d20",c("NW","SO","VI"),sep=".")
ages <- 1:7 # ages for stock
ay <- 2024 # assessment year
substocks <- c("Northwest","South","Viking")

# Source as described in input files for this script. Must be an exact match as these are used for filtering later.
source.now <-"WGNSSK 2024"
source.prev <-"WGNSSK 2023"

# Compare forecast assumptions ------------------------------------------------------------

for (s in 1:3){
  # read in and format forecast assumptions
  dat <- read.xlsx("output/Forecast_assumptions.xlsx", sheetName = substocks[s])
  dat <- dat[dat$Year >(ay-2),]
  dat$Type <- factor(dat$Type,levels=c("Data year","Intermediate year","Advice year"))
  dat$Variable <- factor(dat$variable,levels=c("SSB","Fbar","Total.catch","Recruitment"))
  
  # plot
  png(paste0("report/plots/Fig 4.2", s-1,"e Compare_forecast_assumptions_",substocks[s],".png"),width = 11, height = 7, units = "in", res = 600)
  p1 <- ggplot(dat,aes(x=Year,y=value,colour=Source,shape=Type))+geom_point(size=3)+facet_wrap(~Variable,scales="free_y")+
    theme_bw()+labs(x="",y="",colour="",shape="")+ scale_shape_manual(values=c(16, 2, 0))+ggtitle(stock.name[s])
  print(p1)
  dev.off()
  
  tab <- dat[dat$Year %in% 2023:2024, -c(3,4)]
  tab <- dcast(tab, Variable+Year~Source,value.var="value")
  tab <- tab[c(7,8,5,1:4),]
  tab <- tab[,c(1,2,4,3)]
  
  if (s==1){
    write.xlsx(tab,file = "report/tables/Tab 4.17 Compare_forecast_assumptions.xlsx",row.names=F, sheetName = substocks[s])
  }else{
    write.xlsx(tab,file = "report/tables/Tab 4.17 Compare_forecast_assumptions.xlsx",row.names=F, sheetName = substocks[s], append=TRUE)
  }
}

# Compare forecast stock weights ------------------------------------------------------------

for (s in 1:3){
  # read in stock weights used in the forecasts
  dat <- read.xlsx("output/Forecast_stockwts.xlsx", sheetName = substocks[s])
  
  # plot
  png(paste0("report/plots/Fig 4.2", s-1, "b Compare_forecast_stockwts_",substocks[s],".png"),width = 11, height = 7, units = "in", res = 600)
  p1 <- ggplot(data=dat, aes(x=Age, y=Weight,colour=Source, group=Source)) +
    facet_wrap(~Year,nrow = 2)+
    geom_line() +geom_point(size=2) + theme_bw()+ labs(colour="",x="",y="Mean weight (kg)") +
    theme(axis.title=element_text(size=8),axis.text=element_text(size=8),
          legend.text=element_text(size=9)) + ggtitle(stock.name[s])
  print(p1)
  dev.off()
}

# Compare forecast selectivites ------------------------------------------------------------

for (s in 1:3){
  # read in selectivities used in forecast
  dat <- read.xlsx("output/Forecast_selectivity.xlsx", sheetName = substocks[s])
  
  # plot
  png(paste0("report/plots/Fig 4.2", s-1, "d Compare_forecast_selectivity_",substocks[s],".png"),width = 11, height = 7, units = "in", res = 600)
  p1 <- ggplot(data=dat, aes(x=Age, y=Selectivity,colour=Source, group=Source)) +
    facet_wrap(~Year,nrow = 2)+
    geom_line() + geom_point(size=2)+theme_bw()+ labs(colour="",x="",y="Selectivity") +
    theme(axis.title=element_text(size=8),axis.text=element_text(size=8),
          legend.text=element_text(size=9)) + ggtitle(stock.name[s])
  print(p1)
  dev.off()
}

# Compare forecast numbers at age ------------------------------------------------------------

for (s in 1:3){
  # read in N at age from the forecasts
  dat <- read.xlsx("output/Forecast_N_at_age.xlsx", sheetName = substocks[s])
  
  # plot
  png(paste0("report/plots/Fig 4.2", s-1, "a Compare_forecast_stock_numbers-at-age_",substocks[s],".png"),width = 11, height = 7, units = "in", res = 600)
  p1 <- ggplot(dat,aes(x=Age,y=N,group=interaction(Type,Source),colour=Source,shape=Type))+
    geom_line() + geom_point(size=3) + labs(colour="",y="Abundance (thousands)",shape="")+
    facet_wrap(~Year,nrow=2)+theme_bw()+scale_shape_manual(values=c(16, 2, 0))  + ggtitle(stock.name[s])
  print(p1)
  dev.off()
  
  # Read in N at age tables - as output from assessment
  dat.prev <- read.xlsx("output/Prev_assessment_N_at_age.xlsx", sheetName = substocks[s])
  dat.now <- read.xlsx("output/Now_assessment_N_at_age.xlsx", sheetName = substocks[s])
  colnames(dat.prev) <- colnames(dat.now) <- c("Year",ages)
  
  # Combine with forecast N at age
  # Forecast values ovewrite the assessment output here to account for when forecast values (e.g. median recruitment) are different to the assessment output
  tmp.prev <- dat[dat$Source %in% source.prev,]
  tmp.prev <- dcast(tmp.prev,Year~Age,value.var="N")
  dat.prev <- rbind(dat.prev[dat.prev$Year < min(tmp.prev$Year),],tmp.prev)
  
  tmp.now <- dat[dat$Source %in% source.now,]
  tmp.now <- dcast(tmp.now,Year~Age,value.var="N")
  dat.now <- rbind(dat.now[dat.now$Year < min(tmp.now$Year),],tmp.now)
  
  # Find ratio
  comp.yrs <- intersect(dat.now$Year,dat.prev$Year)
  rat.n <- dat.now[dat.now$Year %in% comp.yrs,] / dat.prev[dat.prev$Year %in% comp.yrs,]
  rat.n$Year <- comp.yrs
  
  rep <- data.frame(Age=c(as.character(1:6),"7+"),
                    round(t(rat.n[(nrow(rat.n)-9):nrow(rat.n),])[-1,],2))
  colnames(rep) <- c('Age', rat.n$Year[(nrow(rat.n)-9):nrow(rat.n)])
  
  # write out csv file
  if(s==1){
    write.xlsx(rat.n,file = "report/tables/Tab 4.16a Compare_forecast_ratio_stock_numbers-at-age.xlsx",row.names=F, sheetName = substocks[s])
    write.xlsx(rep, file = "report/tables/Tab 4.16a Compare_forecast_ratio_stock_numbers-at-age.xlsx",row.names=F, sheetName = paste0(substocks[s],"_report"), append=TRUE)
  }else{
    write.xlsx(rat.n,file = "report/tables/Tab 4.16a Compare_forecast_ratio_stock_numbers-at-age.xlsx",row.names=F, sheetName = substocks[s],append=TRUE)
    write.xlsx(rep, file = "report/tables/Tab 4.16a Compare_forecast_ratio_stock_numbers-at-age.xlsx",row.names=F, sheetName = paste0(substocks[s],"_report"), append=TRUE)
  }
}

# Compare forecast biomass at age ------------------------------------------------------------

for (s in 1:3){
  # read in biomass at age from the forecasts
  dat <- read.xlsx("output/Forecast_B_at_age.xlsx", sheetName = substocks[s])
  
  # plot
  png(paste0("report/plots/Fig 4.2", s-1, "c Compare_forecast_stock_biomass-at-age_",substocks[s],".png"),width = 11, height = 7, units = "in", res = 600)
  p1 <- ggplot(dat,aes(x=Age,y=B,group=interaction(Type,Source),colour=Source,shape=Type))+
    geom_line() + geom_point(size=3) + labs(colour="",y="Biomass (tonnes)",shape="")+
    facet_wrap(~Year,nrow=2)+theme_bw()+scale_shape_manual(values=c(16, 2, 0))  + ggtitle(stock.name[s])
  print(p1)
  dev.off()
  
  # Read in B at age tables - as output from assessment (N * stock weights)
  dat.prev <- read.xlsx("output/Prev_assessment_B_at_age.xlsx", sheetName = substocks[s])
  dat.now <- read.xlsx("output/Now_assessment_B_at_age.xlsx", sheetName = substocks[s])
  colnames(dat.prev) <- colnames(dat.now) <- c(ages,"Year")
  
  # Combine with forecast B at age
  # Forecast values ovewrite the assessment output here to account for when forecast values (e.g. median recruitment) are different to the assessment output
  tmp.prev <- dat[dat$Source %in% source.prev,]
  tmp.prev <- dcast(tmp.prev,Year~Age,value.var="B")
  dat.prev <- rbind(dat.prev[dat.prev$Year < min(tmp.prev$Year),],tmp.prev)
  
  tmp.now <- dat[dat$Source %in% source.now,]
  tmp.now <- dcast(tmp.now,Year~Age,value.var="B")
  dat.now <- rbind(dat.now[dat.now$Year < min(tmp.now$Year),],tmp.now)
  
  # Find ratio
  comp.yrs <- intersect(dat.now$Year,dat.prev$Year)
  rat.b <- dat.now[dat.now$Year %in% comp.yrs,] / dat.prev[dat.prev$Year %in% comp.yrs,]
  rat.b$Year <- comp.yrs
  
  rep <- data.frame(Age=c(as.character(1:6),"7+"),
                    round(t(rat.b[(nrow(rat.b)-9):nrow(rat.b),])[-8,],2))
  colnames(rep) <- c('Age', rat.b$Year[(nrow(rat.b)-9):nrow(rat.b)])
  
  # write out csv file
  if (s==1){
    write.xlsx(rat.b,file = "report/tables/Tab 4.16b Compare_forecast_ratio_stock_biomass-at-age.xlsx",row.names=F, sheetName = substocks[s])
    write.xlsx(rep, file = "report/tables/Tab 4.16b Compare_forecast_ratio_stock_biomass-at-age.xlsx",row.names=F, sheetName = paste0(substocks[s],"_report"), append=TRUE)
  }else{
    write.xlsx(rat.b,file = "report/tables/Tab 4.16b Compare_forecast_ratio_stock_biomass-at-age.xlsx",row.names=F, sheetName = substocks[s], append=TRUE)
    write.xlsx(rep, file = "report/tables/Tab 4.16b Compare_forecast_ratio_stock_biomass-at-age.xlsx",row.names=F, sheetName = paste0(substocks[s],"_report"), append=TRUE)
  }
}

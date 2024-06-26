
library(stockassessment)
library(multiStockassessment)
library(reshape2)
library(ggplot2)

## Load the multi stock fit
load("model/model.Rdata")

rp <- attr(fit,"m_rep")
pl <- attr(fit,"m_pl")

################################################################################

## Calculate

getLandings_SSF <- function(stock, season, fleet){
    ## NOTE: This assumes the seasons are quarters! It will only work for NS cod, not in general without modification for the time of seasons
    logN <- multiStockassessment:::splitParameter(pl$logN)[[stock]]
    CIF <- rp$mort.CIF_F_breakpoints[[stock]][,,fleet,season]
    if(season == 1){
        logSurvival <- 0
    }else{
        HazardPerSeason <- rp$mort.Hazard_breakpoints[[stock]][,,seq_len(season-1)]
        CumulativeHazardPerSeason <- HazardPerSeason * 0.25
        logSurvival <- - apply(CumulativeHazardPerSeason,1:2,sum)
    }
    CatchNum <- exp(logN) * exp(logSurvival) * CIF
    LF <- t(fit[[stock]]$data$landFrac[,,fleet])
    LW <- t(fit[[stock]]$data$landMeanWeight[,,fleet])
    LandingWeightAge <- CatchNum[,1:ncol(LF)] * LF * LW
    colSums(LandingWeightAge)
}

getLandings <- function(stock,season){
    Reduce("+",lapply(which(fit[[1]]$data$fleetTypes %in% c(0,7)), getLandings_SSF, stock=stock,season=season))
}

## Q1
names(fit)
NWQ1 <- getLandings(1,1)
SOQ1 <- getLandings(2,1)
VIQ1 <- getLandings(3,1)

## Q2-4
NWQ234 <- Reduce("+",lapply(2:4,function(s)getLandings(1,s)))
SOQ234 <- Reduce("+",lapply(2:4,function(s)getLandings(2,s)))
VIQ234 <- Reduce("+",lapply(2:4,function(s)getLandings(3,s)))

# Plots
Q1 <- data.frame(cbind(NWQ1,SOQ1,VIQ1), Year=names(NWQ1), Quarter="1")
colnames(Q1)[1:3] <- substr(colnames(Q1)[1:3], 1, 2)
Q1 <- melt(Q1, id.vars = c("Year","Quarter"))
Q234 <- data.frame(cbind(NWQ234,SOQ234,VIQ234), Year=names(NWQ234), Quarter="234")
colnames(Q234)[1:3] <- substr(colnames(Q234)[1:3], 1, 2)
Q234 <- melt(Q234, id.vars = c("Year","Quarter"))
plot <- rbind(Q1, Q234)

ggplot(plot[plot$Year>=1995,], aes(x = Year, y = value, fill = Quarter)) +
  geom_bar(stat = "identity", position = 'stack') +
  facet_grid(.~variable) +
  scale_y_continuous("Landings",
                     #labels = scales::percent,
                     expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("report/plots/Seasonal landings.png", width = 30, height=15, units='cm')

ggplot(plot[plot$Year>=1995,], aes(x = Year, y = value, fill = Quarter)) +
  geom_bar(stat = "identity", position = 'fill') +
  facet_grid(.~variable) +
  scale_y_continuous("Landings",
                     #labels = scales::percent,
                     expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("report/plots/Seasonal proportions.png", width = 30, height=15, units='cm')

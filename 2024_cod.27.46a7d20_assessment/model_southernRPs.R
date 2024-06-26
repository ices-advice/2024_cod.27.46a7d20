library(stockassessment)
library(msy)
library(FLCore)

mkdir("model/refpts/southern")

load("output/cod46a7d20_FLstockObjects_assessment_2024.Rdata")
load("model/refpts/sigmas.Rdata")

## extract subcomponents ----
southern <- stks[[2]]

## rename stocks ----
southern@name <- "Southern"

## Fbar ranges ---- 
range(southern, c("minfbar", "maxfbar")) <- c(2, 4) 

## Landings obligation adjustments ----
# Use wanted catch for ages 1-2: Fish below MCRS cannot be marketed
# Use all catch for ages 3+: Landing obligation
southern@landings.n[3:7,,,,]<-southern@landings.n[3:7,,,,]+southern@discards.n[3:7,,,,]
southern@landings.wt[3:7,,,,]<-southern@catch.wt[3:7,,,,]
southern@discards.n[3:7,,,,]<-0

southern@landings <- apply(southern@landings.n *  southern@landings.wt,2,sum)
southern@discards <- apply(southern@discards.n *  southern@discards.wt,2,sum)

## Flim function ----
FlimPredict<-function(SIM,BTarget){
  data.95<-SIM$rbp
  x.95 <- data.95[data.95$variable == "Spawning stock biomass", ]$Ftarget
  b.95 <- data.95[data.95$variable == "Spawning stock biomass", ]$p50
  par(mfrow = c(1, 1), mar = c(5, 4, 2, 1), mgp = c(3, 1, 0))
  b.lm <- loess(b.95 ~ x.95, span = 0.2)
  b.lm.pred <- data.frame(x = seq(min(x.95), max(x.95),length = 100000), y = rep(NA, 1000))
  b.lm.pred$y <- predict(b.lm, newdata = b.lm.pred$x)
  
  plot(b.lm.pred$x, b.lm.pred$y, type="l" ,ylim = c(0, max(b.95, na.rm = TRUE)),xlab = "Total catch F", ylab = "Median SSB")
  Flim<-approx(x=b.lm$fitted, y=b.lm$x, xout=BTarget)$y
  segments(0,BTarget,Flim,BTarget, col = "red", lty = 2)
  segments(Flim,0,Flim,BTarget, col = "red", lty = 2)
  text(.5*max(x.95), .5*max(b.95), paste("Flim=",round(Flim,2)), cex = 3)
  print(Flim)
}

## SR fit 1997 Type 2 ----------------------------------------------------------

set.seed(1)

southern_1997 <- window(southern, start = 1997, end = 2024)

FIT_1997_s_type2 <- eqsr_fit(southern_1997, nsamp = 2000, models = "Segreg")
eqsr_plot(FIT_1997_s_type2, n=2e4)
save(FIT_1997_s_type2, file = "model/refpts/southern/FIT_1997_s_type2.RData")

## EqSim 1997 Type 2 5yr -------------------------------------------------------

rm(list = ls.str(mode = 'numeric'))

## Variance params ----
Fcv <- 0.212
Fphi <- 0.423

sigma <- ifelse(sigmas[[2]]<0.2, 0.2, sigmas[[2]])

## Blim and Bpa ----
Blim <- FIT_1997_s_type2$sr.det$b
Bpa <- Blim * exp(1.645 * sigma)

### FMSY ----
set.seed(5)

SIM_1997_s_type2_FMSY_5yr <-
  eqsim_run(
    FIT_1997_s_type2,
    bio.years = c(2019, 2023), # 5 year biols
    bio.const = FALSE,
    sel.years = c(2021, 2023), # 3 year selects
    sel.const = FALSE,
    Fcv = Fcv,
    Fphi = Fphi,
    Btrigger = 0,
    Blim = Blim,
    Fscan = seq(0, 1.2, len = 61),
    verbose = FALSE
  )

eqsim_plot(SIM_1997_s_type2_FMSY_5yr, catch = "FALSE")

taf.png("model/refpts/southern/eqsim_Fmsy.png")
eqsim_plot_range(SIM_1997_s_type2_FMSY_5yr, type="median")
dev.off()

eqsim_plot_range(SIM_1997_s_type2_FMSY_5yr, type='ssb')
SIM_1997_s_type2_FMSY_5yr$Refs2

save(SIM_1997_s_type2_FMSY_5yr, file = "model/refpts/southern/SIM_1997_s_type2_FMSY_5yr.RData")

## Check Fbar below Fmsy ----
tail(fbar(window(southern_1997, end = 2023)),5) <= SIM_1997_s_type2_FMSY_5yr$Refs2["lanF", "medianMSY"] ## CHECK REFS!!

## Bpa and Btrig ----
Btrig <- Bpa 

### FP0.5 ----
set.seed(6)

SIM_1997_s_type2_Fp.05_5yr <-
  eqsim_run(
    FIT_1997_s_type2,
    bio.years = c(2019, 2023), # 5 year biols
    bio.const = FALSE,
    sel.years = c(2021, 2023), # 3 year selects
    sel.const = FALSE,
    Fcv = Fcv,
    Fphi = Fphi,
    Btrigger = Btrig,
    Blim = Blim,
    Fscan = seq(0, 1.2, len = 61),
    verbose = FALSE
  )
eqsim_plot(SIM_1997_s_type2_Fp.05_5yr, catch = "FALSE")

taf.png("model/refpts/southern/eqsim_Fpa.png")
eqsim_plot_range(SIM_1997_s_type2_Fp.05_5yr, type="median")
dev.off()

save(SIM_1997_s_type2_Fp.05_5yr, file = "model/refpts/southern/SIM_1997_s_type2_Fp.05_5yr.RData")

## Flim ----
set.seed(7)

SIM_1997_s_type2_5yr_Flim <-
  eqsim_run(
    FIT_1997_s_type2,
    bio.years = c(2019, 2023),
    bio.const = FALSE,
    sel.years = c(2021, 2023),
    sel.const = FALSE,
    Fcv = 0,
    Fphi = 0,
    Btrigger = 0,
    Blim = Blim,
    Fscan = seq(0, 1.2, len = 61),
    verbose = FALSE
  )

save(SIM_1997_s_type2_5yr_Flim, file = "model/refpts/southern/SIM_1997_s_type2_5yr_Flim.RData")

Flim <- FlimPredict(SIM_1997_s_type2_5yr_Flim, Blim)

## Final refPoints ----
S_1997_type2_5yr_refPoints <- list(MSYbtrigger = Btrig,
                                    Fmsy = SIM_1997_s_type2_FMSY_5yr$Refs2["lanF", "medianMSY"],
                                    Fmsyupper = SIM_1997_s_type2_FMSY_5yr$Refs2["lanF", "Medupper"],
                                    Fmsylower = SIM_1997_s_type2_FMSY_5yr$Refs2["lanF", "Medlower"],
                                    Blim = Blim,
                                    Bpa = Bpa,
                                    Flim = Flim,
                                    Fp.05 = SIM_1997_s_type2_Fp.05_5yr$Refs2["catF", "F05"],
                                    Fpa = SIM_1997_s_type2_Fp.05_5yr$Refs2["catF", "F05"],
                                    biol.years = dim(SIM_1997_s_type2_FMSY_5yr$ibya$weca)[2],
                                    sel.years = dim(SIM_1997_s_type2_FMSY_5yr$ibya$sel)[2])

S_1997_type2_5yr_refPoints <- as.data.frame(S_1997_type2_5yr_refPoints)
save(S_1997_type2_5yr_refPoints, file = "model/refpts/southern/S_1997_type2_5yr_refPoints.RData")

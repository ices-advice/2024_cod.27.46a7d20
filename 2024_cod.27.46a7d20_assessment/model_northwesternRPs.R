library(stockassessment)
library(msy)
library(FLCore)

mkdir("model/refpts/northwestern")

load("output/cod46a7d20_FLstockObjects_assessment_2024.Rdata")
load("model/refpts/sigmas.Rdata")

## extract subcomponents ----
northwestern <- stks[[1]]

## rename stocks ----
northwestern@name <- "Northwestern"

## Fbar ranges ----
range(northwestern, c("minfbar", "maxfbar")) <- c(2, 4) 

## Landings obligation adjustments ----
# Use wanted catch for ages 1-2: Fish below MCRS cannot be marketed
# Use all catch for ages 3+: Landing obligation
northwestern@landings.n[3:7,,,,]<-northwestern@landings.n[3:7,,,,]+northwestern@discards.n[3:7,,,,]
northwestern@landings.wt[3:7,,,,]<-northwestern@catch.wt[3:7,,,,]
northwestern@discards.n[3:7,,,,]<-0

northwestern@landings <- apply(northwestern@landings.n *  northwestern@landings.wt,2,sum)
northwestern@discards <- apply(northwestern@discards.n *  northwestern@discards.wt,2,sum)

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


## SR fit 1997 -----------------------------------------------------------------

northwestern_1997 <- window(northwestern, start = 1997, end = 2023)

set.seed(5)

## Type 2: estimate breakpoint of segemented regression internally in EqSim ----
FIT_1997_nw_type2 <- eqsr_fit(northwestern_1997, nsamp = 2000, models = "Segreg")
save(FIT_1997_nw_type2, file = "model/refpts/northwestern/FIT_1997_nw_type2.RData")
eqsr_plot(FIT_1997_nw_type2, n=2e4)


## EqSim Type 2 1997 -----------------------------------------------------------

rm(list = ls.str(mode = 'numeric'))

## Variance params ----
Fcv <- 0.212
Fphi <- 0.423

sigma <- ifelse(sigmas[[1]]<0.2, 0.2, sigmas[[1]])

## Blim ----

Blim <- FIT_1997_nw_type2$sr.det$b

## FMSY ----
set.seed(6)

SIM_1997_nw_type2_FMSY <-
  eqsim_run(
    FIT_1997_nw_type2,
    bio.years = c(2014, 2023),
    bio.const = FALSE,
    sel.years = c(2021, 2023),
    sel.const = FALSE,
    Fcv = Fcv,
    Fphi = Fphi,
    Btrigger = 0,
    Blim = Blim,
    Fscan = seq(0, 1.2, len = 61),
    verbose = FALSE
  )

eqsim_plot(SIM_1997_nw_type2_FMSY, catch = "FALSE")

taf.png("model/refpts/northwestern/eqsim_Fmsy.png")
eqsim_plot_range(SIM_1997_nw_type2_FMSY, type="median")
dev.off()

eqsim_plot_range(SIM_1997_nw_type2_FMSY, type='ssb')
SIM_1997_nw_type2_FMSY$Refs2

save(SIM_1997_nw_type2_FMSY, file = "model/refpts/northwestern/SIM_1997_nw_type2_FMSY.RData")

## Check Fbar against Fmsy ----
tail(fbar(window(northwestern, end = 2023)),5) <= SIM_1997_nw_type2_FMSY$Refs2["lanF", "medianMSY"] ## CHECK REFS!!

## Bpa and Btrigger ----
Bpa <- Blim * exp(1.645 * sigma)
Btrig <- Bpa 

## FP0.5 ----
set.seed(7)

SIM_1997_nw_type2_Fp.05 <-
  eqsim_run(
    FIT_1997_nw_type2,
    bio.years = c(2014, 2023),
    bio.const = FALSE,
    sel.years = c(2021, 2023),
    sel.const = FALSE,
    Fcv = Fcv,
    Fphi = Fphi,
    Btrigger = Btrig,
    Blim = Blim,
    Fscan = seq(0, 1.2, len = 61),
    verbose = FALSE
  )
eqsim_plot(SIM_1997_nw_type2_Fp.05, catch = "FALSE")

taf.png("model/refpts/northwestern/eqsim_Fpa.png")
eqsim_plot_range(SIM_1997_nw_type2_Fp.05, type="median")
dev.off()

save(SIM_1997_nw_type2_Fp.05, file = "model/refpts/northwestern/SIM_1997_nw_type2_Fp.05.RData")

## Flim ----
set.seed(8)

SIM_1997_nw_type2_Flim <-
  eqsim_run(
    FIT_1997_nw_type2,
    bio.years = c(2014, 2023),
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

save(SIM_1997_nw_type2_Flim, file = "model/refpts/northwestern/SIM_1997_nw_type2_Flim.RData")

Flim <- FlimPredict(SIM_1997_nw_type2_Flim, Blim)

## Final refpoints ----
NW_1997_type2_refPoints <- list(MSYbtrigger = Btrig,
                               Fmsy = SIM_1997_nw_type2_FMSY$Refs2["lanF", "medianMSY"],
                               Fmsyupper = SIM_1997_nw_type2_FMSY$Refs2["lanF", "Medupper"],
                               Fmsylower = SIM_1997_nw_type2_FMSY$Refs2["lanF", "Medlower"],
                               Blim = Blim,
                               Bpa = Bpa,
                               Flim = Flim,
                               Fp.05 = SIM_1997_nw_type2_Fp.05$Refs2["catF", "F05"],
                               Fpa = SIM_1997_nw_type2_Fp.05$Refs2["catF", "F05"],
                               biol.years = dim(SIM_1997_nw_type2_FMSY$ibya$weca)[2],
                               sel.years = dim(SIM_1997_nw_type2_FMSY$ibya$sel)[2])

NW_1997_type2_refPoints <- as.data.frame(NW_1997_type2_refPoints)

save(NW_1997_type2_refPoints, file = "model/refpts/northwestern/NW_1997_type2_refPoints.RData")

library(stockassessment)
library(msy)
library(FLCore)

mkdir("model/refpts/viking")

load("output/cod46a7d20_FLstockObjects_assessment_2024.Rdata")
load("model/refpts/sigmas.Rdata")

## extract subcomponents ----
viking <- stks[[3]]

## Fbar ranges ---- 
range(viking, c("minfbar", "maxfbar")) <- c(2, 4)

## Landings obligation adjustments ----
# Use wanted catch for ages 1-2: Fish below MCRS cannot be marketed
# Use all catch for ages 3+: Landing obligation
viking@landings.n[3:7,,,,]<-viking@landings.n[3:7,,,,]+viking@discards.n[3:7,,,,]
viking@landings.wt[3:7,,,,]<-viking@catch.wt[3:7,,,,]
viking@discards.n[3:7,,,,]<-0

viking@landings <- apply(viking@landings.n *  viking@landings.wt,2,sum)
viking@discards <- apply(viking@discards.n *  viking@discards.wt,2,sum)

## FLSR objects ----
vsr <- as.FLSR(viking)

## Fix SegReg breakpoint ----
SegregFix <- function (ab, ssb) {
  log(ifelse (ssb>=B, ab$a*B, ab$a*ssb))
}

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

## SR fit 1997 Type M ----------------------------------------------------------

rm(list = ls.str(mode = 'numeric'))

viking_1997 <- window(viking, start = 1997, end = 2024)

vsr_1997 <- as.FLSR(viking_1997)
data <- model.frame(FLQuants(vsr_1997, c("ssb", "rec")))
med <- median(data$ssb)
lessMed <- data[data$ssb<med,]
recAvg <- mean(data$rec)
big <- lessMed$ssb[lessMed$rec>recAvg]
Blim_v_1997 <- mean(big)

B <- Blim_v_1997

plot(an(ssb(vsr_1997))[1:26], an(rec(vsr_1997))[1:26], xlim = c(0, max(an(ssb(vsr_1997)))), ylim = c(0, max(an(rec(vsr_1997)))),
     xlab="SSB", ylab="Recruitment", bty = "l", cex = 1.5)
points(an(ssb(vsr_1997))[27], an(rec(vsr_1997))[27], col='red', bty = "l", cex = 1.5)
abline(v = B, col='blue')
abline(v = med)
abline(h = recAvg)

set.seed(5)
FIT_1997_v_typeM <- eqsr_fit(viking_1997, nsamp = 2000, models = "SegregFix")
eqsr_plot(FIT_1997_v_typeM, n=2e4)

## EqSim 1997 Type M -----------------------------------------------------------

# rm(list = ls.str(mode = 'numeric'))

## Variance parameters ----
Fcv <- 0.212
Fphi <- 0.423

sigma <- ifelse(sigmas[[3]]<0.2, 0.2, sigmas[[3]])

## Blim ----
Blim <- Blim_v_1997

### FMSY ----
set.seed(6)

SIM_1997_v_typeM_FMSY <-
  eqsim_run(
    FIT_1997_v_typeM,
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

eqsim_plot(SIM_1997_v_typeM_FMSY, catch = "FALSE")

taf.png("model/refpts/viking/eqsim_Fmsy.png")
eqsim_plot_range(SIM_1997_v_typeM_FMSY, type="median")
dev.off()

eqsim_plot_range(SIM_1997_v_typeM_FMSY, type='ssb')
SIM_1997_v_typeM_FMSY$Refs2

save(SIM_1997_v_typeM_FMSY, file = "model/refpts/viking/SIM_1997_v_typeM_FMSY.RData")

## Check Fbar against Fmsy ----
tail(fbar(window(viking_1997, end = 2023)),5) <= SIM_1997_v_typeM_FMSY$Refs2["lanF", "medianMSY"] ## CHECK REFS!!

## Bpa and Btrig ----
Bpa <- Blim * exp(1.645 * sigma)
Btrig <- Bpa 

### FP0.5 ----
set.seed(7)

SIM_1997_v_typeM_Fp.05 <-
  eqsim_run(
    FIT_1997_v_typeM,
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
eqsim_plot(SIM_1997_v_typeM_Fp.05, catch = "FALSE")

taf.png("model/refpts/viking/eqsim_Fpa.png")
eqsim_plot_range(SIM_1997_v_typeM_Fp.05, type="median")
dev.off()

save(SIM_1997_v_typeM_Fp.05, file = "model/refpts/viking/SIM_1997_v_typeM_Fp.05.RData")

## Flim ----
set.seed(8)

SIM_1997_v_typeM_Flim <-
  eqsim_run(
    FIT_1997_v_typeM,
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

save(SIM_1997_v_typeM_Flim, file = "model/refpts/viking/SIM_1997_v_typeM_Flim.RData")

Flim <- FlimPredict(SIM_1997_v_typeM_Flim, Blim)

## Ref points ----
V_1997_typeM_refPoints <- list(MSYbtrigger = Btrig,
                              Fmsy = SIM_1997_v_typeM_FMSY$Refs2["lanF", "medianMSY"],
                              Fmsyupper = SIM_1997_v_typeM_FMSY$Refs2["lanF", "Medupper"],
                              Fmsylower = SIM_1997_v_typeM_FMSY$Refs2["lanF", "Medlower"],
                              Blim = Blim,
                              Bpa = Bpa,
                              Flim = Flim,
                              Fp.05 = SIM_1997_v_typeM_Fp.05$Refs2["catF", "F05"],
                              Fpa = SIM_1997_v_typeM_Fp.05$Refs2["catF", "F05"],
                              biol.years = dim(SIM_1997_v_typeM_FMSY$ibya$weca)[2],
                              sel.years = dim(SIM_1997_v_typeM_FMSY$ibya$sel)[2])

V_1997_typeM_refPoints <- as.data.frame(V_1997_typeM_refPoints)
save(V_1997_typeM_refPoints, file = "model/refpts/viking/V_1997_typeM_refPoints.RData")

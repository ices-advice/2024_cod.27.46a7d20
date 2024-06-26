## Extract results of interest, write TAF output tables

## Before: sw.Rdata, maturity.Rdata, Q1models.Rdata, Q34models.Rdata, M1M2.txt, No_magic_M.csv
## After: assessment input files

library(icesTAF)
library(surveyIndex)

mkdir("output")

load("model/sw.Rdata")
load("model/maturity.Rdata")
load("model/Q1models.Rdata")
models1 <- models
load("model/Q34models.Rdata")
models34 <- models
M_2020 <- read.table(file = "boot/data/M1M2.txt", header = TRUE)
M_2023 <- read.csv(file = 'boot/data/No_Magic_M.csv')

## 1 Stock weights
for (a in unique(WAA.mean$SubArea)){
  swt <-subset(WAA.mean,SubArea==a)
  #output
  outfile=paste0("output/sw_",a,".dat")
  cat(paste0(a," substock weights bootstrapped mean. L-Wt params: 2011-2022. Wts updated: 1983 onwards"),"\n",file=outfile)
  cat(c(1,2),sep="\t","\n",file=outfile,append=TRUE)
  cat(c(min(swt$Year),max(swt$Year)),sep="\t","\n",file=outfile,append=TRUE)
  cat(c(1,7),sep="\t","\n",file=outfile,append=TRUE)
  cat(1,"\n",file=outfile,append=TRUE)
  write.table(swt[,as.character(1:7)],file=outfile,row.names=FALSE,col.names=FALSE,append=TRUE)
  write.taf(swt[,as.character(1:7)], file=paste0("output/sw_",a,".csv"))
}

## 2 Maturity
for (i in 1:length(mo)){
  file=paste0("output/mo_", names(mo)[i], ".dat")
  cat(paste0(names(mo)[i], " cod maturity"),"\n",file=file)
  cat(c(1,6),sep="\t","\n",file=file,append=TRUE)
  cat(c(1983,(1983+nrow(mo[[i]])-1)),sep="\t","\n",file=file,append=TRUE)
  cat(c(1,(1+ncol(mo[[i]])-1)),sep="\t","\n",file=file,append=TRUE)
  cat(1,"\n",file=file,append=TRUE)
  write.table(mo[[i]],file=file,row.names=FALSE,col.names=FALSE,append=TRUE)
  write.taf(mo[[i]], file=paste0("output/mat_", names(mo)[i], ".csv"))
}

## 3 Indices

# Q1
indices1 = lapply(models1,function(x) x$idx )
los1 = lapply(models1,function(x) x$lo )
ups1 = lapply(models1,function(x) x$up )

for (i in 2:4){
  # Indices
  exportSI(models1[[i]]$idx,as.numeric(colnames(models1[[i]]$idx)),as.numeric(rownames(models1[[i]]$idx)),toy=0.25/2,file=paste0("output/", names(models1)[i], "Q1.txt"),nam=paste0(names(models1)[i], " cod indices Q1"))
  
  # sds
  sdQ1 = (log(ups1[[i]]) - log(los1[[i]]))/4
  rownames(sdQ1) = rownames(models1[[i]]$idx)
  colnames(sdQ1) = colnames(models1[[i]]$idx)
  write.table(sdQ1,file=paste0("output/", names(models1)[i], "Q1sds.txt"))
}

# Q34
indices34 = lapply(models34,function(x) x$idx )
los34 = lapply(models34,function(x) x$lo )
ups34 = lapply(models34,function(x) x$up )

# Index
exportSI(models34[[1]]$idx[,-1],as.numeric(colnames(models34[[1]]$idx)[-1]),as.numeric(rownames(models34[[1]]$idx)),toy=1.5/2,file="output/q34index.txt",nam="NS cod indices Q3 & Q4")

# Recruitment indices
for (i in 2:4){
  file=paste0("output/", names(models34)[i], "34recruits.txt")
  cat(paste0(names(models34)[i], " recruitment index"), "\n", file = file)
  cat(range(as.numeric(as.character(as.numeric(rownames(models34[[i]]$idx))+1))), "\n", file = file, 
      append = TRUE)
  cat("1 1 ", rep(0, 2), "\n", file = file, append = TRUE)
  cat(c(1,1), "\n", file = file, append = TRUE)
  write.table(round(cbind(1, models34[[i]]$idx[,1]), 4), file = file, row.names = FALSE, 
              col.names = FALSE, append = TRUE)
}

# Index sds
sdSIQ34 = (log(ups34[[1]][,-1]) - log(los34[[1]][,-1]))/4
rownames(sdSIQ34) = rownames(models34[[1]]$idx)
colnames(sdSIQ34) = colnames(models34[[1]]$idx)[-1]
write.table(sdSIQ34,file="output/q34sds.txt")

# Recruitment sds
for (i in 2:4){
  sdR = as.data.frame((log(ups34[[i]][,1]) - log(los34[[i]][,1]))/4)
  rownames(sdR) = as.character(as.numeric(rownames(models[[i]]$idx))+1)
  colnames(sdR) = colnames(models[[i]]$idx)[1]
  write.table(sdR,file=paste0("output/", names(models34)[i], "34Rsds.txt"))
}

## 4 Natural mortality

# 2020 Key run
row.names(M_2020) <- M_2020[,1]
M_2020 <- M_2020[,-(1:2)] # We don't use age 0
M_2020 <- rbind(as.matrix(M_2020[-(1:9),]),
           matrix(NA, nrow = length(2020:2024), ncol = ncol(M_2020)))
M_2020 <- cbind(M_2020, M_2020[,rep(ncol(M_2020), length(11:15))])

file="output/nm_2020.dat"
cat("Natural Mortality 2020 key run","\n",file=file)
cat(c(1,5),sep="\t","\n",file=file,append=TRUE)
cat(c(1983,(1983+nrow(M_2020)-1)),sep="\t","\n",file=file,append=TRUE)
cat(c(1,15),sep="\t","\n",file=file,append=TRUE)
cat(1,"\n",file=file,append=TRUE)
write.table(M_2020,file=file,row.names=FALSE,col.names=FALSE,append=TRUE)

# 2023 cod run
row.names(M_2023) <- M_2023[,1]
M_2023 <- M_2023[,-(1:2)] # We don't use age 0
M_2023 <- rbind(as.matrix(M_2023[-(1:9),]),
                matrix(NA, nrow = length(2023:2024), ncol = ncol(M_2023)))
M_2023 <- cbind(M_2023, M_2023[,rep(ncol(M_2023), length(11:15))])

file="output/nm_2023.dat"
cat("Natural Mortality 2023 key run for cod","\n",file=file)
cat(c(1,5),sep="\t","\n",file=file,append=TRUE)
cat(c(1983,(1983+nrow(M_2023)-1)),sep="\t","\n",file=file,append=TRUE)
cat(c(1,15),sep="\t","\n",file=file,append=TRUE)
cat(1,"\n",file=file,append=TRUE)
write.table(M_2023,file=file,row.names=FALSE,col.names=FALSE,append=TRUE)

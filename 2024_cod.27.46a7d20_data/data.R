## Preprocess data, write TAF data tables

## Before: codData050424.Rdata & associated files
## After: Q1preprocessed.Rdata, Q34preprocessed.Rdata, mat_data.Rdata

library(icesTAF)

mkdir("data")

source.taf("data_prepareSurveyQ1.R")
source.taf("data_prepareSurveyQ34.R")
source.taf("data_prepareMaturity.R")
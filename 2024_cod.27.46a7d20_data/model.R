## Run analysis, write model results

## Before: Q1preprocessed.Rdata, Q34preprocessed.Rdata, mat_data.Rdata
## After: sw.Rdata, maturity.Rdata, Q1models.Rdata, Q1retro.Rdata, Q34models.Rdata, Q34retro.Rdata

library(icesTAF)

mkdir("model")

source.taf("model_stockWeights.R")
source.taf("model_maturity.R")

# Takes a long time to run.
# I switch to WSL and run from the scripts.
# source.taf("model_indicesQ1.R")
# source.taf("model_indicesQ34.R")
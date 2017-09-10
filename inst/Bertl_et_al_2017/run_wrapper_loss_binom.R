#### settings ####

settings = commandArgs(trailingOnly=T)
cv = as.numeric(settings[1]) # how many CV slices?
cv_index = as.numeric(settings[2]) # which CV slice?
mi_index = as.numeric(settings[3]) # which MI slice?
model_index = as.numeric(settings[4]) # which model?
modelfile = settings[5] # model file
datafolder = settings[6] # data folder
resultsfolder = settings[7] # results folder
per_obs = as.logical(settings[8]) # results folder
nested_samples = as.logical(settings[9]) # results folder

print(nested_samples)

# packages
.libPaths("/project/PCAWG/nullmodel/R_packages")
require("multinomutils")

wrapper_loss_binom(cv, cv_index, mi_index, model_index, modelfile, datafolder, resultsfolder, per_obs=per_obs, nested_samples = nested_samples)

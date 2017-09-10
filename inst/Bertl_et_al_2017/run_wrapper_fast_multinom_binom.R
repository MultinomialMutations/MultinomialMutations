#### settings ####

# settings from the command line
settings = commandArgs(trailingOnly=T)
model.index = as.numeric(settings[1]) # number of the model in the model file
mi = as.numeric(settings[2]) # mi replicate
datafolder = settings[3] # data folder
modelfile = settings[4] # model file
resultsfolder = settings[5] # results folder
contrasts = settings[6]
nested_samples = as.logical(settings[7])
VC = as.logical(settings[8])
loglik = as.logical(settings[9])
startfile = settings[10] # starting values (can be NULL)
if(startfile=="NULL") startfile=NULL
libloc = settings[11] # location of packages (fastmultinom)

# packages
.libPaths(libloc)
require("multinomutils")

wrapper_fast_multinom_binom(model.index, mi, datafolder, modelfile, resultsfolder, startfile, contrasts, nested_samples, VC, loglik)

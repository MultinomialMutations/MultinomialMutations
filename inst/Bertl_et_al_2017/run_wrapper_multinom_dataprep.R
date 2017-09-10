settings = commandArgs(trailingOnly=T)
cv = as.numeric(settings[1]) # which CV slice?
m = as.numeric(settings[2]) # how-many fold imputation?
impute = as.logical(settings[3])
datafolder = settings[4] # data folder
datafile = settings[5] # data file (incl. number, excl. ".txt")
impfolder = settings[6] # subfolder for imputed datasets
data_source = settings[7] # fredriksson or pcawg
libfolder = settings[8] # folder for packages

.libPaths(libfolder)
library("multinomutils")

wrapper_multinom_dataprep(cv, m, impute, datafolder, datafile, impfolder, data_source)

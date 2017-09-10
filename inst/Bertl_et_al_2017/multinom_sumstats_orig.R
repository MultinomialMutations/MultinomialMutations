settings = commandArgs(trailingOnly=T)
cv = as.numeric(settings[1])
datafolder = settings[2] # data folder
datafile = settings[3] # data file 
libfolder = settings[4] # folder for packages

.libPaths(libfolder)
library("multinomutils")

sumstats_orig(inputfile=paste0(datafolder, datafile, cv-1), outputfile=paste0(datafolder, datafile, cv-1, "_sumstats.Rdata"), explanatory=c("genomicSeg", "replication_timing", "cancer_type", "LocalMutRate4"), no_pdf_output=T)
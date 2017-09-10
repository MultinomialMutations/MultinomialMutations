### NOTE! This code can produce wrong results when the model file contains more models than used here! ################################################################################

#### settings ####

settings = commandArgs(trailingOnly=T)
cv = as.numeric(settings[1]) # how many CV slices?
mi = as.numeric(settings[2]) # how many MI slices?
mo = as.numeric(settings[3]) # how many models?
modelfile = settings[4] # model file
resultsfolder=settings[5] # results folder
outputfile = settings[6] # name for output files (tex and csv)

load(modelfile)
require("xtable")

# read loss values and compute mean and standard error for each model
lossmatrix = matrix(NA, nrow=cv*mi, ncol=mo)

for(mi_index in 1:mi){
  for(cv_index in 1:cv){
    for(mo_index in 1:mo){
      est.loss = NA
      try(load(file=paste0(resultsfolder, "loss_mi", mi_index, "cv", cv_index, "model", mo_index, ".Rdata")), silent=T) # -> est.loss
      lossmatrix[cv_index + (mi_index-1)*cv,mo_index] = est.loss
      
    }
  }
}

meanloss = colMeans(lossmatrix, na.rm=T)
seloss = apply(lossmatrix, 2, sd, na.rm=T)


#### output matrix ####

# csv:
modelmat = cbind(models, mean = meanloss, se = seloss) 
write.csv(modelmat, file=paste0(resultsfolder, outputfile, ".csv"))

# latex table: 
modelmat[modelmat==0]=NA # for xtable: NAs will not be printed
alignvec = paste0("r|", paste(rep("c", ncol(models)), collapse=""), "|rr")
displayvec = c(rep("d", ncol(models)+1), "f", "f")
print(xtable(modelmat, display=displayvec, align=alignvec), file=paste0(resultsfolder, outputfile, ".tex"))

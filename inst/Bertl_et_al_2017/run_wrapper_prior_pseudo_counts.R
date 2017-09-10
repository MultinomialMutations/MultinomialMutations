### settings ####

settings = commandArgs(trailingOnly=T)
cv = as.numeric(settings[1]) # which CV partition?
mi = as.numeric(settings[2]) # which MI replicate?
datafolder = settings[3] # data folder
dirichletfolder = settings[4] # where to save the file with the dirichlet prior
rm_folder = settings[5]
rm_cancer = as.logical(settings[6]) # which cancer types to remove (or NULL)
libfolder = settings[7] # folder that contains packages

.libPaths(libfolder)
library("multinomutils")

### adding a Dirichlet prior ###

if(rm_cancer) rm_cancer="c('KICH', 'LGG', 'PRAD', 'THCA')"

wrapper_prior_pseudo_counts(cv, mi, datafolder, dirichletfolder, rm_folder, rm_cancer)

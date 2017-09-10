#!/bin/sh
# to be run on GenomeDK using bash

# Multiple imputation and computation of the categories for multinomial regression and some explanatory variables
# For the CV slices and for the joint dataset (created here)



### settings: ###

# jobname
jobname="PCAWG"

# data folders and files:
datafolder=/project/PCAWG/faststorage/count_final/
impfolder=imputed/ # subfolder for the prepared (and imputed) datasets 
datafile=set # CV slices, numbering starting from 0
datafileall=setall # name of the joint dataset
impfolderall=imputed_all/ # subfolder for the prepared (and imputed) joint datasets 

# package locations
libfolder=/project/PCAWG/nullmodel/R_packages/ # where to find R-packages

# number of CV partitions:
cv=1
# number of MIs:
mi=0
# imputation?
impute=F
cv_index=1


### preparations ###

# load R
source /com/extra/R/3.1.0/load.sh

# make folders
mkdir -p ${datafolder}/${impfolder}
mkdir -p ${datafolder}/${impfolderall}
mkdir -p ${datafolder}/dataprep_log


### summary statistics on the original file
qx -v -p -t 01:00:00 -m 32g --queue=normal -i=*.R "Rscript multinom_sumstats_orig.R 1 $datafolder $datafileall $libfolder > ${datafolder}/dataprep_log/dataprep_${jobname}_sumstats_${cv_index}.Rout " | sbatch --job-name=sumstats_${jobname}_cv${cv_index} --out=${datafolder}/dataprep_log/sumstats_${jobname}_cv${cv_index}.out --mail-user=guo@cs.au.dk --mail-type=END -A PCAWG

# creation of multiple imputation replicates and categories

qx -v -p -t 3:00:00 -m 64g --queue=normal -i=*.R "Rscript run_wrapper_multinom_dataprep.R 1 $mi $impute $datafolder $datafileall $impfolderall 'fredriksson' $libfolder > ${datafolder}/dataprep_log/dataprep_${jobname}_all.Rout " | sbatch --job-name=dataprep_${jobname}_all --out=${datafolder}/dataprep_log/dataprep_${jobname}_all.out --mail-user=guo@cs.au.dk --mail-type=END -A PCAWG


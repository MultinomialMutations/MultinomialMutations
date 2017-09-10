#!/bin/sh

### estimating a set of models ###
# for the goal of parameter estimation (not model selection).
# The package fastmultinom (by Johanna) is used.

# to be run on genomeDK using bash


### settings: ###


# data folder:
datafolder=/project/PCAWG/faststorage/count_final/imputed_dirichlet_all_rm/1/
# results folder:
resultsfolder=/project/PCAWG/faststorage/count_final/model/
# file that contains all models:
modelfile=/home/guo/PCAWG/faststorage/count_final/final_model_without_sample.RData
# contrasts
contrasts=all_sum
# sample_id nested in cancer type
nested_samples=T
# compute VC matrix?
VC=T
# compute log likelihood?
loglik=T
# file that contains the starting values:
startfile=NULL
# folder that contains the R-packages that are used
libloc=/project/PCAWG/nullmodel/R_packages/
# jobname
jobname=final_all


# number of models:
mo=2
# number of MIs:
mi=1


### preparations ###

# make resultsfolder
mkdir -p $resultsfolder
# make subfolder for log files
mkdir -p $resultsfolder/logfiles

# load R
source /com/extra/R/3.1.0/load.sh



### multinomial model: estimation on all mi replicates: ###

for (( mo_index = 1; mo_index <= $mo; mo_index++ ))
do
	for (( mi_index = 1; mi_index <= $mi; mi_index++ ))
	do
		qx -v -p -t 02:00:00 -m 64g --queue=normal -i=*.R "Rscript run_wrapper_fast_multinom.R $mo_index $mi_index $datafolder $modelfile $resultsfolder $contrasts $nested_samples $VC $loglik $startfile $libloc > ${resultsfolder}logfiles/fast_multinom_${jobname}_mo${mo_index}_mi${mi_index}.Rout " | sbatch --job-name=final_${jobname}_${mo_index}_${mi_index} --out=${resultsfolder}/logfiles/fast_multinom_${jobname}_mo${mo_index}_mi${mi_index}.out --err=${resultsfolder}logfiles/fast_multinom_${jobname}_mo${mo_index}_mi${mi_index}.err -A PCAWG
	done
done 

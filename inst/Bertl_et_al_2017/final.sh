#!/bin/sh

### estimating a set of models ###
# for the goal of estimation (not model selection).
# The package fastmultinom (by Johanna) is used.

# to be run on genomeDK using bash


### settings: ###


# data folder:
datafolder=/project/PCAWG/faststorage/count_final/imputed_dirichlet_all/1/
# results folder:
resultsfolder=/project/PCAWG/faststorage/count_final/model/
# file that contains all models:
#####
modelfile=/project/PCAWG/faststorage/count_final/final_model.Rdata
# contrasts
contrasts=treatment
# sample_id nested in cancer type
nested_samples=F
# compute VC matrix?
VC=T
# compute log likelihood?
loglik=T
# file that contains the starting values:
startfile=NULL
# folder that contains the R-packages that are used
libloc=/project/PCAWG/nullmodel/R_packages/
# jobname
jobname=pcawg
# cp (guo)
cp /project/PCAWG/faststorage/count_data_Fredriksson_whole_genome_Jan2016/old/imputed_dirichlet_all/1/sample_id.txt /project/PCAWG/faststorage/count_final/imputed_dirichlet_all/1/sample_id.txt
# list of sample_ids
samplelist=/project/PCAWG/faststorage/count_final/imputed_dirichlet_all/1/sample_id.txt
#samplelist=/project/PCAWG/faststorage/count_data_Fredriksson_whole_genome_Jan2016/imputed_dirichlet_all/1/s.txt


# number of model:
mo=1
# number of MI:
mi=1
# (there is no loop over the models and MIs here!)


### preparations ###

# make resultsfolder
mkdir -p $resultsfolder
# make subfolder for log files
mkdir -p $resultsfolder/logfiles

# load R
source /com/extra/R/3.1.0/load.sh


### multinomial model: estimation on all mi replicates: ###

samplenum=0
# open the file to read parameter values and names from:
exec 3< $samplelist

while read -u 3 sample  # read file line by line, in each one extracting this variable
do 
	mkdir -p ${resultsfolder}${sample}	
	samplenum=$((samplenum+1))

	qx -v -p -t 01:00:00 -m 4g --queue=normal -i=*.R "Rscript run_wrapper_fast_multinom.R $mo $mi ${datafolder}${sample}/ $modelfile $resultsfolder${sample}/ $contrasts $nested_samples $VC $loglik $startfile $libloc > ${resultsfolder}logfiles/final_$jobname_mo${mo}_mi${mi}_${samplenum}.Rout " | sbatch --job-name=final_${jobname}_${samplenum} --out=${resultsfolder}logfiles/fast_multinom_${jobname}_${samplenum}.out --err=${resultsfolder}logfiles/fast_multinom_${jobname}_${samplenum}.err -A PCAWG
done 

# close file
exec 3<&- 

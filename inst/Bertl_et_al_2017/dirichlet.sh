#!/bin/sh
# to be run on GenomeDK using bash


### settings: ###

# jobname
jobname="pcawg"

# data folders and files for the joint dataset:
datafolder=/project/PCAWG/faststorage/count_final/imputed/
datafolderall=/project/PCAWG/faststorage/count_final/imputed_all/
dirichletfolderall=/project/PCAWG/faststorage/count_final/imputed_dirichlet_all/ # folder for the datasets with Dirichlet prior
rmfolderall=/project/PCAWG/faststorage/count_final/imputed_dirichlet_all_rm/ # folder for the datasets with Dirichlet prior

# delete cancer types with few mutations?
rm_cancer=F

# package locations
libfolder=/project/PCAWG/nullmodel/R_packages/ # where to find R-packages

# number of CV partitions:
cv=5
# number of MIs:
mi=1


### preparations ###

# make resultsfolders
mkdir -p $dirichletfolderall
mkdir -p ${datafolder}../dirichlet_log

# load R
source /com/extra/R/3.1.0/load.sh

# for the joint dataset

mkdir -p ${dirichletfolderall}1
mkdir -p ${rmfolderall}1

for (( mi_index = 1; mi_index <= $mi; mi_index++ ))
do
	qx -v -p -t 12:00:00 -m 64g --queue=express -i=*.R "Rscript run_wrapper_prior_pseudo_counts.R 1 $mi_index $datafolderall $dirichletfolderall $rmfolderall $rm_cancer $libfolder > ${datafolder}../dirichlet_log/dirichlet_${jobname}_all_mi${mi_index}.Rout " | sbatch --job-name=dirichlet_${jobname}_all_mi${mi_index} --out=${datafolder}../dirichlet_log/dirichlet_${jobname}_all_mi${mi_index}.out --mail-user=guo@cs.au.dk --mail-type=END -A PCAWG
done

# This is how the example dataset cancermutations was created (data(cancermutations)).

setwd("/home/johanna/Documents/ThemenAarhus/CancerDriversMutations/R_packages_github/multinomutils")
library(tools)
library(devtools)
library(multinomutils)

# use system.file to find the raw dataset that was installed along with the package:
location = system.file("extdata", "set0", package = "multinomutils")
count.raw = read.table(file=location, header = T, as.is=T)

count.complete = count_table_prep_multinom(count.raw, m=5, impute=F)
cancermutations = as.data.frame(count.complete$imp[[1]])[,-6] # remove the superfluous column i.cancer_type

# make factors: 
cancermutations$sample_id = as.factor(cancermutations$sample_id)
cancermutations$cancer_type = as.factor(cancermutations$cancer_type)
cancermutations$neighbors = as.factor(cancermutations$neighbors)

# reshuffle all lines (to avoid that subsets of the data only contain very few factor levels)
reshuffle = sample.int(nrow(cancermutations))
cancermutations = cancermutations[reshuffle, ]

### Find out what the best way of compression is:
# save(cancermutations, file="/home/johanna/Documents/ThemenAarhus/CancerDriversMutations/R_packages_github/multinomutils/data/cancermutations.rda")
# tools::resaveRdaFiles("/home/johanna/Documents/ThemenAarhus/CancerDriversMutations/R_packages_github/multinomutils/data/cancermutations.rda", compress="auto")
# tools::checkRdaFiles("/home/johanna/Documents/ThemenAarhus/CancerDriversMutations/R_packages_github/multinomutils/data/cancermutations.rda")
### xz is the best.

devtools::use_data(cancermutations, compress="xz", overwrite = T)



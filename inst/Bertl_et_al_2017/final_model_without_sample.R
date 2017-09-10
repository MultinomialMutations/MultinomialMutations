samplevec = c("strong*sample_id")

cancervec = c("strong*cancer_type", "strong*neighbors*cancer_type", "strong*expression*cancer_type", "strong*expression_dummy*cancer_type", "strong*phyloP*cancer_type", "strong*replication_timing*cancer_type")

simplevec = c("strong*neighbors", "strong*phyloP", "strong*replication_timing", "strong*expression", "strong*expression_dummy")

models = matrix(F, ncol=length(samplevec) + length(cancervec) + length(simplevec), nrow=4)
colnames(models) = c(samplevec, cancervec, simplevec)

# with neighbors
models[1,] = c(F, T, T, T, T, T, T, rep(F, 5)) 

# without neighbors
models[2,] = c(T, F, T, T, T, T, T, rep(F, 5)) 

# without neighbors and sample_id, but cancer type
models[3,] = c(F, T, F, T, T, T, T, rep(F, 5)) 

# model with phyloP, neighbors, replication timing and expression for the figure (without interactions with cancer types or samples)
models[4,] = c(rep(F, 7), rep(T, 5)) 

subsets=NULL

intercept=rep(0, 4)

save(models, subsets, intercept, file="final_model_without_sample.RData")

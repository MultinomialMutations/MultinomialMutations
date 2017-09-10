

### Example with DNase1_peak dummy variable ###


# Model with DNase1_peak dummy
# (DNase1 is not in the included dataset chrom21, so I just create a variable with that name here)
names(chrom21)[3] = "DNase1"
chrom21$DNase1_peak = as.factor(as.numeric(chrom21$DNase1>0))

fit4 = fast_multinom(cbind(NO, I, VA, VG) ~ 0 + DNase1_peak*DNase1, data = chrom21, refLevel=1, loglik=T, predictions=F, VC=T)
# doesn't work, because of a bug in glm4 that overspecifies the model (or rather, the model is overspecified, but glm would remove a superfluous term)

# How does it work with glm?
glm(cbind(NO, I) ~ 0 + DNase1_peak*DNase1, data=chrom21, family=binomial)
glm(cbind(NO, I) ~ 0 + DNase1_peak + DNase1_peak:DNase1, data=chrom21, family=binomial)
# same results

# work around the glm4 bug by computing the variables explicitly:
# chrom21$DNase1 is already zero for DNase1_peak==0, so it is the same as DNase1*DNase1_peak (if * means multiplication in the mathematical sense)
fit3 = fast_multinom(cbind(NO, I, VA, VG) ~ 0 + DNase1_peak + DNase1, data = chrom21, refLevel=1, loglik=T, predictions=F, VC=T)

glm(cbind(NO, I) ~ 0 + DNase1_peak + DNase1, data=chrom21, family=binomial)


### Example with missing values ###

# The example dataset doesn't contain missing values, so they are randomly added here:
chrom21.missing = chrom21
chrom21.missing$replication_timing[sample.int(n=nrow(chrom21), size=1000)]=NA

fit.NA = fast_multinom(cbind(NO, I, VA, VG) ~ replication_timing + cancer_type, data = chrom21.missing, refLevel=1, loglik=T, predictions=T, VC=T)
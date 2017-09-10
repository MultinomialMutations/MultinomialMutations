#' @rdname wrapper_loss

wrapper_loss_binom = function(cv, cv_index, mi_index, model_index, modelfile, datafolder, resultsfolder, per_obs, nested_samples=T){

  # load data
  # load the training data (complement of CV slice i)
  cvslices=(1:cv)[-cv_index]
  training = fread(input=paste0(datafolder, cvslices[1], "/imp", mi_index, ".txt"), header = T)
  for(i in 2:(cv-1)){
    training = rbind(training, read.table(file=paste0(datafolder, cvslices[i], "/imp", mi_index, ".txt"), header = T))
  }
  # aggregate rows with the same combination of categories (using data.table)
  training[, zero:= NULL]
  byvec = colnames(training)
  byvec = byvec[!(byvec %in% c("NO", "VA", "VG", "I", "YES"))]
  setkeyv(training, byvec) # using setkeyv for a character vector of column names
  training <- training[, .("NO" = sum(NO), "VA" = sum(VA), "VG" = sum(VG), "I" = sum(I), "YES" = sum(YES)), by=key(training)]

  # make dummies factors
  # make strong a factor
  training$strong = as.factor(training$strong)
  # make the other binary variables factors if they exist
  if(exists("Cgi", where = training)) training$Cgi = as.factor(training$Cgi)
  if(exists("simple_repeat", where = training)) training$simple_repeat = as.factor(training$simple_repeat)
  if(exists("DNase1_peak", where = training)) training$DNase1_peak = as.factor(training$DNase1_peak)
  if(exists("expression_dummy", where = training)) training$expression_dummy = as.factor(training$expression_dummy)

  # set nested treatment contrasts for sample_id (nested within cancer_type)
  if(nested_samples){
    contrasts(training$cancer_type) = contr.treatment(nlevels(training$cancer_type))
    nesting = nested_treatment_contrasts(outer.factor=training$cancer_type, inner.factor=training$sample_id)
    how.many = nlevels(training$sample_id) - nlevels(training$cancer_type)
    contrasts(training$sample_id, how.many=how.many) = nesting
  }


  # load model
  load(modelfile)

  if(sum(models[model_index,])>0){
    multinomformula = formula(paste("cbind(NO, YES) ~ ", paste(intercept), "+", paste(colnames(models)[models[model_index,]], collapse="+")))
  } else {
    multinomformula = formula("cbind(NO, YES) ~ 1") # empty model
  }

  # fit model:
  fit.fast_multinom = fast_multinom(multinomformula, data = training, family = binomial(link="logit"), VC=F, subsetmatrix=subsets[[model_index]], parallel=F)
  # save model fit
  save(fit.fast_multinom, file=paste0(resultsfolder, "fit_mi", mi_index, "cv", cv_index, "model", model_index, ".Rdata"))

  #### estimate loss ####

  # load data
  validation = read.table(file=paste0(datafolder, cv_index, "/imp", mi_index, ".txt"), header = T)

  # make dummies factors
  # make strong a factor
  validation$strong = as.factor(validation$strong)
  # make the other binary variables factors if they exist
  if(exists("Cgi", where = validation)) validation$Cgi = as.factor(validation$Cgi)
  if(exists("simple_repeat", where = validation)) validation$simple_repeat = as.factor(validation$simple_repeat)
  if(exists("DNase1_peak", where = validation)) validation$DNase1_peak = as.factor(validation$DNase1_peak)
  if(exists("expression_dummy", where = validation)) validation$expression_dummy = as.factor(validation$expression_dummy)

  # estimate loss
  est.loss = deviance_loss(fit.fast_multinom, validation, per_obs)

  # save loss
  save(est.loss, file=paste0(resultsfolder, "loss_mi", mi_index, "cv", cv_index, "model", model_index, ".Rdata"))

}

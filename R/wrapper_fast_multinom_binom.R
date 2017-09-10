#' @rdname wrapper_fast_multinom

wrapper_fast_multinom_binom = function(model.index, mi, datafolder, modelfile, resultsfolder, startfile, contrasts="all_sum", nested_samples=T, VC=T, loglik=T){

  # load data
  countmat = read.table(paste0(datafolder, "imp", mi, ".txt"), header=T)

  # make strong a factor
  countmat$strong = as.factor(countmat$strong)
  # make the other binary variables factors if they exist
  if(exists("Cgi", where = countmat)) countmat$Cgi = as.factor(countmat$Cgi)
  if(exists("simple_repeat", where = countmat)) countmat$simple_repeat = as.factor(countmat$simple_repeat)
  if(exists("DNase1_peak", where = countmat)) countmat$DNase1_peak = as.factor(countmat$DNase1_peak)
  if(exists("expression_dummy", where = countmat)) countmat$expression_dummy = as.factor(countmat$expression_dummy)

  # nested treatment contrasts for the samples
  if(nested_samples){
    nesting = nested_treatment_contrasts(outer.factor=countmat$cancer_type, inner.factor=countmat$sample_id)
    how.many = nlevels(countmat$sample_id) - nlevels(countmat$cancer_type)
    contrasts(countmat$sample_id, how.many=how.many) = nesting
    # corrected bug here
  }

  if(contrasts=="all_sum"){
    # set sum contrasts for cancer_type, sample_id (nested within cancer_type), genomicSeg, neighbors
    contrasts(countmat$cancer_type) = contr.sum(nlevels(countmat$cancer_type))
    contrasts(countmat$neighbors) = contr.sum(nlevels(countmat$neighbors))

    if(exists("genomicSeg", where = countmat)) contrasts(countmat$genomicSeg) = contr.sum(nlevels(countmat$genomicSeg))

    if(nested_samples){
      nesting = nested_sum_contrasts(outer.factor=countmat$cancer_type, inner.factor=countmat$sample_id)
      how.many = nlevels(countmat$sample_id) - nlevels(countmat$cancer_type)
      contrasts(countmat$sample_id, how.many=how.many) = nesting
    }
  }

  # load models
  load(modelfile)
  # load starting values
  if(is.null(startfile)){
    startvec=NULL
  } else {
    load(paste0(resultsfolder, startfile))
    startvec = startlist[[model.index]]
  }
  # formula
  if(sum(models[model.index,])>0){
    multinomformula = formula(paste("cbind(NO, YES) ~ ", paste(intercept[model.index]), "+", paste(colnames(models)[models[model.index,]], collapse="+")))
  } else {
    multinomformula = formula(paste("cbind(NO, YES) ~ 1"))
  }


  #### estimation ####
  final.fast_multinom = fast_multinom(multinomformula, data = countmat, family = binomial(link="logit"), VC=VC, subsetmatrix=subsets[[model.index]], loglik=loglik, parallel=F)


  #### save results ####

  warn = warnings() # in case glm4 produced warnings

  save(final.fast_multinom, warn, file=paste0(resultsfolder, "final_fastmultinom_model", model.index, "mi", mi, ".Rdata"))
}

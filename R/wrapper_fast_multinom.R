#' Wrapper for the function fast_multinom
#'
#' The function \code{wrapper_fast_multinom} estimates the regression coefficients of a multinomial logistic model with \code{\link{fast_multinom}}. This wrapper was used in our analysis in Bertl et al. (2007) (see References). The function \code{wrapper_fast_multinom_binom} uses a binomial model instead.
#'
#' A dataset similar to the example data set \code{cancermutations} is read. The dummy variables \code{strong}, \code{Cgi}, \code{simple_repeat}, \code{DNase1_peak}, \code{expression_dummy} are changed to factors, optionally, sum contrasts are set to cancer_type, genomicSeg and sample_id, and nested contrasts are set to create a nesting of sample_id in cancer_type.
#'
#' Then, the multinomial regression model indexed by \code{model.index} is obtained from a model file and estimated with the function \code{\link{fast_multinom}}.
#'
#' The scripts that were used to run this function and that show all settings used in Bertl et al. (2007) are available in this package in the folder \code{inst/Bertl_et_al_2017}. The pre-processed data can be downloaded from figshare.
#'
#' @param model.index integer. Number of the model in the model matrix.
#' @param mi integer. Number of multiple imputation replicate of the data set.
#' @param datafolder character. Folder, where the data is saved. See data(cancermutations) for an example. The filename is \code{paste0("imp", mi, ".txt")}, where \code{mi} is the multiple imputation replicate.
#' @param modelfile character. File that contains the models in the form of a matrix (see examples).
#' @param resultsfolder character. Where to save the results?
#' @param startfile character. File that contains starting values for the parameter estimation.
#' @param contrasts character. What contrasts should be set for the variables cancer_type, genomicSeg and sample_id. Default: 'all_sum'. All other settings will use the default contrasts, i. e. treatment contrasts.
#' @param nested_samples Logical. Should the contrasts be defined such that the samples are nested within the cancer types? Default: T
#' @param VC Logical. Should the VC matrix be computed? Default: T
#' @param loglik Logical. Should the log-likelihood be computed? Default: T
#'
#' @return There is no output. The estimates are saved along with (potential) warning messages.
#'
#' @seealso \code{\link[multinomutils::fast_multinom]{fast_multinom}}
#'
#' @examples
#'
#' @rdname wrapper_fast_multinom
#'
#' @references
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#'
#' @author Johanna Bertl

wrapper_fast_multinom = function(model.index, mi, datafolder, modelfile, resultsfolder, startfile, contrasts="all_sum", nested_samples=T, VC=T, loglik=T){

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
    if(nlevels(countmat$cancer_type)>1) contrasts(countmat$cancer_type) = contr.sum(nlevels(countmat$cancer_type))
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
    multinomformula = formula(paste("cbind(NO, I, VA, VG) ~ ", paste(intercept[model.index]), "+", paste(colnames(models)[models[model.index,]], collapse="+")))
  } else {
    multinomformula = formula(paste("cbind(NO, I, VA, VG) ~ 1"))
  }


  #### estimation ####
  final.fast_multinom = fast_multinom(multinomformula, data = countmat, family = binomial(link="logit"), VC=VC, subsetmatrix=subsets[[model.index]], loglik=loglik, parallel=F)


  #### save results ####

  warn = warnings() # in case glm4 produced warnings

  save(final.fast_multinom, warn, file=paste0(resultsfolder, "final_fastmultinom_model", model.index, "mi", mi, ".Rdata"))
}

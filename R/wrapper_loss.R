#' Wrapper to estimate the deviance loss by cross-validation
#'
#' The function \code{wrapper_loss} estimates the deviance loss in a multinomial regression model by leave-one-out cross validation using \code{\link{fast_multinom}} and \code{\link{deviance_loss}}. This wrapper was used in our analysis in Bertl et al. (2007) (see References). The function \code{wrapper_loss_binom} uses a binomial model instead.
#'
#' This function estimates a multinomial regression model on the joint set of all cross validation pieces that dataset has been divided into except \code{cv_index}. Then, the deviance loss is estimated on the dataset \code{cv_index}. In a further step, the function wrapper_average_loss should be used for averaging over the loss estimates.
#'
#' The data is prepared and the regression is estimated as in \code{\link{wrapper_fast_multinom}}. As the contrasts are irrelevant for prediction, they cannot be set here. By default, nested contrasts are used for the sample to avoid overspecifying the model (because this is not handled correctly by the function \code{\link[MatrixModels]{glm4}}, see \code{\link{fast_multinom}} for details. The option nested_samples allows to remove the nesting, if the cancer_type is not part of the model.
#'
#' The scripts that were used to run this function and that show all settings used in Bertl et al. (2007) are available in this package in the folder \code{inst/Bertl_et_al_2017}. The pre-processed data can be downloaded from figshare.
#'
#' @param cv integer. Number of pieces the dataset has been divided into for cross validation.
#' @param cv_index integer. Which cross-validation slice is currently used?
#' @param mi_index integer. Number of multiple imputation replicate.
#' @param model_index integer. Number of the model in the model matrix.
#' @param modelfile character. File that contains the models in the form of a matrix (see examples).
#' @param datafolder character. Folder that contains the dataset at the location \code{paste0(datafolder, cvslices[1], "/imp", mi_index, ".txt")}.
#' @param resultsfolder character. Where to save the estimated loss and the estimated regression model. Note that the VC matrix is not saved.
#' @param per_obs logical. If per_obs==T, the loss is normalized by the total number of observations (sum of all counts), so it is the mean loss.
#' @param nested_samples logical. Default=T. Are the samples nested in the cancer types?
#'
#' @return There is no output. The regression coefficients and the loss estimate are saved.
#'
#' @seealso \code{\link{fast_multinom}}, \code{\link{deviance_loss}}, \code{\link{wrapper_fast_multinom}}
#'
#' @examples
#'
#' @rdname wrapper_loss
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#'
#' @author Johanna Bertl
#'
wrapper_loss = function(cv, cv_index, mi_index, model_index, modelfile, datafolder, resultsfolder, per_obs, nested_samples=T){

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
    multinomformula = formula(paste("cbind(NO, I, VA, VG) ~ ", paste(intercept), "+", paste(colnames(models)[models[model_index,]], collapse="+")))
  } else {
    multinomformula = formula("cbind(NO, I, VA, VG) ~ 1") # empty model
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

#' ### TESTING THE PREDICTION FUNCTION ###
#'
#' ########################################################################
#' # Basic example ########################################################
#' ########################################################################
#' data(chrom21)
#'
#' # the APOBEC signature is only relevant for transitions and transversions to a G:C basepair -- construct the corresponding subset of parameters for the 3 binomial models:
#' subs = matrix(T, ncol=3, nrow=4)
#' subs[2,2] = F
#'
#' # fit the multinomial model
#' fit = fast_multinom(cbind(NO, I, VA, VG) ~ strong + apobec + cancer_type, data = chrom21, refLevel=1, VC=T, subsetmatrix=subs, predictions=T)
#'
#' # predictions on the data that was used for fitting (only available, because predictions=T in the function fast_multinom):
#' head(predict.fast_multinom(fit))
#'
#' # predict on new data (with fewer factor levels):
#' set.seed(123)
#' new = droplevels(chrom21[sample.int(nrow(chrom21), 5),])
#' pred = predict.fast_multinom(fit, new)
#' 
#'
#' #########################################################################
#' # Test specifying different contrasts ###################################
#' #########################################################################
#'
#' data(chrom21)
#'
#' # the APOBEC signature is only relevant for transitions and transversions to a G:C basepair -- construct the corresponding subset of parameters for the 3 binomial models:
#' subs = matrix(T, ncol=3, nrow=4)
#' subs[3,2] = F
#'
#' # fit the multinomial model with default (treatment) contrasts
#' fit = fast_multinom(cbind(NO, I, VA, VG) ~ strong + apobec + cancer_type, data = chrom21, refLevel=1, VC=T, subsetmatrix=subs, predictions=F)
#'
#' set.seed(123)
#' subsample = sample.int(nrow(chrom21), 100)
#' new = chrom21[subsample,]
#' pred = predict.fast_multinom(fit, new)
#'
#' mM = model.Matrix(fit$formulae$I, data = new, sparse=T, xlev = fit$xlevels$I, contrasts.arg = fit$contrasts$I)
#'
#' # sum contrasts
#' chrom21.sc = chrom21
#' contrasts(chrom21.sc$cancer_type) = contr.sum(nlevels(chrom21.sc$cancer_type))
#'
#' # fit the same multinomial model
#' fit.sc = fast_multinom(cbind(NO, I, VA, VG) ~ strong + apobec + cancer_type, data = chrom21.sc, refLevel=1, VC=T, subsetmatrix=subs, predictions=F)
#'
#'
#' # predict on the same data and compare:
#' new.sc = chrom21.sc[subsample,]
#' pred.sc = predict.fast_multinom(fit.sc, new.sc)
#' ### there shouldn't be a warning!!! ###
#' identical(pred.sc, pred)
#' max(abs(pred.sc - pred))
#' # not identical, but practically the same (differences of the order of 10^(-20))
#'
#'
#'
#' ########################
#'
#' #### The examples below might not work any more, because I changed the above examples, and I think they depend on them. #################################################################################
#'
#' # What if a factor in the new data has fewer levels than in the original dataset?
#' new.drop = droplevels(new)
#' pred.drop = predict.fast_multinom(fit.sc, new.drop)
#' # Okay.
#'
#' # Test if an interaction term btw. a factor and a numerical variable (which is not a dummy) works correctly with different contrasts:
#' # Treatment contrasts:
#' contrasts(chrom21$cancer_type) = contr.treatment(nlevels(chrom21$cancer_type))
#' fit.treat = fast_multinom(cbind(NO, I, VA, VG) ~ replication_timing*cancer_type, data = chrom21, refLevel=1, VC=T, subsetmatrix=NULL)
#'
#' pred.treat = predict.fast_multinom(fit.treat, new)
#'
#' # Treatment contrasts with a different base:
#' contrasts(chrom21$cancer_type) = contr.treatment(nlevels(chrom21$cancer_type), base = 2)
#' fit.treat2 = fast_multinom(cbind(NO, I, VA, VG) ~ replication_timing*cancer_type, data = chrom21, refLevel=1, VC=T, subsetmatrix=NULL)
#'
#' pred.treat2 = predict.fast_multinom(fit.treat2, new)
#'
#' # Comparison
#' identical(pred.treat, pred.treat2)
#' max(abs(pred.treat - pred.treat2))
#' # Practically the same. Largest difference: 10^{-20}
#'
#' # Sum contrasts:
#' contrasts(chrom21$cancer_type) = contr.sum(nlevels(chrom21$cancer_type))
#' fit.sum = fast_multinom(cbind(NO, I, VA, VG) ~ replication_timing*cancer_type, data = chrom21, refLevel=1, VC=T, subsetmatrix=NULL)
#' ### warning (the warning stems from model.Matrix used within glm4)
#'
#' pred.sum = predict.fast_multinom(fit.sum, new)
#' ### same warning
#' max(abs(pred.treat - pred.sum))
#' ### not the same as pred.treat. Largest difference: 10^{-5}
#'
#' ### Comparisons between the results of model.Matrix and model.matrix for this model show differences, but this is not the case for the 2 versions of treatment contrasts.
#'
#' ### Conclusion: sum.contrasts are not handled correctly when involved in interactions (only tested interactions with a continuous variable). Treatment contrasts are handled correctly, no matter what the baseline is.
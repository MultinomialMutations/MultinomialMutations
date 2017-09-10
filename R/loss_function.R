#' Deviance loss function for multinomial regression
#'
#' Computing the deviance loss function for a multinomial regression estimate
#'
#' The deviance loss function for categorical data is estimated as described for example in Hastie et al, p. 221.
#'
#' @param fit an object of class fast_multinom, fit on the training data.
#' @param validation a validation dataset.
#' @param per_obs logical. If per_obs==T, the loss is normalized by the total number of observations (sum of all counts), so it is the mean loss.
#' 
#' @seealso \code{\link{fast_multinom}}
#'
#'
#' @examples
#' data(cancermutations)
#'
#' # the APOBEC signature is only relevant for transitions and transversions to a G:C basepair -- construct the corresponding subset of parameters for the 3 binomial models:
#' subs = matrix(T, ncol=3, nrow=4)
#' subs[3,2] = F
#'
#' # fit the multinomial model on the first half of the data
#' split = nrow(cancermutations)/2
#' cancermutations_est = cancermutations[1:split,]
#' cancermutations_loss = cancermutations[(split+1):(split*2),]
#' fit = fast_multinom(cbind(NO, I, VA, VG) ~ strong + apobec + cancer_type, data = cancermutations_est, refLevel=1, VC=T, subsetmatrix=subs)
#'
#' # estimation of the loss function on second half of the data:
#' deviance_loss(fit, cancermutations_loss)
#'
#' @author Johanna Bertl
#'
#' @references
#'
#' Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#' Hastie, T.; Tibshirani, R. & Friedman, J. The Elements of Statistical Learning Springer New York Inc., 2001


deviance_loss = function(fit, validation, per_obs = F){

  # predicted probabilities for all categories
  test.pred = predict(fit, newdata = validation)

  # only the counts of the validation set
  counts = validation[colnames(test.pred)]

  # loss estimate
  lossest = -2 * sum(counts*log(test.pred))

  # if per_obs == T: compute the loss estimate per observation
  if(per_obs) lossest = lossest/sum(counts)

  return(lossest)
}

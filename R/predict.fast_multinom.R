#' Prediction method for fast_multinom fits
#'
#' Obtains predictions from a multinomial model fitted using fast_multinom.
#'
#' Only output on the level of probabilities are generated. This is equivalent to type="response" in a logistic regression model estimated with glm with family=binomial(link="logit").
#'
#' Note that the factor levels and contrasts defined in the fast_multinom object are used. This allows to run prediction on a dataset with fewer levels and different contrasts (in this case, the contrasts are overwritten and there is a warning. There can be the same warning, when the contrasts are not changed -- I don't know why, it is issued by \code{\link[MatrixModels]{model.Matrix}}.)
#'
#' @param object Object of class fast_multinom.
#' @param newdata Optionally, a data frame for which to predict.
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#' @examples
#' 
#' data(cancermutations)
#'
#' # the APOBEC signature is only relevant for transitions and transversions to a G:C basepair -- construct the corresponding subset of parameters for the 3 binomial models:
#' subs = matrix(T, ncol=3, nrow=4)
#' subs[3,2] = F
#'
#' # fit the multinomial model
#' fit = fast_multinom(cbind(NO, I, VA, VG) ~ strong + apobec + cancer_type, data = cancermutations, refLevel=1, VC=T, subsetmatrix=subs, predictions=T)
#'
#' # predictions on the data that was used for fitting (only available, because predictions=T in the function fast_multinom):
#' head(predict.fast_multinom(fit))
#'
#' # predict on new data (with fewer factor levels):
#' set.seed(123)
#' new = droplevels(cancermutations[sample.int(nrow(cancermutations), 5),])
#' pred = predict.fast_multinom(fit, new)
#'
#'
#' @author Johanna Bertl


predict.fast_multinom <- function(object, newdata=NULL){

  # Predictions that were part of the fit
  if (is.null(newdata)) {
    y <- object$predictions
    if(is.null(y)) warning("The predictions on the original data couldn't be retreived, because they weren't saved along with the fit (predictions=F in the function call fast_multinom).")
  } else {
  # Predictions on new data

    if(is.data.table(newdata)){
      warning("The class of the new data is changed from data.table to data.frame.")
      data = as.data.frame(data)
    }

    M = length(object$formulae)


    # preparing data structures
    n = nrow(newdata)
    lin.pred = matrix(NA, ncol=M, nrow=n)

    # newdata also has to contain output columns in order to obtain the model matrix with the function model.Matrix (below). I just add them here. (Should find a better solution! This matrix might be large.)
    left.side = as.character(attr(terms(object$formulae[[1]]), "variables")[2])
    refname.regexp = regexpr("(?<=, )\\w+(?=\\)$)", left.side, perl=T)
    # regular expression: perl=T is used, because lookaround is not part of the R regexpr. Finds any "word characters" (letters, digits, underscores) between ", " and ")" at the end of the word.
    refname = regmatches(left.side, refname.regexp)
    if(length(regmatches)==0) {stop("Can't find out the name of the reference category")}
    outputnames = c(refname, names(object$coefficients))
    for(name in outputnames){
      if(!name %in% colnames(newdata)){
        currnames = names(newdata)
        newdata = cbind(newdata, numeric(n))
        names(newdata) = c(currnames, name)
      }
    }

    for(m in 1:M){
      mm = model.Matrix(object$formulae[[m]], data = newdata, xlev = object$xlevels[[m]], contrasts.arg = object$contrasts[[m]], sparse=T)
      lin.pred[,m] = as.numeric(mm %*% object$coefficients[[m]]) # (as.numeric is not a very elegant solution!)
    }

    # probabilities
    logits = exp(lin.pred)
    phi0 = 1/(1 + rowSums(logits))

    # predictions for the joint model (phi matrix)
    y = cbind(phi0, logits*phi0)
    colnames(y) = c(all.vars(object$formulae[[1]])[2], names(object$coefficients))

  }

  return(y)

}

#' Fast multinomial logistic regression on large and sparse datasets.
#'
#' \code{fast_multinom} estimates an M-category multinomial logistic regression
#' model by M-1 separate binomial logistic regression models and estimates the
#' VC matrix according to Begg and Gray (1984).
#'
#' The separate regression models are estimated with the function \code{\link[MatrixModels]{glm4}} from package \code{MatrixModels}, which is very fast for sparse design
#' matrices (usually models including many factors and interactions with factors). For not sparse design matrices, a different function like
#' \code{\link[speedglm]{speedglm}} or \code{\link[biglm]{bigglm}} might be preferable, but this is not
#' implemented here.
#'
#' The computation of the VC matrix follows Begg & Gray (1984). As very large and sparse matrices are involved, the \code{Matrix} package is used here. The computation of the VC matrix is very fast, but computations can be sped up by setting VC=F. This can be useful if the model is only used for prediction or estimation of a loss function.
#'
#' The function \code{\link[VGAM]{vglm}} from package \code{VGAM} provides a more flexible interface to multinomial logistic regression, has much
#' more functionalities and a more extensive output, and also allows the
#' specification of joint parametes and many different (custom) link functions.
#' It is certainly preferable for many problems, but when it reaches its
#' boundary in terms of runtime and memory, this function can be an alternative.
#'
#' If the data is a data.table, it is internally changed to a data.frame. The code of this function is developed for data.frames and data.table is not compatible with this (see e. g. FAQ 1.1 in the data.table package for examples).
#'
#' @section Caution!:
#' Note that there is a bug in \code{\link[MatrixModels]{glm4}} that produces erroneous results when contrasts other than the default treatment contrasts are used for a factor that is involved in an interaction with a numeric variable (see \url{http://r-forge.r-project.org/tracker/index.php?func=detail&aid=6192&group_id=61&atid=294} for my bug report). The function doesn't find these errors, so only use treatment contrasts for that case. If the numeric variable is a dummy variable, it is enough to change it to a factor.
#'
#' Note that the function \code{\link[MatrixModels]{model.Matrix}} that is internally used by \code{\link[MatrixModels]{glm4}} has a bug, s. th. models with interaction terms that are specified with ":" (partial interaction) instead of "*" (full interaction) are often not handled correctly and result in model matrices that don't have full rank. The function doesn't catch these errors, but usually the underlying Cholesky decomposition will stop with an error message. If unsure about the specified model, create the model matrix first with \code{model.Matrix(..., sparse=T, ...)} and compare its rank with its dimension using the function \code{\link[Matrix]{rankMatrix}}.
#'
#' Note that it is currently impossible to run the function on data with missing values. glm4 allows standard handling of missing values (with na.action), but I have not implemented handling missing data (and as a consequence, changing sizes of predicted vs observed data, and potentially different numbers of missing values in the binary submodels). I set na.action=na.fail in glm4.
#'
#'
#' @seealso \code{\link[VGAM]{vglm}}, \code{\link[MatrixModels]{glm4}}, \code{\link{predict.fast_multinom}}
#'
#'
#' @references Begg, C. B. & Gray, R. Calculation of polychotomous logistic regression parameters using individualized regressions. Biometrika, 1984, 71, 11-18
#'
#' Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#' @param multinomformula A formula of the form cbind(category1, category2, ...,
#'   categoryM) ~ explanatory1 + explanatory2 + ... + explanatory3).
#' @param data Dataframe.
#' @param refLevel Reference level.
#' @param family I don't think I implemented anything except binomial with logit link.
#' @param subsetmatrix Logical matrix with M-1 columns and the same number of
#'   rows as there are explanatory variables (including the intercept, even if it is 0 or -1). F
#'   indicates that an explanatory variable is not included in the submodel.
#' @param VC logical Should the variance-covariance-matrix of the estimates be
#'   calculated?
#' @param predictions logical Should the predictions (fitted values) for the data be calculated?
#' @param loglik logical Should the log likelihood be calculated? Note that the predictions have to be computed to obtain the log likelihood. Setting loglik=T and predictions=F saves space when saving the results, but it doesn't save time.
#' @param parallel logical. Should the estimation of the M-1 submodels be
#'   distributed between M-1 cores? Not implemented yet.
#'
#' @return An object of class fast_multinom (an S3 class, a list, actually) that consists of the following elements:
#' \describe{
#'  \item{formulae}{A named list that contains the formulae for the separate logistic regression models.}
#'  \item{coefficients}{A named list that contains the estimated regression coefficients for each of the M-1 submodels. }
#'  \item{VC}{Variance covariance matrix of the coefficients, if VC=T, otherwise NULL. Some type of \code{Matrix} object.}
#'  \item{gamma}{NULL (the mean fraction of missing information, only computed by the function \code{\link{pooling}}).}
#'  \item{predictions}{A matrix of predicted probabilities obtained from the current dataset, if predictions=T, otherwise NULL. The columns correspond to the M categories and each line corresponds to one observation.}
#'  \item{log.lik}{Log-likelihood of the multinomial model, if loglik=T, otherwise NULL.}
#'  \item{contrasts}{A named list that contains a named list of the contrasts of each factor in each of the M-1 submodels. This is needed for prediction with \code{\link{predict.fast_multinom}}. The contrasts can be specified either as a contrast matrix or with a character string with the name of the function to make them (usually, when the default contrasts are used).}
#'  \item{xlevels}{A named list that contains a named list of the factor levels for each factor in each of the M-1 submodels. This is needed for prediction with \code{\link{predict.fast_multinom}}, when the new data set contains less levels.}
#'  \item{n}{Number of observations (sum over the weights).}
#' }
#'
#'
#' @examples
#' # Load the internal example dataset:
#' data(cancermutations)
#'
#' # the APOBEC signature is only relevant for transitions and transversions to a G:C basepair -- construct the corresponding subset of parameters for the 3 binomial models:
#' subs = matrix(T, ncol=3, nrow=4)
#' subs[3,2] = F
#'
#' # fit the multinomial model
#' fit = fast_multinom(cbind(NO, I, VA, VG) ~ strong + apobec + cancer_type, data = cancermutations, refLevel=1, loglik=T, predictions=F, VC=T, subsetmatrix=subs)
#'
#' # The same without intercept
#' fit2 = fast_multinom(cbind(NO, I, VA, VG) ~ 0 + strong + apobec + cancer_type, data = cancermutations, refLevel=1, loglik=T, predictions=F, VC=T, subsetmatrix=subs)
#'
#' # Simpler model without subsetting and continuous variable
#' fit3 = fast_multinom(cbind(NO, I, VA, VG) ~ replication_timing, data = cancermutations, refLevel=1, loglik=T, predictions=F, VC=T, subsetmatrix=NULL)
#'
#'
#' @author Johanna Bertl
#'

fast_multinom <-
function(multinomformula, data, refLevel=1, family = binomial(link="logit"), subsetmatrix = NULL, VC=F, predictions = F, loglik=F, parallel = F){

  #### preparations ####

  # changing the class of the data to data frame, if it is a data table

  if(is.data.table(data)){
    #warning("The class of the data is changed from data.table to data.frame.")
    data = as.data.frame(data)
  }

  # preparing the formulae

  # left sides
  output = gsub(" ", "", as.character(attr(terms(multinomformula), "variables")[2]))
  categories = unlist(strsplit(substr(output, 7, nchar(output)-1), ",")) ### only works with cbind in the formula
  # number of models
  M = length(categories) - 1

  # right sides
  # variables = c(attr(terms(multinomformula), "intercept"), attr(terms(multinomformula), "term.labels"))
  # Theoretically, this is a good solution, because it uses the decomposed formula where full interactions specified with "*" have been written out in detail using the main effects and interactions with ":". Problem: glm4 uses model.Matrix, which handles interactions with ":" erroneously. So I have to replace this by string manipulation:
  variables = unlist(strsplit(as.character(multinomformula)[3], "+", fixed=T))
  # add 1 for the intercept if it isn't there, but the formula contains an intercept
  if(attr(terms(multinomformula), "intercept") == 1 & !(grepl("\\<1\\>", as.character(multinomformula)[3]))) {
    # the regular expression finds a "word" that only consists of 1. This is save, because as.character inserts blanks between the elements of the formula (e. g. strong+1+apobec becomes strong + 1 + apobec)
    variables = c(1, variables)
  }

  # make a full subsetmatrix if none is given
  if(is.null(subsetmatrix)){
    subsetmatrix = matrix(T, ncol=M, nrow=length(variables))
  }
  # check if the number of rows of the submset matrix is correct.
  if(nrow(subsetmatrix) != length(variables)) stop("The subset matrix doesn't have the same number of rows as there are terms in the formula. Maybe you forgot the intercept.")

  # go on with the right sides
  formula.rightside = character(M)
  for(m in 1:M)
  {
    formula.rightside[m] = paste(variables[subsetmatrix[,m]], collapse="+")
  }

  # the whole formulae
  logisticformula = vector("list", M)
  for(m in 1:M)
  {
    logisticformula[[m]] = formula(paste0("cbind(", (categories[-refLevel])[m], ", ", categories[refLevel], ") ~ ", formula.rightside[m]))
  }

  # number of cases and weights
  count.cols = which(names(data) %in% categories)
  weights = rowSums(as.matrix(data[,count.cols]))
  n = sum(weights)


  #### Estimation ####

  coef.list = vector("list", M) #### should be a matrix (with NAs) and contain the names!!!
  design.list = vector("list", M)
  #xlevels = vector("list", M)
  contrast.list = vector("list", M)
  levels.list = vector("list", M)
  pred.mat = matrix(NA, nrow=nrow(data), ncol=M)

  names(coef.list) = categories[-refLevel]

  if(parallel) {
    stop("Parallel processing is not implemented yet. Set parallel=F.")
  } else {

    for(m in 1:M) {
      fit = glm4(formula(logisticformula[[m]]), data = data, family = family, sparse=T, na.action=na.fail)
      coef.list[[m]] = fit@pred@coef
      names(coef.list[[m]]) = fit@pred@X@Dimnames[[2]]
      design.list[[m]] = fit@pred@X
      pred.mat[,m] = fitted(fit)
#       factornames = names(fit@pred@X@contrasts)
#       lfn = length(factornames)
#       if(lfn>0){
#         xlevels[[m]] = vector("list", lfn)
#         for(fnindex in 1:lfn){
#           #if(is.data.table(data)) {
#            # xlevels[[m]][[fnindex]] = levels(as.factor(data[,factornames[fnindex], with=F]))
#           #} else {
#             xlevels[[m]][[fnindex]] = levels(as.factor(data[,factornames[fnindex]]))
#           #}
#         }
#         names(xlevels[[m]]) = factornames
#       }
      contrast.list[[m]] = fit@pred@X@contrasts

      factornames = names(fit@pred@X@contrasts)
      levels.list[[m]] = vector("list", length(factornames))
      if(length(factornames)>0){
        for(i in 1:length(factornames)){
          levels.list[[m]][[i]] = levels(data[,which(factornames[i] == colnames(data))])
        }
        names(levels.list[[m]]) = factornames
      }



      rm(fit)
    }
  }

  #### Computation of the VC matrix, predictions and log likelihood ###

  log.lik=NULL
  predictions.multi=NULL
  sigma=NULL

  if(VC | predictions | loglik){

    # predictions for the separate binomial models:
    logits = pred.mat/(1 - pred.mat)
    phi0 = 1/(1 + rowSums(logits))

    # predictions for the joint model (phi matrix) -- not needed for the computation of VC
    if (loglik | predictions) {
      # predictions for the joint model
      predictions.multi = cbind(phi0, logits * phi0)
      colnames(predictions.multi) = c(categories[refLevel], categories[-refLevel])

      if(loglik) {
        log.lik = sum(data[colnames(predictions.multi)]*log(predictions.multi))
      }

      if (!predictions) {
        predictions.multi = NULL
      }

    }

    if(VC){

      # upscaling, because unlike in Begg & Gray, we have multiple observations per line
      phi0 = phi0 * weights

      # the matrix underline(S)
      S = bdiag(design.list)
      #print(class(S))

      # constructing the matrix G
      # G = ( diag(vec(theta)) x (EINS kronecker diag (phi0)) x diag(vec(theta)) ) * diag(vec(1/theta))
      # where x is matrix multiplication and * element-wise multiplication

      theta.vec = as.numeric(pred.mat)
      theta.mat = Diagonal(x=theta.vec) # sparse diagonal matrix
      phi.mat = matrix(1, 3, 3) %x% Diagonal(x=phi0)
      G = theta.mat %*% phi.mat %*% theta.mat #crossprod(crossprod(theta.mat, phi.mat), theta.mat)
      diagvec = Matrix::diag(G)/theta.vec
      Matrix::diag(G) = diagvec # correcting the diagonal

      # final W matrix:
      W = t(S) %*% G %*% S /n #crossprod(S, (G %*% S)) / n

      # H matrix:
      pd = cumsum(unlist(lapply(design.list, ncol))) # indices of the block matrices along the diagonal
      H =  make.bdiag(W, pd) # make.bdiag function in make.bdiag.R
      Hinv = solve(H)

      # final VC matrix:
      sigma = (Hinv %*% W %*% Hinv) / n

    }

  }


  #### Output ####

  # names for the lists in the output that represent the submodels
  names(logisticformula) = names(coef.list)
  names(contrast.list) = names(coef.list)
  names(levels.list) = names(coef.list)

  output = list(multinomformula = multinomformula, formulae = logisticformula, coefficients = coef.list, VC = sigma, gamma = NULL, predictions = predictions.multi, contrasts = contrast.list, xlevels = levels.list, log.lik = log.lik, n=n)
  class(output) = "fast_multinom"

  output
}

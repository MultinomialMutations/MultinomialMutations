#' Pooling for multinomial regression
#'
#' The function \code{pool_multinom} pools multinomial regression estimates obtained on multiple imputation data.
#'
#' Pooling of the VC matrix is implemented following Schafer, p. 113-114.
#'
#' @param reslist list of fast_multinom objects.
#' @param VC logical. Should the variance covariance matrix be pooled? Note that for this all entries of reslist need to have a not null VC component. When VC=T, the fraction of missing information (gamma) is also estimated.
#'
#' @return object of class fast_multinom (see \code{\link{fast_multinom}} for details).
#'
#' @references J L Schafer. Analysis of Incomplete Multivariate Data. Chapman & Hall, 1996
#'
#' @examples
#'
#' # use a subset of the raw dataset that was installed along with the package:
#' location = system.file("extdata", "set0", package = "multinomutils")
#' count_table = fread(file=location)
#' count_table = count_table[count_table$cancer_type %in% c("KICH", "LGG")]
#' count_table_imp = count_table_prep_multinom(count_table, 2)
#'
#' # fit the same model to each MI replicate
#' fitlist = vector("list", 2)
#' for(i in 1:2){
#'  countmat = count_table_imp[[1]][[i]]
#'
#'  # set factors
#'  countmat$strong = as.factor(countmat$strong)
#'  countmat$cancer_type = as.factor(countmat$cancer_type)
#'  countmat$sample_id = as.factor(countmat$sample_id)
#'  countmat$neighbors = as.factor(countmat$neighbors)
#'
#'  # set contrasts
#'  contrasts(countmat$cancer_type) = contr.sum(nlevels(countmat$cancer_type))
#'  nesting = nested_treatment_contrasts(outer.factor=countmat$cancer_type, inner.factor=countmat$sample_id)
#'  how.many = nlevels(countmat$sample_id) - nlevels(countmat$cancer_type)
#'  contrasts(countmat$sample_id, how.many=how.many) = nesting
#'
#'  # estimation
#'  fitlist[[i]] = fast_multinom(cbind(NO, I, VA, VG) ~ 0 + strong*sample_id + strong*replication_timing*cancer_type + strong*neighbors*cancer_type, data = countmat, refLevel=1, VC=T, loglik=T, subsetmatrix=NULL)
#' }
#'
#' # pooling
#' pooled.fit = pooling(fitlist, VC=T)
#'
#'
#'
#' @author Johanna Bertl

pooling <- function(reslist, VC=F){

  # check if the VC matrices are available for pooling
  if (VC) {
    testvec = sapply(reslist, function(x) is.null(x$VC))
    if (sum(testvec)>0) stop("The VC matrix cannot be pooled, because it is not available for all estimates.")
  }

  # pooling the means

  M = length(reslist[[1]]$coefficients)
  coef.list = vector("list", M)
  names(coef.list) = names(reslist[[1]]$coefficients)
  for(i in 1:M){
    coef.list[[i]] = rowMeans(coefficients_to_matrix(reslist, i))
  }

  # pooling the VC matrices (Schafer, p. 113 + 114)

  if (VC) {

    # number of multiple imputations
    mi = length(reslist)

    # Q_bar
    Q_bar = unlist(coef.list)
    k = length(Q_bar)

    # U_bar
    U_sum = reslist[[1]]$VC
    for(i in 2:mi) U_bar = U_sum + reslist[[i]]$VC
    U_bar = U_bar/mi

    # B
    Qmat = unlist(reslist[[1]]$coefficients)
    for(i in 2:mi){
      Qmat = cbind(Qmat, unlist(reslist[[i]]$coefficients))
    }
    # either B or solve(U_bar) below doesn't work. Here, I try to make sure that B works by implementing it using the Matrix package logic.
    Qmat = Matrix(Qmat)
    Qcent = Qmat - rowMeans(Qmat)
    # instead of: B = var(t(Qmat))
    B = tcrossprod(Qcent)/(mi-1)

    # gamma_bar (mean fraction of missing information, estimated under the assumption that gamma_bar = gamma_i, i=1, ..., k; see Li et al, 1991) and r_1
    lambda_bar = sum(diag(B%*%solve(U_bar)))/k
    gamma_bar = 1/(1 + lambda_bar)
    r_1 = (1 + 1/mi) * lambda_bar

    # pooled VC matrix (T_tilde)
    sigma = (1 - r_1)* U_bar

  } else {
    sigma = NULL
    gamma_bar = NULL
  }

  # pooling the likelihoods (if available)

  testvec2 = sapply(reslist, function(x) is.null(x$log.lik))
  if (sum(testvec2)==0) {
    log.lik.vec = numeric(length(reslist))
    for(i in 1:(length(reslist))) log.lik.vec[i] = reslist[[i]]$log.lik
    log.lik = log(mean(exp(log.lik.vec)))
  } else {log.lik=NULL}

  # output: object of class fast_multinom

  output = list(multinomformula = reslist[[1]]$multinomformula, formulae = reslist[[1]]$formulae, coefficients = coef.list, VC = sigma, gamma = gamma_bar, predictions = NULL, contrasts = reslist[[1]]$contrasts, xlevels = reslist[[1]]$xlevels, log.lik=log.lik, n = reslist[[1]]$n)

  class(output) = "fast_multinom"
  output

}


# list_of_results is a list that contains fast_multinom objects
# i is the index of the coefficient vector that should be put into a matrix using rbind

coefficients_to_matrix <- function(list_of_results, i){
  mat = numeric()
  for(j in 1:length(list_of_results)) mat = cbind(mat, list_of_results[[j]]$coefficients[[i]])
  mat
}

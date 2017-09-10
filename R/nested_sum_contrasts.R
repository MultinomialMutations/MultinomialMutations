#' Sum contrasts for a nested factor
#'
#' The function creates sum contrasts for a factor whose levels are nested
#' within the levels of another factor.
#'
#' The outer and inner factors need to have the same length. It can be data
#' that is used in a regression model, so it is allowed to have multiple
#' instances of the inner factor in the inner.factor, but it can also be
#' shorter factors that just describe the nesting.
#'
#' The contrast matrix can be used directly on the inner factor using the
#' funtion contrasts, but it is usually necessary to specify the argument
#' how.many=nlevels(inner.factor) - nlevels(outer.factor).
#'
#' The function \code{\link[psych]{superMatrix}} is used to make a block diagonal matrix.
#'
#'
#' @usage nested_sum_contrasts(outer.factor, inner.factor)
#'
#' @param outer.factor A factor.
#' @param inner.factor A factor whose levels are nested in the levels of the
#' outer.factor.
#'
#' @return An object of class \code{matrix} with \code{nlevels(inner.factor)}
#' rows and \code{nlevels(inner.factor) - nlevels(outer.factor)} columns.
#'
#' @note The nesting of the input vectors is not checked.
#'
#'
#' @author Johanna Bertl
#'
#' @seealso \code{\link{contr.sum}}, \code{\link{contrasts}}, \code{\link{nested_treatment_contrasts}}
#'
#' @keywords contrasts regression
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#' @examples
#'
#'
#' # Sum contrasts for the cancer type and the sample ID, which is nested in the cancer type.
#'
#' data("cancermutations")
#'
#' contrasts(cancermutations$cancer_type) = contr.sum(nlevels(cancermutations$cancer_type))
#' nesting = nested_sum_contrasts(outer.factor=cancermutations$cancer_type, inner.factor=cancermutations$sample_id)
#' how.many = nlevels(cancermutations$sample_id) - nlevels(cancermutations$cancer_type)
#' contrasts(cancermutations$sample_id, how.many=how.many) = nesting
#' 
#' # Visualize the contrast matrix: 
#' image(t(nesting[273:1,]))

nested_sum_contrasts <- function(outer.factor, inner.factor){

  # coerce to data.frame:
  dat = data.frame(outer.factor, inner.factor)
  # make rows of the data.frame unique, i. e. each row contains one inner.factor level and the corresponding outer.factor level
  nesting = unique(dat)
  # sort by outer.factor
  nesting = nesting[do.call(order, nesting),]

  # frequencies of the outer factor (sorted by the outer factor level)
  freq = as.vector(table(nesting$outer.factor))

  # sum contrasts for the inner factor: for each level of the outer factor, a sum contrast matrix with the number of occurrences of the level; combined in a block matrix
  nlev = length(freq)
  contrast.list = vector("list", nlev)
  for(i in 1:nlev) contrast.list[[i]] = contr.sum(freq[i])
  contrast.mat = superMatrix(contrast.list)

  # resorting such that the rows of the contrast matrix are sorted by the levels of the inner factor
  contrast.mat = contrast.mat[rank(nesting$inner.factor),]
  colnames(contrast.mat) = paste0("nsc_", levels(nesting$inner.factor)[rowSums(contrast.mat)>0]) ### nsc for nested sum contrasts

  return(contrast.mat)
}

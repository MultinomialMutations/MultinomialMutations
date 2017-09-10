#' Treatment contrasts for a nested factor
#'
#' The function creates treatment contrasts for a factor whose levels are nested
#' within the levels of another factor.
#'
#' The outer and inner factors need to have the same length. It can be the data
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
#' @usage nested_treatment_contrasts(outer.factor, inner.factor)
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
#' @seealso \code{\link{contr.treatment}}, \code{\link{contrasts}}, \code{\link{nested_sum_contrasts}}
#'
#' @keywords contrasts regression
#'
#' @examples
#'
#'
#' # Treatment contrasts for the cancer type and the sample ID, which is nested in the cancer type.
#'
#' data("cancermutations")
#' small = droplevels(cancermutations[cancermutations$cancer_type %in% c("KICH", "LGG"),])
#' 
#' nesting = nested_treatment_contrasts(outer.factor=small$cancer_type, inner.factor=small$sample_id)
#' how.many = nlevels(small$sample_id) - nlevels(small$cancer_type)
#' contrasts(small$sample_id, how.many=how.many) = nesting
#' 
#' # Visualize contrast matrix: 
#'  
#' image(t(nesting[18:1,]))
#' 
#' 
#' # Test estimation and prediction
#' 
#' fit = fast_multinom(cbind(NO, I, VA, VG) ~ cancer_type + sample_id, data = small, refLevel=1, VC=F, subsetmatrix=NULL)
#' 
#' newdata = small[sample.int(nrow(small), 100),] 
#' pred = predict.fast_multinom(fit, newdata)
#' 
#' # Comparison to estimation with sample_id only and standard treatment contrasts
#' 
#' small.alt = small
#' contrasts(small.alt$sample_id) = contr.treatment(nlevels(small.alt$sample_id))
#' newdata.alt = newdata
#' contrasts(newdata.alt$sample_id) = contr.treatment(nlevels(small.alt$sample_id))
#' 
#' fit2 = fast_multinom(cbind(NO, I, VA, VG) ~ sample_id, data = small.alt, refLevel=1, VC=F, subsetmatrix=NULL)
#' 
#' pred2 = predict.fast_multinom(fit2, newdata.alt)
#' 
#' max(abs(pred - pred2))
#' # Nearly the same. 

nested_treatment_contrasts <- function(outer.factor, inner.factor){

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
  for(i in 1:nlev) contrast.list[[i]] = contr.treatment(freq[i])
  contrast.mat = superMatrix(contrast.list)

  # resorting such that the rows of the contrast matrix are sorted by the levels of the inner factor
  contrast.mat = contrast.mat[rank(nesting$inner.factor),]
  colnames(contrast.mat) = paste0("ntc_", levels(nesting$inner.factor)[rowSums(contrast.mat)>0]) ### ntc for nested treatment contrasts

  return(contrast.mat)
}

#' Adding pseudo counts to the number of mutations
#'
#' The function \code{add_counts} adds pseudo counts to the number of mutations to incorporate prior information and avoid counts of zero.
#'
#' The function uses functionalities from the packages data.table.
#'
#' Let c_1, ..., c_n be a set of categorical variables and v_1, ..., v_m the remaining (categorical or continuous) explanatory variables in the dataset, except for the two variables sample_id and cancer_type.
#'
#' First, the function checks if there is a positive count for each mutation type, denoted by nI, nVA and nVG for each combination of c_1 x ... x c_n x sample_id. If there isn't, pseudo counts are added to nNO, nI, nVA and nVG.
#'
#' The pseudo counts are obtained from nNO, nI, nVA and nVG for each combination of c_1 x ... x c_n x v_1 x ... x x_m x cancer_type and added to the observed counts for each combination c_1 x ... x c_n x v_1 x ... x x_m x sample_id. The pseudo counts and the observed counts are weighted equally. The new sum nNO + nI + nVA + nVG is the same as originally, so the number of sites of a specific category is preserved and the size of the genome doesn't change.
#'
#' If 'make_integers=T is used, this only holds approximately, because after adjusting to the number of observed sites, the ceiling is used (to avoid zero counts). This avoids very small non-integer counts that can induce the same numerical problems as zero counts, but on the other hand it increases the number of counts and can cause quite substantial biases.
#'
#' Note that the number of mutations per sample is not preserved.
#'
#'
#' @param x data frame with a set of count columns named NO, I, VA, VG (and YES).
#' @param reference reference data frame of the same format. If reference=NULL, x is used. Using a different reference data frame than is currently not implemented.
#' @param categorical names of the categorical variables. This vector should not include cancer_type and sample_id.
#' @param make_integers logical. Should the final output contain integer counts only?
#'
#' @return A data frame (or data table) of the exact same format as the input table x with an additional logical column 'zero' (indicating the addition of pseudocounts because of a zero mutation count).
#'
#' @seealso \code{\link{add_counts_pres_mut}} -- does the same, but preserving the number of mutations per sample. 
#'
#' @author Johanna Bertl
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#' 
#' 
#' @examples
#' 
#' # Adding a prior to the example data
#' 
#' data(cancermutations)
#' newdata = add_counts(cancermutations, categorical=c("strong", "neighbors"), make.integers=T)
#' 
#' # Looking at a sample with few mutations to see the effect of the imputation: 
#' 
#' sample02 = cancermutations[cancermutations$sample_id=="GBM_TCGA_02_2483_01A" & cancermutations$strong==1 & cancermutations$neighbors=="TG",]
#' new02 = newdata[newdata$sample_id=="GBM_TCGA_02_2483_01A" & newdata$strong==1 & newdata$neighbors=="TG",]
#' 
#' # number of mutations before adding the prior:
#' sum(sample02$YES)
#' # number of mutations after adding the prior:
#' sum(new02$YES)
#' 


add_counts = function(x, reference = NULL, categorical, make.integers = F){

  if(!is.null(reference)){stop("Using a different reference data frame than 'x' is currently not implemented. Set 'reference=NULL'.")}
  if(any(c("cancer_type", "sample_id") %in% categorical)){stop("The vector 'categorical' should not include cancer_type and sample_id, because they are used in a different way.")}

  x = as.data.table(x)

  # set key to categorical variables
  setkeyv(x, c("sample_id", categorical))

  # compute sum of mutations for each category
  cat.sum = x[,list(sumI = sum(I), sumVA = sum(VA), sumVG = sum(VG)), by=c("sample_id", categorical)]
  # find out if the minimal count is zero
  cat.sum[, zero:= (sumI == 0) | (sumVA == 0) | (sumVG == 0)]
  # remove the sum columns
  cat.sum[, sumI:= NULL]
  cat.sum[, sumVA:= NULL]
  cat.sum[, sumVG:= NULL]

  # add the column zero to the data.table x
  x.new = x[cat.sum]

  # prepare the reference dataset
  all.names = names(x)[-which(names(x) %in% c("sample_id", "cancer_type", "NO", "I", "VA", "VG", "YES"))]
  setkeyv(x.new, c("cancer_type", all.names))
  reference = x.new[,list(sumNO = sum(NO), sumI = sum(I), sumVA = sum(VA), sumVG = sum(VG)), by=c("cancer_type", all.names)]
  # make the sum columns into numeric columns (to avoid integer overflows later)
  sumcols <- c("sumNO", "sumI", "sumVA", "sumVG")
  reference[, (sumcols) := lapply(.SD, as.numeric), .SDcols=sumcols]

  # join the reference dataset onto x.new
  x.joined = x.new[reference]

  ### adjust number of sites
  # add a column with the total number of sites of this category
  x.joined[,nsites := NO + I + VA + VG]
  x.joined[,nsites := as.numeric(nsites)]
  # add a column with the total number of sites of this category per cancer type
  x.joined[,sumsites:= sumI + sumNO + sumVA + sumVG]
  # adjust the sum columns such that the total site count is not more than within the sample
  x.joined[,sumI:= sumI*nsites/sumsites]
  x.joined[,sumNO:= sumNO*nsites/sumsites]
  x.joined[,sumVA:= sumVA*nsites/sumsites]
  x.joined[,sumVG:= sumVG*nsites/sumsites]

  # add the pseudo counts
  if(!make.integers){
    x.joined[,I:=as.numeric(I)]
    x.joined[,NO:=as.numeric(NO)]
    x.joined[,VA:=as.numeric(VA)]
    x.joined[,VG:=as.numeric(VG)]
    # this is necessary, because they were integer columns. Assigning them non-integer values just cuts what's after the comma (stupid data.table problem).
    x.joined[zero==T, I:= (I + sumI)/2]
    x.joined[zero==T, NO:= (NO + sumNO)/2]
    x.joined[zero==T, VA:= (VA + sumVA)/2]
    x.joined[zero==T, VG:= (VG + sumVG)/2]
  } else {
    x.joined[zero==T, I:= ceiling((I + sumI)/2)]
    x.joined[zero==T, NO:= ceiling((NO + sumNO)/2)]
    x.joined[zero==T, VA:= ceiling((VA + sumVA)/2)]
    x.joined[zero==T, VG:= ceiling((VG + sumVG)/2)]
  }

  # (dividing by 2 is not really necessary for the multinomial regression, but if the raw counts in the dataset are also used for some kinds of analyses, this makes sure that they are approximately the same as originally.)

  # remove superfluous columns
  x.joined[,c("sumNO", "sumI", "sumVA", "sumVG", "nsites", "sumsites"):=NULL]
  # add the column YES
  x.joined[,YES:=I + VA + VG]

  return(x.joined)

}



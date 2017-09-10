#' Wrapper for the preparation of mutation data for multinomial regression
#'
#' This function is a wrapper for the function \code{\link{count_table_prep_multinom}} that prepares the mutation data for multinomial (and binomial) regression. This wrapper was used in our analysis in Bertl et al. (2007) (see References).
#'
#' The scripts that were used to run this function and that show all settings used in Bertl et al. (2007) are available in this package in the folder \code{inst/Bertl_et_al_2017}. The pre-processed data can be downloaded from figshare.
#'
#'
#' @param cv integer. Which cross-validation piece should be handled?
#' @param m integer. Number of multiple imputation replicates. 0 means no imputation.
#' @param impute logical. Impute missing values or remove them?
#' @param datafolder character. The data is read from the file \code{paste0(datafolder, datafile, cv-1)}.
#' @param datafile character.
#' @param impfolder character. Folder where the imputed data should be saved.
#' @param data_source character. Either "fredriksson" or "pcawg".
#'
#' @return There is no return value. The reformatted and imputed data is saved along with the number of missing values and total mutation count.
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#'
#' @author Johanna Bertl
#'
#' @seealso \code{\link{count_table_prep_multinom}}, \code{\link{cancermutations}}
#'
wrapper_multinom_dataprep = function(cv, m, impute, datafolder, datafile, impfolder, data_source){

  # read file
  chrom21 = read.table(file=paste0(datafolder, datafile, cv-1), header = T, as.is=T)

  # data preparation, multiple imputation
  chrom21.imp = count_table_prep_multinom(chrom21, m, impute, DNase1_dummy=T, expression_dummy=T, data_source=data_source)

  # saving files
  system(paste0("mkdir -p ", datafolder, impfolder, cv))
  m = length(chrom21.imp[[1]])
  for(j in 1:m){
    # save imputed data to file
    write.table(chrom21.imp[[1]][[j]], file=paste0(datafolder, impfolder, cv, "/imp", j, ".txt"), quote=F, row.names=F)
    # save number of missing values
    missingval = chrom21.imp[[2]]
    total_count = chrom21.imp[[3]]
    prop_missing = missingval/total_count
    save(missingval, total_count, file = paste0(datafolder, impfolder, cv, "/missing_summary.Rdata"))
  }

}

#' Wrapper to add prior pseudo counts to the mutation data
#'
#' This function is a wrapper for the function \code{\link{add_counts}} that adds prior dirichlet counts to the mutation data for multinomial (and binomial) regression. This wrapper was used in our analysis in Bertl et al. (2007) (see References).
#'
#' The scripts that were used to run this function and that show all settings used in Bertl et al. (2007) are available in this package in the folder \code{inst/Bertl_et_al_2017}. The pre-processed data can be downloaded from figshare.
#'
#' @param cv integer. Number of the cross-validation piece that should be handled.
#' @param mi integer. Number of the multiple imputation replicate that should be handled.
#' @param datafolder character. The dataset from the file \code{paste0(datafolder, "/", cv, "/imp", mi, ".txt")} is used.
#' @param dirichletfolder character. Where to save the dataset after the dirichlet counts were added.
#' @param rm_folder character. Where to save the dataset after the dirichlet counts were added and the cancer types (below) removed.
#' @param rm_cancer character vector. Cancer types that should be deleted (because they contain too few mutations) (or NULL).
#'
#' @return There is no return value. The data is saved in the specified folder(s).
#'
#' @seealso \code{\link{add_counts}}
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#' @author Johanna Bertl
wrapper_prior_pseudo_counts = function(cv, mi, datafolder, dirichletfolder, rm_folder, rm_cancer=NULL){

  count = read.table(file=paste0(datafolder, "/", cv, "/imp", mi, ".txt"), header = T, as.is=T)

  new_count = add_counts(count, categorical=c("strong", "neighbors"), make.integers=F)

  write.table(new_count, file=paste0(dirichletfolder, "/", cv, "/imp", mi, ".txt"), quote=F, row.names=F)

  if(!is.null(rm_cancer)){

    rm_cancer = eval(parse(text = rm_cancer))

    dat = as.data.table(new_count)
    setkey(dat, cancer_type)
    newdat = subset(dat, !(cancer_type%in%rm_cancer))

    write.table(newdat, paste0(rm_folder, "/", cv, "/imp", mi, ".txt"), quote=F, row.names=F)

  }
}



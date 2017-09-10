#' Computing simple summary statistics on the original data
#'
#' \code{sumstats_orig} computes simple summary statistics on the original data (see example) and outputs them in a pdf file for a quick check. \code{sumstats_orig_knit_only} is for a two-stage procedure where the sumstats are computed first with \code{sumstats_orig} and the pdf file is created afterwards.
#'
#' R Markdown is used to create the pdf in \code{sumstats_orig}. The option \code{no_pdf_output=T} can be used to compute the summary statistics without creating the pdf. In that case, the output is saved in an Rdata file. This is useful when pandoc is not installed. Subsequently, the function \code{sumstats_orig_knit_only} can be used to knit the pdf document from the Rdata file.
#'
#' The set of columns that is always assumed to be part of the data set: sample_id, from, to, count, left, right. Additional explanatory variables can vary, so they have to be set by the parameter \code{explanatory}.
#'
#' @param inputfile Character. Original data file.
#' @param outputfile Character. Where to save the output (pdf file or Rdata file for \code{no_pdf_output=T}).
#' @param explanatory Character vector with the names of the additional explanatory variables. See Details for the set of default variables.
#' @param no_pdf_output Logical. Instead of the pdf output, the computed summary statistics are saved in an Rdata file.
#' @param sumstatsfile Rdata file that contains the sumstats that were produced with the option no_pdf_output=T.
#'
#' @return Nothing is returned. The output is a pdf or an Rdata file.
#'
#' @author Johanna Bertl
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#' @examples
#'
#' # use system.file to find the raw dataset that was installed along with the package:
#' location = system.file("extdata", "set0", package = "multinomutils")
#'
#' # immediately producing a pdf file
#' sumstats_orig(location, outputfile="~/sumstats.pdf", explanatory=c("phyloP", "replication_timing", "expression"))
#'
#' # producing an Rdata file with the summary statistics ...
#' sumstats_orig(location, outputfile="~/sumstats.Rdata", explanatory=c("phyloP", "replication_timing", "expression"), no_pdf_output=T)
#' # ... and the pdf file afterwards
#' sumstats_orig_knit_only("~/sumstats.Rdata", outputfile="~/sumstats_knit_only.pdf")

sumstats_orig = function(inputfile, outputfile, explanatory, no_pdf_output=F){

  # read file
  dat = read.table(file=inputfile, header = T, as.is=T) # change as.is=T to F?

  print(head(dat))

  if(!no_pdf_output){
    # knit the sum stats
    render(system.file("rmd/sumstats_original_data.Rmd", package="multinomutils"), output_file=outputfile)
  } else {
    sumstats_orig_no_pdf_output(dat, explanatory, outputfile)
  }

}


#' @rdname sumstats_orig

sumstats_orig_knit_only = function(sumstatsfile, outputfile){

  # open sumstatsfile
  load(sumstatsfile)

  # knit the sum stats
  render(system.file("rmd/sumstats_original_data_knit_only.Rmd", package="multinomutils"), output_file=outputfile)

}


# internal functions to compute summary statistics
sum.numeric = function(x) sum(as.numeric(x))
sum.na = function(x) sum(is.na(x))

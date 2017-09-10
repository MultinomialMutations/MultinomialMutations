#' Count table for multinomial regression
#'
#' \code{count_table_prep_multinom} imputes missing values and reshapes the raw count table for multinomial logistic regression on the strand-symmetric mutation model.
#'
#' If impute=T and m>=1, missing data is imputed m times. If impute=T and m=0, the missing data is not touched and just kept as NA. When impute=F, the value of m is irrelevant. In this case, only the complete cases are output (using the function \link[stats]{complete.cases}).
#'
#' The packages data.table and reshape2 are used for efficient and fast handling of the large mutation datasets. Multiple imputation is handled by the function \link{mimp}.
#'
#' Note that this function uses a few small functions from small_dataprep_functions.R.
#'
#' The Fredriksson and PCAWG datasets are handled in the same way, apart from the location information that is removed from the cancer type in the PCAWG set.
#'
#'
#' @param count_table data frame. Raw count table. See examples.
#' @param m integer. Number of multiple imputations. 0 means no imputation.
#' @param impute logical. Impute missing values or remove them? See details.
#' @param strong logical. Should the variable strong be computed?
#' @param CpG logical.
#' @param apobec logical.
#' @param neighbors logical.
#' @param DNase1_dummy logical. Should a dummy variable be computed for DNase1 peaks? If yes, NAs in the original variable DNase1 are replaced by zero. This means that DNase1 := DNase1*I(DNase1_dummy==1).
#' @param expression_dummy logical. Should a dummy variable be computed for expression measure available? If yes, NAs in the original variable expression are replaced by zero.
#' @param data_source character. Either "fredriksson" or "pcawg".
#'
#' @author Johanna Bertl & Malene Juul
#'
#' @seealso \code{\link{cancermutations}}
#'
#' @return A list that consists of the following elements:
#' \describe{
#'  \item{imputed}{a list of length min(m, 1) of imputed or complete data.tables (data.frames)}
#'  \item{missing}{a named vector giving the number of sites with missing values}
#'  \item{total_count}{an integer value giving the total number of sites (to compute proportions of missing sites)}
#' }
#'
#' @references Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 \url{http://www.biorxiv.org/content/early/2017/06/21/122879}
#'
#' @examples
#'
#' # This is how the example dataset cancermutations was created (data(cancermutations)).
#'
#' # use system.file to find the raw dataset that was installed along with the package:
#' location = system.file("extdata", "set0.Rdata", package = "multinomutils")
#' count.raw = load(file=location)
#'
#' # data preparation with imputation (on a subset of the data, for speed -- this still takes a few minutes!)
#' set.seed(1234)
#' count.raw.sub = count.raw[sample.int(nrow(count.raw), 1000),]
#' count.imp = count_table_prep_multinom(count.raw.sub, 2)
#' # Note that imputation only works if there is more than one cancer type in the dataset.
#'
#' # data preparation without imputation, but with expression dummy variable:
#' count.noimp = count_table_prep_multinom(count.raw, 0, expression_dummy=T)+
#'
#' # complete cases only
#' count.complete = count_table_prep_multinom(count.raw, m=5, impute=F)
#' # Note that this doesn't work with a very small subset of the data where after removal of the missing cases not all 4 mutation types exist.
#' # This is similar to the example dataset cancermutations. The code to create this dataset is in data-raw/cancermutations.R


count_table_prep_multinom = function(count_table, m, impute=T, strong=T, CpG=T, apobec=T, neighbors=T, DNase1_dummy=F, expression_dummy=F, data_source="fredriksson"){

  count_table = as.data.table(count_table)

  #### create variables for regression ####

  # sample IDs and cancer types:
  # replace "-" by "_" in the sample and cancer columns (to avoid problems in the formula), remove sampling location from the cancer type, add cancer name to sample (for easier interpretation of the regression coefficients)

  # create small table with the mapping btw cancer types and samples:
  setkey(count_table, sample_id, cancer_type)
  samples <- count_table[, .(s = 0), by = key(count_table)] ### can I remove the column s? It's entirely useless, but I didn't find a way to get this data frame without an aggregating function.
  samples[, s:= NULL]
  set.cancer.type.pcawg(dat = samples, cancer_type_column = "cancer_type", data_source=data_source)
  set.sample.id(dat = samples, sample_id_column = "sample_id", cancer_type_id_column = "cancer_type", data_source = data_source)

  # merge the new sample ids and cancer types onto the data set
  setnames(count_table, "sample_id", "sample_id_original")
  setkey(samples, sample_id_original)
  setkey(count_table, sample_id_original)
  count_table <- samples[count_table]
  count_table[,sample_id_original:=NULL]
  count_table[,cancer_type_original:=NULL]

  # explanatory variables:
  # 1) strong: G:C site
  if(strong){
    count_table[, strong := 0]
    setkey(count_table, from)
    count_table[c("C","G"), strong := 1]
  }
  # 2) CpG: (left C & ref G) | (ref C & right G)
  if(CpG){
    count_table[, CpG := 0]
    setkey(count_table, left, from)
    count_table[.("C","G"), CpG := 1]
    setkey(count_table, from, right)
    count_table[.("C","G"), CpG := 1]
  }
  # 3) apobec triplets: TCA, TCT, and reverse complements AGA and TGA
  if(apobec){
    setkey(count_table, left, from, right)
    count_table[, apobec:=0]
    count_table[.("T", "C", c("A", "T")), apobec:=1]
    count_table[.(c("A", "T"), "G", "A"), apobec:=1]
  }
  # 4) all left right combinations or their reverse complements (if there is an A or G at the center, the reverse complement is used)
  if(neighbors){
    left.right.neighbors(dataset = count_table, reference_column_id = "from", left_column_id = "left", right_column_id = "right")
  }
  # 5) DNase1 dummy
  if(DNase1_dummy){
    if(exists("DNase1", count_table)){
      count_table[, DNase1_peak:=as.numeric(!is.na(DNase1))]
      count_table[DNase1_peak==0, DNase1:=0]
    } else {
      warning("The DNase1 dummy could not be created, because the dataset doesn't contain DNase1.")
    }
  }
  # 6) expression dummy
  if(expression_dummy){
    if(exists("expression", count_table)){
      count_table[, expression_dummy:=as.numeric(!is.na(expression))]
      count_table[expression_dummy==0, expression:=0]
    } else {
      warning("The expression dummy could not be created, because the dataset doesn't contain expression.")
    }
  }


  # mutation types: no mutation, transition, 2 types of transversions
  # purines:
  set.multinom.mutation.type(dat = count_table)


  # remove columns from the original count table that are not used any more
  count_table[,c("left", "right", "ref", "chr", "from", "to"):=NULL]
  # aggregate categories that come up multiple times
  byvec = colnames(count_table)
  byvec = byvec[byvec!="count"]
  setkeyv(count_table, byvec) # using setkeyv for a character vector of column names
  count_table <- count_table[, .("count" = sum(count)), by=key(count_table)]


  #### imputation ####

  # number of missing values per column of count_table

  namevec = names(count_table)
  namevec_NA = paste0(namevec, "_NA")
  len = length(namevec)

  for(i in 1:len){
    count_table[, paste0(namevec[i], "_NA"):=is.na(get(namevec[i]))*count, with=F]
  }

  missing = numeric(len)
  for (i in 1:len){
    missing[i] = sum(count_table[, namevec_NA[i], with=F])
    count_table[,namevec_NA[i]:=NULL, with=F]
  }
  names(missing) = namevec[1:len]

  total.count = sum(as.numeric(count_table$count))


  # impute (or remove) missing values

  if(impute){

    if (m>0) {
      if(sum(missing)>0){
      # multiple imputation:
        count_table_imp = mimp(x = count_table, to.impute = names(missing)[which(missing>0)], m = m)
      } else {
        count_table_imp = list(count_table)
        m = 1
      }
    } else {
      # no imputation
      count_table_imp = list(count_table)
      m = 1
    }

  } else {
    # complete cases
    count_table_imp = list(count_table[complete.cases(count_table),])
    m = 1
  }


  #### reshape ####

  for(i in 1:m){

    # does this work if there are missing values???

    ######### ------------ change this -----------------------
    # reshape to have two count columns: no mutation and mutation
#    if (all(count_table_imp[[i]][,mut]==0)) {
#      setnames(count_table_imp[[i]], gsub(pattern="count", replace="nomut", x=names(out)))
#    } else {
##### --------------------------------------------------------
      new.cols = names(count_table_imp[[i]])
      new.cols.formula = formula(paste(paste(new.cols[!new.cols %in% c("count", "mut")], collapse=" + ")," ~ mut"))
      count_table_imp[[i]] = data.table(dcast(count_table_imp[[i]], new.cols.formula, value.var="count", fill=0, fun.aggregate=function(x) sum(as.numeric(x))))
      count_table_imp[[i]][, YES:= I + VG + VA]
#    }
  }

  # return results
  return(list(imputed = count_table_imp, missing = missing, total_count = total.count))
}


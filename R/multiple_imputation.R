#' Multiple imputation for count data
#'
#' The function \code{mimp} produces multiply imputed datasets for categorical count data.
#'
#' A count dataset consists of a set of categorical variables and a column of counts (named "count"). It is not necessary that each combination of categories is present in the dataset. The imputation is conducted for each variable independently by sampling from its empirical distribution function.
#'
#' Functionalities from the package \code{\link{data.table}} are used to speed up computations.
#'
#' There are sophisticated and flexible R packages for multiple imputation available. I wrote this function because the functions from package \code{cat}, e. g. \code{\link[cat]{em.cat}}, can't handle the large counts that we have in our dataset (it uses the integer class for communicating with an underlying Fortran function, so it can cause integer overflows) and the package \code{\link[mice]{mice}} doesn't support count data.
#'
#'
#' @param x data frame with ncols-1 categorical columns (do not need to be factors, but note that numeric vectors are handled as if they were categorical) and one column with counts, named "count".
#' @param to.impute character vector of names of the columns where the missing values should be imputed (if there are any, this is tested in the function)
#' @param m number of multiple imputations
#'
#' @return If there are no missing values in the specified columns, there is no output and the function throws a warning. Otherwise, a list of \code{m} complete data frames. Note that they have the same number of counts, but they will usually not have the same number of rows as the original data frame.
#'
#' @author Johanna Bertl
#'
#' @seealso \link[cat]{imp.cat}, \link[mice]{mice}
#'
#' @examples
#'
#' # use a subset of the raw dataset that was installed along with the package:
#' location = system.file("extdata", "set0", package = "multinomutils")
#' count_table = fread(file=location) 
#' count_table = count_table[count_table$cancer_type %in% c("KICH", "LGG")]
#' count_table[,left:=NULL]
#' count_table[,right:=NULL]
#' 
#' # number of missing values per column of count_table
#' missing = apply(count_table, 2, FUN = function(x) sum(is.na(x)))
#' names(missing) = names(count_table)
#'
#' # 3-fold multiple imputation
#' count_table_imp = mimp(x = count_table, to.impute = names(missing[missing>0]), m = 3)
#'

mimp = function(x, to.impute, m){

  x = as.data.table(x)
  if(!is.character(to.impute)) {stop("to.impute has to be a character vector.")}

  # prepare a column for each variable in to.impute with T/F for missing values, because this is needed later (if the overall column is only used for this check, maybe it can be skipped)
  setkeyv(x, cols=to.impute)
  to.impute.logical = paste0(to.impute, ".NA")  
  x[, (to.impute.logical) := F]
  for(i in 1:length(to.impute)){
    x[is.na(get(to.impute[i])), to.impute.logical[i] := T, with=F]
  }

  # index of incomplete cases:
  ### is apply a good idea?? --------------------------------------------
  x[, incomplete.index := apply(x[, to.impute.logical, with=F], 1, any)]
  ################ ------------------------------------------------------
  # are there any incomplete cases at all? If not, exit.
  incomplete.count = sum(x[,incomplete.index])
  if(incomplete.count==0) stop("No incomplete cases in any of the specified columns.")

  # list that contains m copies of the incomplete cases
  incomplete.list = vector("list", m)
  for(j in 1:m) incomplete.list[[j]] = subset(x, incomplete.index)

  # empty data structures:
  # list of indices of missing values, each entry of the list corresponds to one entry in to.impute
  missing.list = vector("list", length(to.impute))
  # vector of the number of missing values, each entry corresponds to one entry in to.impute
  #!# improve this with data.table functionality? ----------------
  missing.count = colSums(x[, to.impute.logical, with=F])
  # --------------------------------------------------------------

  # Multiple imputation variable by variable.
  # After imputation of a variable, the size of all data frames in incomplete.list can change: The number of combinations of categories changes and can be different for each of the m imputed sets. This requires the recomputation of the index of missing values for each variable in each imputed set.
  for(i in 1:length(to.impute)){

    # list of indices of missing values in the current variable (empty)
    missing.list[[i]] = vector("list", m)

    if(missing.count[i]>0){

      
      #### empirical frequencies of the observed values (exkl. NAs): ####
      
      setkeyv(x, to.impute.logical[i])
      ### -------------- see comment below --------------------------------------
      freq.raw = x[, sum(count), by=eval(to.impute[i])]
      freq.raw = freq.raw[as.vector(!is.na(freq.raw[,1, with=F]))]
      freq = unlist(freq.raw[, 2, with=F])
      values = unlist(freq.raw[, 1, with=F])
      ### originally, I used as.numeric, and that might have something to do with the NA category ...:
      # freq = tapply(as.numeric(x$count), x[,to.impute[i]], sum)
      ### It would probably be smarter to exclude the NA values when computing freq.raw, but I didn't manage. ------------------------------------------------------------

      ## NOTE: The complete cases are not needed anymore until merging with the imputed ones. If memory usage is a problem, they can be discarded (or saved to disc) now.

      # multiple imputation:
      for(j in 1:m){

        missing.list[[i]][[j]] = which(is.na(incomplete.list[[j]][,to.impute[i], with=F]))

        # vector of imputed values.
        # preliminary: matrix of imputed values. Note the format: columns - missing values, rows: counts of each category
        imputed.matrix = sapply(incomplete.list[[j]]$count[missing.list[[i]][[j]]], rmultinom.size, n=1, prob=freq, simplify=T)
        # change to vector
        imputed.vec = as.vector(imputed.matrix)

        # new matrix that contains in cols: 1) index of na value, 2) imputed value, 3) imputed count (zero counts are deleted)
        new.mat = data.table(index = rep(missing.list[[i]][[j]], each=length(values)), values = values, count = imputed.vec)
        new.mat = new.mat[new.mat[,count] !=0, ]

        # add an index to the incomplete data matrix
        incomplete.list[[j]][,index:=1:nrow(incomplete.list[[j]])]

        # join with the incomplete data matrix
        #setkey(incomplete.list[[j]], index)
        setkey(new.mat, index)
        setkey(incomplete.list[[j]], index)
        incomplete.list[[j]] = merge(incomplete.list[[j]], new.mat, by="index", all.x=T, suffixes=c("", ".y"), allow.cartesian=T)
        ### alternative (is supposed to be faster, but not as convenient, because it joins "the wrong way" and I don't have as much control on the column names):
        # incomplete.list[[j]] = new.mat[incomplete.list[[j]], allow.cartesian=T] 
        ### Is there a way to change the data.table incomplete.list[[j]] by reference?
        
        # Get the imputed values into the right columns
        # counts:
        incomplete.list[[j]][get(to.impute.logical[i]), count := count.y]
        incomplete.list[[j]][, count.y := NULL]
        # values:
        incomplete.list[[j]][get(to.impute.logical[i]), to.impute[i] := values, with=F ]
        incomplete.list[[j]][, values := NULL] 

      } # loop over j in 1:m

    } else {
      warning(paste0("There are no missing values in column ", to.impute[i], "."))
    }

  } # loop over i in 1:length(to.impute)

  # putting the complete dataset together: rbind of complete cases and imputed cases, then summing up over duplicated lines (lines that contain the same combinations of categories)

  # remove unneccessary columns from the whole data
  x[, to.impute.logical := NULL, with=F]
  
  for(j in 1:m){
    
    # remove unnecessary columns
    incomplete.list[[j]][, index:=NULL]
    incomplete.list[[j]][, to.impute.logical := NULL, with=F]
    
    # join with complete part of the data
    incomplete.list[[j]] = rbind(x[incomplete.index==F,], incomplete.list[[j]])
    
    # remove the last unnecessary column
    incomplete.list[[j]][, incomplete.index:=NULL]
    
    # aggregate rows with identical combinations of values
    byvec = colnames(incomplete.list[[j]])
    byvec = byvec[byvec!="count"]
    setkeyv(incomplete.list[[j]], byvec) # using setkeyv for a character vector of column names
    incomplete.list[[j]] <- incomplete.list[[j]][, .("count" = sum(count)), by=key(incomplete.list[[j]])]
    
    # old version using plyr: 
    # incomplete.list[[j]] = ddply(.data=incomplete.list[[j]], .variables = names(incomplete.list[[j]])[-countcol], .fun=sumfun)
  }

  incomplete.list
}


#####################################################################################
# convenience functions:

# union of a list of vectors: (l is a list)
unionSeveral <- function(l) { Reduce(union, l) }

# change order of arguments in rmultinom and take care of sizes that are larger than the largest possible integer (approx. 2*10^9, see help(integer), rmultinom uses integers):
rmultinom.size = function(size, n, prob){
  if(size>2*10^9){
    # the random number generation is divided into multiple runs of the function rmultinom each one with maximum size 2*10^9
    nruns = size %/% (2*10^9) + 1
    sizevec = c(rep(2*10^9, nruns-1), size %% (2*10^9))
    resmat = apply(as.matrix(sizevec), 1, rmultinom.size2, n=n, prob=prob)
    t(t(rowSums(resmat)))
  } else {
    rmultinom(n = n, size = size, prob = prob)
  }
}

# change order of arguments in rmultinom without taking care of the largest integer (used internally in rmultinom.size)
rmultinom.size2 = function(size, n, prob){
  rmultinom(n = n, size = size, prob = prob)
}

# sum over the counts (used in ddply)
# (avoids integer overflow, but I don't know if it can occur here)
sumfun = function(x) sum(as.numeric(x$count))

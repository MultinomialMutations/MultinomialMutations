# The two main functions (strong.variable and left.right.neighbors) do not return the data object, as they change them by reference...


# - function to create variable strong (reference = C/G)
strong.variable <- function(dataset, reference_column_id) {

  dataset[, strong := 0]
  setkeyv(dataset, reference_column_id)
  dataset[c("C","G"), strong := 1]

}


# - function to reverse complement a vector of bases
rev.compl <- function(base_vec) {

  out <- base_vec
  out[base_vec=="A"] <- "T"
  out[base_vec=="C"] <- "G"
  out[base_vec=="G"] <- "C"
  out[base_vec=="T"] <- "A"

  return(out)
}


# - function to create neighbor column (reverse complement if reference G/A)
left.right.neighbors <- function(dataset, reference_column_id, left_column_id, right_column_id) {

  # - set left and right reverse complement to the original left and right
  dataset[, left.rev.compl := get(left_column_id)]
  dataset[, right.rev.compl := get(right_column_id)]

  # - if the reference is G/A, change left_rc and right_rc to reverse complements
  setkeyv(dataset, reference_column_id)
  dataset[.(c("A", "G")), left.rev.compl := rev.compl(get(right_column_id))]
  dataset[.(c("A", "G")), right.rev.compl := rev.compl(get(left_column_id))]

  # - create neighbor variable
  dataset[, neighbors := paste0(left.rev.compl, right.rev.compl)]

  # - clean up
  dataset[, right.rev.compl := NULL][, left.rev.compl := NULL]

}

# function to reformat the cancer type for the pcawg dataset (remove the sequencing location)
set.cancer.type.pcawg = function(dat, cancer_type_column, data_source){

  if (data_source=="pcawg") {

    # keep original cancer type for mapping
    dat[, cancer_type_original:=get(cancer_type_column)]

    # remove sequencing location from cancer type
    dat[, cancer_type := apply(dat, 1, FUN = function(x) {
      gsub(pattern = "-..?", replacement = "", x = x[cancer_type_column])
    })]

  }

}

# - function to reformat the sample_id (replace "-" by "_" and add cancer type)
set.sample.id <- function(dat, sample_id_column, cancer_type_id_column, data_source) {

    # keep original sample id
    dat[, sample_id_original := get(sample_id_column)]

    # update sample id
    dat[, sample_id := apply(dat, 1, FUN = function(x) {
      paste(x[cancer_type_id_column],
            gsub(pattern = "-+", replacement = "_", x = strsplit(x[sample_id_column], "_")),
            sep="_")
    })]

}


# - function to set the mutation type
set.multinom.mutation.type <- function(dat) {

  # mutation types: no mutation, transition, 2 types of transversions

  # purines:
  dat[, pur_from := 0]
  dat[, pur_to:=0]

  dat[from %in% c("A", "G"), pur_from := 1]
  dat[to %in% c("A", "G"), pur_to := 1]

  # transitions:
  dat[, mut := "NO"]
  dat[from != to & pur_from == pur_to, mut := "I"]

  # transversions:
  dat[from != to & pur_from != pur_to & to %in% c("G", "C"), mut := "VG"] # to G or C
  dat[from != to & pur_from != pur_to & to %in% c("A", "T"), mut := "VA"] # to A or T

  # remove temporary columns pur_from and pur_to
  dat[, c("pur_from", "pur_to") := NULL]
}





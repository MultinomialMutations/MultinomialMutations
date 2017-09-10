sumstats_orig_no_pdf_output = function(dat, explanatory, outputfile){

  #### Basic sum stats ####

  # Total number of counts
  total = sum(as.numeric(dat$count))
  # Total number of mutations
  mut = sum(dat$count[dat$from != dat$to])
  # Mutation rate
  mut.rate = mut/total
  # number of samples
  nsamples = nlevels(as.factor(dat$sample_id))
  # number of genomic positions
  npos = total/nsamples
  # numbers of lines with missing values
  lines.missing = apply(dat, 2, sum.na)


  #### User specified explanatory variables ####

  nexp = length(explanatory)
  tab = vector("list", length=nexp)
  tab.prop = vector("list", length=nexp)
  missing.prop = vector("list", length=nexp)

  for(i in 1:nexp){

    # "explanatory variable:", var))
    # distribution: #
    tab[[i]] = unclass(by(dat$count, as.factor(dat[,explanatory[i]]), sum.numeric))
    tab.prop[[i]] = tab[[i]]/total

    # proportion missing: #
    missing.prop[[i]] = sum.numeric(dat$count[is.na(dat[,explanatory[i]])])/total

  }


  #### Default explanatory variables ####

  # left & right
  LR = unclass(by(dat$count, list(as.factor(dat$left), as.factor(dat$right)), sum.numeric))
  LR.prop = LR/total
  # from
  from = unclass(by(dat$count, as.factor(dat$from), sum.numeric))
  from.prop = from/total
  # to
  to = unclass(by(dat$count, as.factor(dat$to), sum.numeric))
  to.prop = to/total

  save(total, mut, mut.rate, nsamples, npos, lines.missing, explanatory, tab, tab.prop, missing.prop, LR, LR.prop, from, from.prop, to, to.prop, file=outputfile)

}

#' # Treatment contrasts for the cancer type and the sample ID, which is nested in the cancer type.
#'
#' data("cancermutations")
#' small = droplevels(cancermutations[cancermutations$cancer_type %in% c("KICH", "LGG"),])
#' 
#' nesting = nested_treatment_contrasts(outer.factor=small$cancer_type, inner.factor=small$sample_id)
#' how.many = nlevels(small$sample_id) - nlevels(small$cancer_type)
#' contrasts(small$sample_id, how.many=how.many) = nesting
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
#' 
#' 
#' # Test estimation and prediction with interactions
#' 
#' fit.inter = fast_multinom(cbind(NO, I, VA, VG) ~ cancer_type*replication_timing + sample_id*replication_timing, data = small, refLevel=1, VC=F, subsetmatrix=NULL)
#' 
#' pred.inter = predict.fast_multinom(fit.inter, newdata)
#' 
#' fit.inter2 = fast_multinom(cbind(NO, I, VA, VG) ~ sample_id*replication_timing, data = small.alt, refLevel=1, VC=F, subsetmatrix=NULL)
#' 
#' pred.inter2 = predict.fast_multinom(fit.inter2, newdata.alt)
#' 
#' max(abs(pred.inter - pred.inter2))
#' # Nearly the same. 
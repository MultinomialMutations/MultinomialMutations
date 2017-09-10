# This function is used by the function fast_multinom.R. 
# It is used to make a block diagonal matrix, using the function bdiag from package Matrix.

make.bdiag <-
function(X, indexvec){
  
  # get the indices where each block begins.
  indexmatrix = cbind(c(1, indexvec[-length(indexvec)]+1), indexvec)
  
  # extract the list of blocks
  nblocks = nrow(indexmatrix)
  blocklist = vector("list", nblocks)
  for(i in 1:nblocks){
    blocklist[[i]] = X[indexmatrix[i,1]:indexmatrix[i,2], indexmatrix[i,1]:indexmatrix[i,2]]
  }
  
  # make the block diagonal matrix (function bdiag from package Matrix)
  bdiag(blocklist)
}

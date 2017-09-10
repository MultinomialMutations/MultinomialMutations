\name{make.bdiag}
\alias{make.bdiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make a block diagonal matrix from a quadratic matrix.
}
\description{
The function makes a block diagonal matrix from quadratic blocks (potentially of different sizes) along the diagonal of a given matrix X. 
}
\usage{
make.bdiag(X, indexvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    A quadratic matrix.
}
  \item{indexvec}{
    A vector of indices that describe where the quadratic blocks end.
}
}
\details{
Block 1 is X[1:indexvec[1], 1:indexvec[1]], block 2 is X[(indexvec[1]+1):indexvec[2], (indexvec[1]+1):indexvec[2]], etc.
}
\value{
A sparse block diagonal matrix of class CsparseMatrix.
}
\author{
Johanna Bertl
}

\seealso{
\code{\link{bdiag}}, \code{\link{diag}} % and the function from the Psych package?
}
\examples{
X = matrix(1:100, 10, 10)
make.bdiag(X, c(5, 8, 10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ block matrix }
\keyword{ diagonal block matrix }
\keyword{ block diagonal matrix }
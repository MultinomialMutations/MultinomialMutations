\name{chrom21}
\alias{chrom21}
\docType{data}
\title{
Chromosome 21 for multinomial regression
}
\description{
Chromosome 21, preformatted for multinomial regression.
}
\usage{data("chrom21")}
\format{
  A data frame with 65290 observations on the following 12 variables.
  \describe{
    \item{\code{phyloP}}{a numeric vector}
    \item{\code{replication_timing}}{a numeric vector}
    \item{\code{fitCons}}{a numeric vector}
    \item{\code{sample_id}}{a factor with levels \code{sam_BLC_TCGA_BL_A13J_01A} \code{sam_BLC_TCGA_BT_A20P_01A} ... \code{sam_UCE_TCGA_EY_A1GW_01A}}
    \item{\code{cancer_type}}{a factor with levels \code{BLCA} \code{BRCA} \code{CRC} \code{GBM} \code{HNSC} \code{KICH} \code{KIRC} \code{LGG} \code{LUAD} \code{LUSC} \code{PRAD} \code{SKCM} \code{THCA} \code{UCEC}}
    \item{\code{strong}}{a numeric vector}
    \item{\code{CpG}}{a numeric vector}
    \item{\code{apobec}}{a numeric vector}
    \item{\code{I}}{a numeric vector}
    \item{\code{NO}}{a numeric vector}
    \item{\code{VA}}{a numeric vector}
    \item{\code{VG}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(chrom21)
## maybe str(chrom21.multinom) ; plot(chrom21.multinom) ...
}
\keyword{datasets}

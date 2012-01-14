\name{evora-package}
\alias{evora-package}
\alias{evora}
\docType{package}
\title{
Epigenetic Variable Outliers for Risk prediction Analysis (EVORA)
}
\description{
EVORA is a package aimed at DNA methylation studies of common diseases
which have been generated with Infinium 27k technology. It is an
adaptive index classification algorithm where a disease "risk" score is
assigned to a DNA methylation profile. Feature selection in EVORA
operates at the level of differential variability, selecting features
(e.g CpGs) according to differential variability and specifically
increased variability in the disease phenotype. The risk score is
computed over an optimal set of features and using an optimal threshold,
both estimated from an internal cross-validation in the training
set. Specifically, the risk score corresponds to the fraction of
selected features which have a methylation value above the
threshold. The threshold is estimated in the COPA (Cancer Outlier
Profile Analysis) basis, since this scale-invariant basis allows better identification of outlier samples.
}
\details{
\tabular{ll}{
Package: \tab evora\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-01-14\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
There are a total of three functions: \code{copa}, \code{DoEVORA} and
\code{PredEVORA} perform the COPA transformation, the training and
testing, respectively. Further details are described in reference.
}
\author{
Andrew E Teschendorff
Maintainer: <a.teschendorff@ucl.ac.uk>
}
\references{
Teschendorff AE, Widschwendter M. Differential variability improves the identification of cancer risk and early detection markers in DNA methylation studies. Submitted.
}
\keyword{multivariate}
\examples{
### See example for DoEVORA
}

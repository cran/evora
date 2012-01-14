\name{PredEVORA}
\alias{PredEVORA}
\title{
Testing EVORA on a test data set
}
\description{
This function uses the output from DoEVORA to predict risk scores in
  an independent test set.
}
\usage{
PredEVORA(evora.o, test.m, testCOPA.m, phenoTEST.v)
}

\arguments{
  \item{evora.o}{
The output from DoEVORA.
}
  \item{test.m}{
The test data set, with rows labeling features and
    columns samples with same identifiers as before.
}
  \item{testCOPA.m}{
The corresponding COPA transformed data matrix.
}
  \item{phenoTEST.v}{
A vector specifying the binary phenotype of
    interest of the test samples.
}
}

\value{
 \item{auc}{A vector giving the Area under the Curve (AUC) of the risk
    prediction, plus the 95 percent confidence interval.}
  \item{riskS}{The vector of EVORA risk scores, one entry for each
    sample in the test set.}
}
\references{
Teschendorff AE, Widschwendter M. Differential variability improves the identification of cancer risk and early detection markers in DNA methylation studies. Submitted.
}
\author{
Andrew E Teschendorff
}



\examples{
### See example for DoEVORA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}


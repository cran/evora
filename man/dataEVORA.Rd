\name{dataEVORA}
\alias{dataEVORA}
\docType{data}
\title{Data for evora package}
\description{
  Synthetic DNA methylation training and test data sets, both defined
  over 1000 features and 50 samples, in each case with a binary
  phenotype of interest. 25 samples represent "cancer" and the other 25
  a normal state. Of the 1000 features, 100 (the first 100 rows) exhibit increased
  variability in cancer due to methylation outliers, mimicking real
  data. Methylation outliers constitute 1/5 of all 25 "cancer" samples,
  i.e for each of the first 100 features, 5 of the 25 cancer samples are
  outliers, modelled as data from a beta distribution with shape
  parameters (5,5) representing a (hemi)methylated state. The rest of
  the data matrix is generated from a beta distribution with shape
  parameters (1,10), representing the unmethylated state.
}
\usage{data(dataEVORA)}
\format{
  This synthetic data set is a list object containing the following
  elements: \code{data}: data is the training data matrix (1000 features, 40 samples),
  \code{pheno}: is a binary phenotype vector for the training data,
  \code{test}: the test data matrix, \code{testpheno}: is the binary
  phenotype vector for the test data
}
\references{Differential variability improves the identification of
  cancer risk and early detection markers in DNA methylation studies. Submitted.}
\keyword{datasets}
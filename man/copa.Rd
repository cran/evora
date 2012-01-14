\name{copa}
\alias{copa}
\title{
Performs COPA transformation on given data matrix
}
\description{
Given a data matrix it performs the COPA transformation on it and returns
a ranked list of features (ranked according to COPA scores) at each
desired centile.
}
\usage{
copa(data.m, lowercent = 0.8, overexp.log = TRUE)
}

\arguments{
  \item{data.m}{
Data matrix. Rows label features, Columns samples.
}
  \item{lowercent}{
This fraction determines the maximum fraction of
    samples defining an outlier subgroup of interest. Specifically, this
    fraction is 1-\code{lowercent}, which by default is 0.2.
}
  \item{overexp.log}{
A logical variable, specifying the directionality
    of the outlier subgroups, i.e whether they exhibit higher or lower
    levels of the molecular quantity. By default this is higher/overexpressed.
}
}

\value{
  \item{copa}{The COPA transformed data matrix.}
  \item{centiles}{A vector specifying the centiles.}
  \item{out}{A list of matrices where each entry corresponds to one
    centile. Each matrix contains the ranked list of features, ranked
    according to copa-score.}
}
\references{
Tomlins SA et al. Recurrent fusion of TMPRSS2 and ETS transcription factor genes in prostate cancer. Science 2005 310(5748):644-8.
}
\author{
Andrew E Teschendorff
}
\examples{
## See example for DoEVORA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}


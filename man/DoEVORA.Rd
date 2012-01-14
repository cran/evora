\name{DoEVORA}
\alias{DoEVORA}
\title{
Training EVORA on data set
}
\description{
This function trains EVORA on a given data set. Briefly, it selects
features based on differential variability, i.e features for which the
variability in one phenotype is different from another. The current
implementation is aimed at cancer (disease) DNA methylation data and
therefore it selects for features that are hypervariable in the cancer
(disease) phenotype. The actual risk scores are computed using an
adaptive index algorithm over the selected features. The internal
parameters, including the number of features to include, are estimated
using an internal cross-validation.
}
\usage{
DoEVORA(train.m, copa.m, pheno.v, copath.v = 1:10,
 ntop.v = seq(100, 1000, 50), nup = 5000)
}
\arguments{
  \item{train.m}{
Data matrix for training: rows label features, columns label samples.
}
  \item{copa.m}{
The corresponding copa-transformed training data matrix,
  i.e the output from \code{copa}.
}
  \item{pheno.v}{
A vector specifying the binary phenotype of
    interest. It must be a vector of 1's and 2's, with 2 labeling the
    phenotype where the increased variability is hypothesized e.g cancer.
}
  \item{copath.v}{
A vector of thresholds on the copa scale to
    consider. By default this runs on an integer scale from 1 to 10,
    which is reasonable for DNA methylation data generated from the 27k
    platform.
}
  \item{ntop.v}{
This vector specifies the range and interval of top
    ranked features to consider in the internal cross-validation
    optimisation procedure.
}
  \item{nup}{
This number specifies the maximum number of top ranked
    features to consider when reranking features in the internal
    cross-validation. Lowering this number speeds the analysis although
    at the expense of increased bias. This number should be
    significantly larger than the maximum value in \code{ntop.v}.
}
}
\value{
 \item{data}{The training data matrix.}
  \item{riskCpG}{The optimal set of differentially variable features
    associated with phenotype of interest, ranked according to the
    differential variability statistic. What is returned is the index,
    i.e the row positions of data matrix labeling the features.}
  \item{th}{Optimal copa threshold to use in computing scores.}
  \item{ranked}{The ranked index vector of all features, ranked
    according to differential increased variability in the disease
    phenotype}
  \item{auc}{The matrix of Area Under the Curve (AUC) values from the
    internal cross-validation. Rows label the number of selected top
    ranked features. Columns labels the COPA threshold. The row and
    column corresponding to maximal AUC is the optimal choice.}
}
\references{
Teschendorff AE, Widschwendter M. Differential variability improves the identification of cancer risk and early detection markers in DNA methylation studies. Submitted.
}

\author{
Andrew E Teschendorff
}

\examples{
### Example
data(dataEVORA);
print(dim(dataEVORA$data));
### perform COPA transformation
copa.m <- copa(dataEVORA$data)$copa;

### Train EVORA
evora.o <- DoEVORA(dataEVORA$data,copa.m,dataEVORA$pheno,ntop.v=seq(50,500,50),
nup=nrow(dataEVORA$data));
print(paste("The optimal AUC in training set is = ",max(evora.o$auc),sep=""));

### How many of the differentially variable features are found?
### (the first 100 are truly differntially variable)
print(paste("Sensitivity=",length(intersect(evora.o$risk,1:100))/100,sep=""));
print(paste("PPV=",length(intersect(evora.o$risk,1:100))/length(evora.o$risk),
sep=""));

### Now test on independent data
copaTEST.m <- copa(dataEVORA$test)$copa;
pred.o <-
PredEVORA(evora.o,dataEVORA$test,copaTEST.m,dataEVORA$testpheno);
print(paste("Test set AUC=",pred.o$auc[2],sep=""));

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}


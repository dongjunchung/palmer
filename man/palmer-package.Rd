\name{palmer-package}
\alias{palmer-package}
\docType{package}

\title{
PALMER (Constrained Biclustering Algorithm to Improve Pathway Annotation Based on the Biomedical Literature Mining)
}

\description{
This package provides functions for fitting palmer,
a constrained biclustering approach that allows to identify indirect relationships among genes
based on the text mining of biomedical literature, which allows researchers to utilize prior biological knowledge to guide identification of gene-gene associations.
}

\details{
\tabular{ll}{
Package: \tab palmer\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2019-02-28\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package contains a main class, \code{palmer}, which represents PALMER model fit.

The main methods for the PALMER is \code{palmer}, which fits constrained biclustering.
}

\author{
Jin Hyun Nam

Maintainer: Dongjun Chung <chungd@musc.edu>
}

\references{
Nam JH, Couch D, Silveira W.A, Yu Z and Chung D (2019) ''PALMER: A constrained biclustering Algorithm to improve pathway annotation based on the biomedical literature mining''.
}

\keyword{ package }

\seealso{
\code{\link{palmer}}
}


\examples{
library(palmer)

# Loading data

data(sdata)
data(pathway)

# Fit PALMER

fit.palmer <- palmer(X=sdata,path=pathway,K=2,L=3,B=100)
fit.palmer

# Extract Gene and GO term cluster and its probability

predict(fit.palmer)

# Visualize the clustering analysis result

plot(fit.palmer)
}

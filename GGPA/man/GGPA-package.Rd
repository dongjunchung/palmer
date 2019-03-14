\name{GGPA-package}
\alias{GGPA-package}
\docType{package}
\title{
graph-GPA
}
\description{
This package provides functions for fitting graph-GPA,
a statistical framework to prioritize GWAS results by integrating pleiotropy.
}
\details{
\tabular{ll}{
Package: \tab GGPA\cr
Type: \tab Package\cr
Version: \tab 0.9.1\cr
Date: \tab 2016-07-31\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains a main class, \code{GGPA}, which represents graph-GPA model fit.
This package contains four main methods, 
\code{GGPA}, \code{assoc}, and \code{plot}.
\code{GGPA} method fits the graph-GPA model
and \code{assoc} method implements association mapping.
\code{plot} method provides a graph representing genetic relationship among phenotypes.
}
\author{
Hang J. Kim and Dongjun Chung

Maintainer: Hang J. Kim <hang.kim@uc.edu>, Dongjun Chung <dongjun.chung@gmail.com>
}
\references{
Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture."
}
\keyword{ package }
\seealso{
\code{\link{GGPA}}, \code{\link{assoc}}, \code{\link{plot}}, 
\code{\linkS4class{GGPA}}.
}
\examples{
library(GGPA)

# load simulation data

data(simulation)
		
# fit graph-GPA model

fit <- GGPA( simulation$pmat, nBurnin=200, nMain=200 )
fit

# association mapping for each phenotype

head(assoc( fit, FDR=0.1, fdrControl="global" ))

# hypothesis testing for 1st and 2nd phenotype pair

head(assoc( fit, FDR=0.1, fdrControl="global", i=1, j=2 ))

# plot phenotype graph

plot(fit)
}
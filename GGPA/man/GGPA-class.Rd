\name{GGPA-class}
\alias{GGPA-class}
\alias{show,GGPA-method}
\alias{plot,GGPA-method}
\alias{plot,GGPA,missing-method}
\alias{fdr,GGPA-method}
\alias{estimates,GGPA-method}
\alias{fdr}
\alias{estimates}
\Rdversion{1.1}
\docType{class}
\title{Class "GGPA" }
\description{
This class represents graph-GPA model fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GGPA", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"list"}, 
    	representing the MCMC draws. }
    \item{\code{summary}:}{Object of class \code{"list"}, 
    	representing the summary statistics. }
    \item{\code{setting}:}{Object of class \code{"list"}, 
    	representing the setting for graph-GPA model fitting. }
    \item{\code{gwasPval}:}{Object of class \code{"matrix"}, 
    		representing the p-value matrix from GWAS data. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "GGPA")}: provide brief summary of the object. }
    \item{plot}{\code{signature(x = "GGPA", y = "missing", pCutoff = 0.5, betaCI = 0.95)}: plot a phenotype graph.
    Nodes i and j are connected if the posterior probability of E_{ij} > \code{pCutoff} and the posterior probability of beta_{ij} > \code{betaCI}. }
    \item{fdr}{\code{signature(object = "GGPA", i=NULL, j=NULL)}: provide local FDR.
    By default (i.e., \code{i=NULL, j=NULL}),
    it returns a matrix of local FDR that a SNP is not associated with each phenotype (i.e., marginal FDR),
		where the order of columns is same as that in input GWAS data.
    If phenotype indices i and j are specified, a vector of corresponding local FDR is provided.
    }
    \item{estimates}{\code{signature(object = "GGPA")}: 
        extract parameter estimates from graph-GPA model fit. }
     }
}
\references{
Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture."
}
\author{ Hang J. Kim, Dongjun Chung }
\seealso{
\code{\link{GGPA}}.
}
\examples{
showClass("GGPA")

\dontrun{		
library(GGPA)
data(simulation)
fit <- GGPA( simulation$pmat )
fit
plot(fit)
head(fdr( fit ))
head(fdr( fit, i=1, j=2 ))
str(estimates( fit ))
}
}
\keyword{classes}

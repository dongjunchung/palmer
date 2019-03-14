\name{GGPA}
\alias{GGPA}
\title{
Fit graph-GPA model
}
\description{
Fit graph-GPA model.
}
\usage{
GGPA( gwasPval, nBurnin=10000, nMain=40000, seedNum=1234, lbPval=1e-10, verbose=1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gwasPval}{p-value matrix from GWAS data,
		where row and column correspond to SNP and phenotype, respectively. }
    \item{nBurnin}{Number of burn-in iterations. Default is 10000.}
    \item{nMain}{Number of main MCMC iterations. Default is 40000.}
    \item{seedNum}{Random seed number. Default is 1234.}
  	\item{lbPval}{Lower bound for GWAS p-value. 
  		Any GWAS p-values smaller than \code{lbPval} are set to \code{lbPval}.
  		Default is 1e-10.}  
  	\item{verbose}{Amount of progress report during the fitting procedure. 
  		Possible values are 0 (minimal output), 1, 2, or 3 (maximal output).
  		Default is 1.} 
}
\details{
\code{GGPA} fits the graph-GPA model. It requires to provide GWAS p-value to \code{gwasPval}. Based on this GGPA fit, \code{assoc} implements association mapping and \code{plot} provides a phenotype graph.
}
\value{
Construct \code{GGPA} class object.
}
\references{
Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture."
}
\author{
Hang J. Kim and Dongjun Chung
}
\seealso{
\code{\link{assoc}}, \code{\linkS4class{GGPA}}.
}
\examples{
\dontrun{
library(GGPA)
data(simulation)
fit <- GGPA( simulation$pmat )
fit
head(assoc( fit, FDR=0.1, fdrControl="global" ))
head(assoc( fit, FDR=0.1, fdrControl="global", i=1, j=2 ))
plot(fit)
}
}
\keyword{models}
\keyword{methods}

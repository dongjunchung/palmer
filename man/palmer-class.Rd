\name{palmer-class}
\alias{palmer-class}
\docType{class}

\title{Class "palmer"}

\description{
This class represents palmer fit.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("palmer", ...)}.
}

\section{Slots}{
\describe{

\item{\code{data}}{Object of class \code{"matrix"}, 
    	representing the raw data for PALMER model. }

\item{\code{init}}{Object of class \code{"list"}, 
    	representing the setting for PALMER model fitting.} 

\item{\code{result}}{Object of class \code{"list"}, 
    	representing the fitted PALMER model.} 
}}


\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "palmer")}: provide brief summary of the object. }
    
    \item{predict}{\code{signature(object = "palmer")}: 
		provide Gene and GO term cluster and its probability.}
		
    \item{plot}{\code{signature(x = "plamer", y=missing)}: 
    provide visualization of constrained biclustersing analysis result.}

}}

\details{

}

\references{

}

\author{ Jin Hyun Nam }

\seealso{
\code{\link{palmer}}.
}

\examples{
showClass("palmer")

data(sdata)
data(pathway)

fit.palmer <- palmer(X=sdata,path=pathway,K=2,L=3,B=100)
fit.palmer
predict(fit.palmer)
plot(fit.palmer)
}

\keyword{classes}
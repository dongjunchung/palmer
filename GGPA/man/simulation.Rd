\name{simulation}
\docType{data}
\alias{simulation}
\title{Simulation dataa for graph-GPA}
\description{
  This is an simulation dataset.
}
\usage{ data(simulation) }
\format{
  \code{simulation} list object containing simulation data (element Y) and its simulation setting (the remaining elements).
}
\examples{
\dontrun{
data(simulation)
head(t(simulation$pmat))
}
}
\keyword{datasets}

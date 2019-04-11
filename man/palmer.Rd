\name{palmer}
\alias{palmer}

\title{PALMER: A Constrained Biclustering Algorithm to Improve Pathway Annotation Based on the Biomedical Literature Mining}

\description{
Fit PALMER model
}

\usage{
palmer(X, path = NULL, K = 3, L = 3, B = 1000)
}

\arguments{
\item{X}{Matrix of binary values with n samples by p variables}

\item{path}{List. Gene list for the prior information (default: NULL)}

\item{K}{The number of gene clusters (default: 3)}

\item{L}{The number of GO term clusters (default: 3)}

\item{B}{The number of bootstrappings (default: 1000)}
}

\details{
\code{palmer} fits the PALMER model. It requires binary matrix to be clusterd and 
list of gene set for the constraints.
}

\value{
Construct \code{palmer} class object.
}

\references{
...
}

\author{ Jin Hyun Nam }

\examples{
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

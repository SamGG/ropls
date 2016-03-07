\name{getWeightMN}

\alias{getWeightMN}
\alias{getWeightMN,opls-method}

\title{
    getWeightMN method for (O)PLS(-DA) models
}

\description{
  (Orthogonal) weights of the (O)PLS(-DA) model
}

\usage{
\S4method{getWeightMN}{opls}(object, orthoL = FALSE)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
  \item{orthoL}{
   Logical: Should the orthogonal weight matrix be returned (default is
  FALSE and the predictive weight matrix is returned)
  }
}

\value{
  Numeric matrix with a number of rows equal to the number of variables
  and a number of columns equal to the number of components
}

\author{
  Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)

sacurine.plsda <- opls(dataMatrix,
                       sampleMetadata[, "gender"])

getWeightMN(sacurine.plsda)

detach(sacurine)

}
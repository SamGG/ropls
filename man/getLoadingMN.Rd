\name{getLoadingMN}

\alias{getLoadingMN}
\alias{getLoadingMN,opls-method}

\title{
    getLoadingMN method for PCA/(O)PLS(-DA) models
}

\description{
  (Orthogonal) loadings of the PCA/(O)PLS(-DA) model
}

\usage{
\S4method{getLoadingMN}{opls}(object, orthoL = FALSE)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
  \item{orthoL}{
   Logical: Should the orthogonal loading matrix be returned (default is
  FALSE and the predictive loading matrix is returned)
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

getLoadingMN(sacurine.plsda)

detach(sacurine)

}
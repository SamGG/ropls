\name{getScoreMN}

\alias{getScoreMN}
\alias{getScoreMN,opls-method}

\title{
    getScoreMN method for PCA/(O)PLS(-DA) models
}

\description{
  (Orthogonal) scores of the (O)PLS(-DA) model
}

\usage{
\S4method{getScoreMN}{opls}(object, orthoL = FALSE)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
  \item{orthoL}{
   Logical: Should the orthogonal score matrix be returned (default is
  FALSE and the predictive score matrix is returned)
  }
}

\value{
  Numeric matrix with a number of rows equal to the number of samples
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

getScoreMN(sacurine.plsda)

detach(sacurine)

}
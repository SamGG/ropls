\name{residuals.opls}

\alias{residuals.opls}
\alias{residuals,opls-method}

\title{
    Residuals method for (O)PLS models
}

\description{
  Returns the residuals from the (O)PLS(-DA) regression models
}

\usage{
\S4method{residuals}{opls}(object, ...)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
  \item{\dots}{
    Currently not used.
  }
}

\value{
  Numeric matrix or vector (same dimensions as the modeled y response);
  if y is a character vector or a factor (in case of classification),
  the residuals equal 0 (predicted class identical to the true
  class) or 1 (prediction error)
}

\author{
  Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)

sacurine.pls <- opls(dataMatrix,
                     sampleMetadata[, "age"])

head(residuals(sacurine.pls))

detach(sacurine)

}
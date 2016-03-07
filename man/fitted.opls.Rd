\name{fitted.opls}

\alias{fitted.opls}
\alias{fitted,opls-method}

\title{Fitted method for 'opls' objects}

\description{
  Returns predictions of the (O)PLS(-DA) model on the training dataset
}

\usage{
\S4method{fitted}{opls}(object, ...)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by the \code{opls} function.
  }
  \item{\dots}{
    Currently not used.
  }
}

\value{
  Predictions (either a vector, factor, or matrix depending on the y
  response used for training the model)
}

\author{
  Philippe Rinaudo and Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)
sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])

fitted(sacurine.plsda)

detach(sacurine)

}
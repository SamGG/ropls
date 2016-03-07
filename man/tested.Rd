\name{tested}

\alias{tested}
\alias{tested,opls-method}

\title{
    Tested method for (O)PLS models
}

\description{
  Returns predictions of the (O)PLS(-DA) model on the out of the box
  samples (when a 'subset' of samples has been selected when training
  the model)
}

\usage{
\S4method{tested}{opls}(object)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
}

\value{
  Predictions (either a vector, factor, or matrix depending on the y
  response used for training the model)
}

\author{
  Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)

testedorMN <- dataMatrix
responseFc <- sampleMetadata[, "gender"]

sacurine.plsda <- opls(testedorMN,
                       responseFc,
                       subset = "odd")

trainVi <- getSubsetVi(sacurine.plsda)

table(responseFc[trainVi], fitted(sacurine.plsda))

detach(sacurine)

}
\name{predict.opls}

\alias{predict.opls}
\alias{predict,opls-method}

\title{
    Predict method for (O)PLS models
}

\description{
  Returns predictions of the (O)PLS(-DA) model on a new dataset
}

\usage{
\S4method{predict}{opls}(object, newdata, ...)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
  \item{newdata}{
    Either a data frame or a matrix, containing numeric columns only,
    with the same number of columns (variables) as the 'x' used
    for model training with 'opls'.
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
  Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)

predictorMN <- dataMatrix
responseFc <- sampleMetadata[, "gender"]

sacurine.plsda <- opls(predictorMN,
                       responseFc,
                       subset = "odd")

trainVi <- getSubsetVi(sacurine.plsda)

table(responseFc[trainVi], fitted(sacurine.plsda))

table(responseFc[-trainVi],
      predict(sacurine.plsda, predictorMN[-trainVi, ]))

detach(sacurine)

}
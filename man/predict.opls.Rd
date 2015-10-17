\name{predict.opls}

\alias{predict.opls}

\title{
    Predict Method for (O)PLS
}

\description{
   This function predicts values based upon a model trained by \code{opls}.
}

\usage{
\method{predict}{opls}(object, newdata, ...)
}

\arguments{
  \item{object}{
    Object of class \code{"opls"}, created by \code{opls}.
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

trainVi <- sacurine.plsda[["subset"]]

table(responseFc[trainVi], fitted(sacurine.plsda))

table(responseFc[-trainVi],
      predict(sacurine.plsda, predictorMN[-trainVi, ]))

detach(sacurine)

}
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
    newdata
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
responseFc <- sampleMetadata[, "genderFc"]

sacurine.plsda <- opls(predictorMN,
                       responseFc,
                       subset = "odd")

trainVi <- sacurine.plsda[["subset"]]

table(responseFc[trainVi], fitted(sacurine.plsda))

table(responseFc[-trainVi],
      predict(sacurine.plsda, predictorMN[-trainVi, ]))

detach(sacurine)

}
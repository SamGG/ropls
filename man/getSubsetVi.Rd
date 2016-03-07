\name{getSubsetVi}

\alias{getSubsetVi}
\alias{getSubsetVi,opls-method}

\title{
    getSubsetVi method for (O)PLS(-DA) models
}

\description{
  Extracts the indices of the samples used for building the model (when
  a subset argument has been specified)
}

\usage{
\S4method{getSubsetVi}{opls}(object)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
}

\value{
  Integer vector with the indices of the samples used for training
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

detach(sacurine)


}
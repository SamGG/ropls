\name{getSummaryDF}

\alias{getSummaryDF}
\alias{getSummaryDF,opls-method}

\title{
    getSummaryDF method for PCA/(O)PLS models
}

\description{
  Summary of model metrics
}

\usage{
\S4method{getSummaryDF}{opls}(object)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
}

\value{
  Data frame
}

\author{
  Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)

sacurine.plsda <- opls(dataMatrix,
                       sampleMetadata[, "gender"])

getSummaryDF(sacurine.plsda)

detach(sacurine)

}
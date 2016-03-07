\name{print.opls}

\alias{print.opls}
\alias{print,opls-method}

\title{Print method for 'opls' objects}

\description{
  Displays information about the dataset and the model.
}

\usage{
\S4method{print}{opls}(x, ...)
}

\arguments{
  \item{x}{
    An S4 object of class \code{opls}, created by the \code{opls} function.
  }
  \item{\dots}{
    Currently not used.
  }
}

\value{
  Invisible.
}

\author{
  Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)
sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])

print(sacurine.plsda)

detach(sacurine)

}
\name{coef.opls}

\alias{coef.opls}
\alias{coef,opls-method}

\title{
    Coefficients method for (O)PLS models
}

\description{
  Coefficients of the (O)PLS(-DA) regression model
}

\usage{
\S4method{coef}{opls}(object, ...)
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
  Numeric matrix of coefficients (number of rows equals the number of
  variables, and the number of columns equals the number of responses)
}

\author{
  Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)

sacurine.plsda <- opls(dataMatrix,
                       sampleMetadata[, "gender"])

head(coef(sacurine.plsda))

detach(sacurine)

}
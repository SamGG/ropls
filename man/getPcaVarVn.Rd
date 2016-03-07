\name{getPcaVarVn}

\alias{getPcaVarVn}
\alias{getPcaVarVn,opls-method}

\title{
    getPcaVarVn method for PCA models
}

\description{
  Variance of the components (score vectors)
}

\usage{
\S4method{getPcaVarVn}{opls}(object)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by \code{opls} function.
  }
}

\value{
  Numeric vector with the same length as the number of components
}

\author{
  Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)

sacurine.pca <- opls(dataMatrix)

getPcaVarVn(sacurine.pca)

detach(sacurine)

}
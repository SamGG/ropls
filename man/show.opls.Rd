\name{show.opls}

\alias{show.opls}
\alias{show,opls-method}

\title{Show method for 'opls' objects}

\description{
  Displays information about the dataset and the model.
}

\usage{
\S4method{show}{opls}(object)
}

\arguments{
  \item{object}{
    An S4 object of class \code{opls}, created by the \code{opls} function.
  }
}

\value{
  Invisible.
}

\author{
  Philippe Rinaudo and Etienne Thevenot (CEA)
}

\examples{

data(sacurine)
attach(sacurine)
sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])

show(sacurine.plsda)

detach(sacurine)

}
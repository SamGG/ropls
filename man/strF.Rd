\name{strF}

\alias{strF}

\title{
  Printed summary of an R object
}

\description{
  Displays the class, mode, size and first...last values of the object

}

\usage{
strF(inpMF, borN = 2, mrkC = ",")
}

\arguments{
  \item{inpMF}{
    Input matrix, dataframe or vector
  }
  \item{borN}{
    Number of border (first and last) rows and columns to display
  }
  \item{mrkC}{
    Big mark separator for summary results
  }
}

\value{
  This function has no output.
}

\author{
Etienne Thevenot (CEA)
}

\seealso{
  \code{\link{str}}
}

\examples{
data(sacurine)
strF(sacurine[['dataMatrix']])
strF(sacurine[['sampleMetadata']])
}


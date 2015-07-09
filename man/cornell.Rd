\name{cornell}

\alias{cornell}

\docType{data}

\title{
  Octane of various blends of gasoline
}

\description{
  Twelve mixture component proportions of the blend are analysed
}

\usage{
data(cornell)
}

\format{
  A data frame with the following parameters:
  \describe{
    \item{\code{num}}{
      mixture number
    }
    \item{\code{x1}}{
      proportion of component 1
    }
    \item{\code{x2}}{
      proportion of component 2
    }
    \item{\code{x3}}{
      proportion of component 3
    }
    \item{\code{x4}}{
      proportion of component 4
    }
    \item{\code{x5}}{
      proportion of component 5
    }
    \item{\code{x6}}{
      proportion of component 6
    }
    \item{\code{x7}}{
      proportion of component 7
      Note: the 7 variables are correlated since they sum up to 1
    }
    \item{\code{y}}{
      octane (quantitative variable)
    }
  }
}

\value{
  Data frame (numeric type only; the first column can be transformed
  into row names) with 12 rows and 9 columns
  corresponding to the 'num'ber of the mixture (column 1), the
  proportion of each of the 7 'x' components within the mixture (columns
  2-8), and the octane indice 'y' (column 9). For details see the 'Format'
  section above.
}

\source{
Tenenhaus (1998), Table 6, page 78.
}

\references{
Tenenhaus (1998). La regression PLS: theorie et pratique. Paris: Editions Technip.
}

\keyword{datasets}

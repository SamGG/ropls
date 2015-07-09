\name{mark}

\alias{mark}

\docType{data}

\title{'mark' Dataset}

\description{
Examination marks obtained by French students in Mathematics, Physics,
French and English
}

\usage{
data(mark)
}

\format{
  A data frame with the following parameters:
  \describe{
    \item{\code{nom}}{
      names of the students
    }
    \item{\code{math}}{
      marks in mathematics
    }
    \item{\code{phys}}{
      marks in physics
    }
    \item{\code{fran}}{
      marks in french
    }
    \item{\code{angl}}{
      marks in english
    }
  }
}

\value{
  Data frame (numeric type except the first column, which can be
  transformed into row names) with 9 rows and 5 columns, corresponding
  to the name of the students (column 1), followed by the marks
  obtained in Maths, Physics, French and English (columns 2-5). For
  details see the 'Format' section above.
}

\source{
'mark' dataset.
}

\references{
Baccini (2010). Statistique Descriptive Multidimensionnelle (pour les nuls).
}

\keyword{datasets}

\name{linnerud}

\alias{linnerud}

\docType{data}

\title{Linnerud Dataset}

\description{
Three physiological and three exercise variables are measured on twenty middle-aged men in a fitness club.
}

\usage{
data(linnerud)
}

\format{
  A data frame with the following parameters:
  \describe{
    \item{\code{num}}{
      subject number
    }
    \item{\code{weight}}{
      weight
    }
    \item{\code{waist}}{
      waist
    }
    \item{\code{pulse}}{
      pulse
    }
    \item{\code{pullUp}}{
      pull-up
    }
    \item{\code{squat}}{
      situp
    }
    \item{\code{jump}}{
      jump
    }
  }
}

\value{
  Data frame (numeric type only; the first column can be transformed
  into row names) with 20 rows and 7 columns
  corresponding to the subject's 'num'ber (column 1), the 3
  physiological variables (columns 2-4), and the 3 exercise variables
  (columns 5-7). For details see the 'Format' section above.
}

\source{
'mixOmics' 'linnerud' dataset.
}

\references{
Tenenhaus (1998). La regression PLS: theorie et pratique. Paris: Editions Technip.
}

\keyword{datasets}

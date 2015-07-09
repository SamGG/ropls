\name{lowarp}

\alias{lowarp}

\docType{data}

\title{
  A multi response optimization data set (LOWARP)
}

\description{
  This example concerns the development of a polymer similar to that
  used in the plastic covering of mobile phones. The desired profile of
  the polymer was low warp and high strength. Four constituents (glas,
  crtp, mica, and amtp) were varied in the polymer formulation by means
  of a 17 run mixture design. For each new polymer, i.e., each new
  experiment in the mixture design, 14 responses relating to both warp
  and strength were measured on the product. The objective of the data
  analysis was to uncover which combination of factors (the four
  ingredients) gave polymers with low warp and high strength. The data
  set contains 10 missing values (NA).
}

\usage{
data(lowarp)
}

\format{
  A data frame with the following parameters:
  \describe{
    \item{\code{num}}{
      mixture number
    }
    \item{\code{glas}}{
      glas constituent
    }
    \item{\code{crtp}}{
      crtp constituent
    }
    \item{\code{mica}}{
      mica constituent
    }
    \item{\code{amtp}}{
      amtp constituent
    }
    \item{\code{wrp1}}{
      warp response 1
    }
    \item{\code{wrp2}}{
      warp response 2
    }
    \item{\code{wrp3}}{
      warp response 3
    }
    \item{\code{wrp4}}{
      warp response 4
    }
    \item{\code{wrp5}}{
      warp response 5
    }
    \item{\code{wrp6}}{
      warp response 6
    }
    \item{\code{wrp7}}{
      warp response 7
    }
    \item{\code{wrp8}}{
      warp response 8
    }
    \item{\code{st1}}{
      strength response 1
    }
    \item{\code{st2}}{
      strength response 2
    }
    \item{\code{st3}}{
      strength response 3
    }
    \item{\code{st4}}{
      strength response 4
    }
    \item{\code{st5}}{
      strength response 5
    }
    \item{\code{st6}}{
      strength response 6
    }    
  }
}

\value{
  Data frame (numeric type only; the first column can be transformed
  into row names) with 17 rows and 19 columns
  corresponding to the subject's 'num'ber (column 1), the 4 constituent
  variables (columns 2-5), the 8 warp responses (columns 6-13), and the
  6 strength responses (columns 14-19). For details see the 'Format'
  section above.
}

\references{
Eriksson et al. (2006). Multi- and Megarvariate Data Analysis. Umetrics
Academy. pp.16, 77, 209.
}

\keyword{datasets}

\name{aminoacids}

\alias{aminoacids}

\docType{data}

\title{Amino-Acids Dataset}

\description{
Quantitative structure property relationship (QSPR)
}

\usage{
data(aminoacids)
}

\format{
  A data frame with the following parameters:
  \describe{
    \item{\code{AA}}{
      amino acid
    }
    \item{\code{PIE}}{
      lipophilicity constant of the AA side chain
    }
    \item{\code{PIF}}{
      lipophilicity constant of the AA side chain
    }
    \item{\code{DGR}}{
      free energy of transfer of an AA side chain from protein interior
      to water
    }
    \item{\code{SAC}}{
      water-accessible surface area of AA's calculated by MOLSV
    }
    \item{\code{MR}}{
      molecular refractivity
    }
    \item{\code{Lam}}{
      polarity parameter
    }
    \item{\code{Vol}}{
      molecular volume of AA's calculated by MOLSV
    }
    \item{\code{DDGTS}}{
      free energy of unfolding of the tryptophane synthase a unit of
      bacteriophage T4 lysosome
    }
  }
}

\value{
  Data frame (numeric type except the first column, which can be
  transformed into row names) with 19 rows and the 9 columns contaning
  information about amino acids. For details see the 'Format' section above.
}

\source{
'aminoacids' dataset.
}

\references{
Wold et al. (2001). PLS-regression: a basic tool of
chemometrics. Chemometrics and Intelligent Laboratory Systems. 58:109-130.
}

\keyword{datasets}

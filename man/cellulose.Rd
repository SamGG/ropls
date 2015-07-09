\name{cellulose}

\alias{cellulose}

\docType{data}

\title{
  NIR-Viscosity example data set to illustrate multivariate calibration
  using PLS, spectral filtering and OPLS
}

\description{
  The data were collected at Akzo Nobel, Ornkoldsvik (Sweden). The raw
  material for their cellulose derivative process is delivered to the
  factory in form of cellulose sheets. Before entering the process the
  cellulose sheets are controlled by a viscosity measurement, which
  functions as a steering parameter for that particular batch.
  In this data set NIR spectra for 180 cellulose sheets were collected
  after the sheets had been sent through a grinding process. Hence the
  NIR spectra were measured on the cellulose raw material in powder
  form. Data are divided in two parts, one used for modeling and one
  part for testing.
}

\usage{
data(cellulose)
}

\format{ 
  A list with the following elements: 1) nirMN: a matrix of 180 samples x 1201 wavelengths in the VIS-NIR
  region, 2) viscoVn: a vector (length = 180) of viscosity of cellulose
  powder, and 3) classVn: a vector (length = 180) of class membership (1 or 2)
}

\value{
  For details see the Format section above.
}

\references{
  Multivariate calibration using spectral data. Simca tutorial. Umetrics.
}

\keyword{datasets}

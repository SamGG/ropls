\name{sacurine}

\alias{sacurine}

\docType{data}

\title{
  Analysis of the human adult urinary metabolome variations with age, body mass index and gender
}

\description{
  Urine samples from 183 human adults were analyzed by liquid
  chromatography coupled to high-resolution mass spectrometry (LTQ
  Orbitrap) in the negative ionization mode. A total of 109 metabolites were
  identified or annotated at the MSI level 1 or 2. After retention time
  alignment with XCMS, peaks were integrated with Quan Browser. After signal
  drift and batch effect correction of intensities, each urine profile was normalized
  to the osmolality of the sample. Finally, the data were log10 transformed.
}

\usage{
data(sacurine)
}

\format{
  A list with the following elements: 1) dataMatrix: a 183 samples x 109
  variables matrix of numeric type corresponding to the intensity
  profiles (values have been log10-transformed); 2) sampleMetadata: a 183 x 3
  data frame, with the volunteers' age ('age', numeric), body mass
  index ('bmi', numeric), and gender ('gender', factor), and 3) variableMetadata: a 109 x 3 data
  frame, with the metabolites' MSI identification level ('msiLevel':
  either 1 or 2), HMDB ID when available ('hmdb', character), chemical
  class according to the 'super class' taxonomy of HMDB
  ('chemicalClass', character)
}

\value{
  List containing the 'dataMatrix' matrix (numeric) of data (samples as
  rows, variables as columns), the 'sampleMetadata' data frame of sample
  metadata, and the variableMetadata data frame of variable metadata. Row
  names of 'dataMatrix' and 'sampleMetadata' are identical. Column names of
  'dataMatrix' are identical to row names of 'variableMetadata'. For
  details see the 'Format' section above.
}

\references{
  Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015). Analysis
  of the human adult urinary metabolome variations with age, body mass
  index and gender by implementing a comprehensive workflow for
  univariate and OPLS statistical analyses. Journal of Proteome
  Research, DOI: 10.1021/acs.jproteome.5b00354
}

\keyword{datasets}

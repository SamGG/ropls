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
  A list with the following elements: 1) profileMN: a 183 samples x 109
  variables [M]atrix of [N]umeric type corresponding to the intensity
  profiles (values have been log10-transformed); 2) sampleDF: a 183 x 3
  [D]ata [F]rame, with the volunteers' age ('ageVn', [V]ector of
  [n]umeric type), body mass index ('bmiVn'), and gender ('genderVc',
  [V]ector of [c]haracter type), and 3) variableDF: a 109 x 3 data
  frame, with the metabolites' MSI identification level ('msiLevelVn':
  either 1 or 2), HMDB ID when available ('hmdbVc'), chemical class according to the 'super class' taxonomy of HMDB ('chemicalClassVc')
}

\value{
  List containing the 'profileMN' matrix (numeric) of data (samples as
  rows, variables as columns), the 'sampleDF' data frame of sample
  metadata, and the variableDF data frame of variable metadata. Row
  names of 'profileMN' and 'sampleDF' are identical. Column names of
  'profileMN' are identical to row names of 'variableDF'. For
  details see the 'Format' section above.
}

\references{
 Thevenot E., Roux A., Ying X., Ezan E. and Junot C. Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. submitted.
}

\keyword{datasets}

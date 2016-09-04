#' Amino-Acids Dataset
#'
#' Quantitative structure property relationship (QSPR)
#'
#'
#' @name aminoacids
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item AA amino acid
#' \item PIE lipophilicity constant
#' of the AA side chain
#' \item PIF lipophilicity constant of the AA
#' side chain
#' \item DGR free energy of transfer of an AA side chain
#' from protein interior to water
#' \item SAC water-accessible
#' surface area of AA's calculated by MOLSV
#' \item MR molecular
#' refractivity
#' \item Lam polarity parameter
#' \item Vol
#' molecular volume of AA's calculated by MOLSV
#' \item DDGTS free
#' energy of unfolding of the tryptophane synthase a unit of bacteriophage T4
#' lysosome
#' }
#' @return Data frame (numeric type except the first column, which can be
#' transformed into row names) with 19 rows and the 9 columns contaning
#' information about amino acids. For details see the 'Format' section above.
#' @references Wold et al. (2001). PLS-regression: a basic tool of
#' chemometrics. Chemometrics and Intelligent Laboratory Systems. 58:109-130.
#' @source 'aminoacids' dataset.
#' @keywords datasets
NULL

#' NIR-Viscosity example data set to illustrate multivariate calibration using
#' PLS, spectral filtering and OPLS
#'
#' The data were collected at Akzo Nobel, Ornkoldsvik (Sweden). The raw
#' material for their cellulose derivative process is delivered to the factory
#' in form of cellulose sheets. Before entering the process the cellulose
#' sheets are controlled by a viscosity measurement, which functions as a
#' steering parameter for that particular batch.  In this data set NIR spectra
#' for 180 cellulose sheets were collected after the sheets had been sent
#' through a grinding process. Hence the NIR spectra were measured on the
#' cellulose raw material in powder form. Data are divided in two parts, one
#' used for modeling and one part for testing.
#'
#'
#' @name cellulose
#' @docType data
#' @format A list with the following elements:
#' \itemize{
#' \item nirMN a matrix of 180 samples x 1201 wavelengths in the VIS-NIR region
#' \item viscoVn a vector (length = 180) of viscosity of cellulose powder
#' \item classVn a vector (length = 180) of class membership (1 or 2)
#' }
#' @return For details see the Format section above.
#' @references Multivariate calibration using spectral data. Simca tutorial.
#' Umetrics.
#' @keywords datasets
NULL

#' Octane of various blends of gasoline
#'
#' Twelve mixture component proportions of the blend are analysed
#'
#'
#' @name cornell
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item num mixture number
#' \item x1 proportion of
#' component 1
#' \item x2 proportion of component 2
#' \item x3 proportion of component 3
#' \item x4 proportion of component 4
#' \item x5 proportion of component 5
#' \item x6 proportion of component 6
#' \item x7 proportion of component 7 Note: the 7 variables are correlated since they sum up to 1
#' \item y octane (quantitative variable)
#' }
#' @return Data frame (numeric type only; the first column can be transformed
#' into row names) with 12 rows and 9 columns corresponding to the 'num'ber of
#' the mixture (column 1), the proportion of each of the 7 'x' components
#' within the mixture (columns 2-8), and the octane indice 'y' (column 9). For
#' details see the 'Format' section above.
#' @references Tenenhaus (1998). La regression PLS: theorie et pratique. Paris:
#' Editions Technip.
#' @source Tenenhaus (1998), Table 6, page 78.
#' @keywords datasets
NULL

#' Food consumption patterns accross European countries (FOODS)
#'
#' The relative consumption of 20 food items was compiled for 16 countries. The
#' values range between 0 and 100 percent and a high value corresponds to a
#' high consumption. The dataset contains 3 missing data.
#'
#'
#' @name foods
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item Country Name of the country
#' \item Gr_CoffeGround Coffee
#' \item Inst_Coffe Instant Coffee
#' \item Tea Tea \item Sweetner Sweetner
#' \item Biscuits Biscuits \item Pa_Soup Powder Soup
#' \item Ti_Soup Tin Soup \item In_Potat Instant Potatoes
#' \item Fro_Fish Frozen Fish
#' \item Fro_Veg Frozen Vegetables
#' \item Apples Apples
#' \item Oranges Oranges
#' \item Ti_Fruit Tin Fruit
#' \item Jam Jam
#' \item Garlic Garlic
#' \item Butter Butter
#' \item Margarine Margarine
#' \item Olive_Oil Olive Oil
#' \item Yoghurt Yoghurt
#' \item Crisp_Brea Crisp Bread
#' }
#' @return Data frame (numeric type except the first column, which can be
#' transformed into row names) with 16 rows and 21 columns, corresponding to
#' the 'Country' (column 1), followed by the consumption of each of the 20 food
#' items (columns 2-21). For details see the 'Format' section above.
#' @references Eriksson et al. (2006). Multi- and Megarvariate Data Analysis.
#' Umetrics Academy. pp.10, 33, 48.
#' @keywords datasets
NULL

#' Linnerud Dataset
#'
#' Three physiological and three exercise variables are measured on twenty
#' middle-aged men in a fitness club.
#'
#'
#' @name linnerud
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item num subject number
#' \item weight weight
#' \item waist waist
#' \item pulse pulse
#' \item pullUp pull-up
#' \item squat situp
#' \item jump jump
#' }
#' @return Data frame (numeric type only; the first column can be transformed
#' into row names) with 20 rows and 7 columns corresponding to the subject's
#' 'num'ber (column 1), the 3 physiological variables (columns 2-4), and the 3
#' exercise variables (columns 5-7). For details see the 'Format' section
#' above.
#' @references Tenenhaus (1998). La regression PLS: theorie et pratique. Paris:
#' Editions Technip.
#' @source 'mixOmics' 'linnerud' dataset.
#' @keywords datasets
NULL

#' A multi response optimization data set (LOWARP)
#'
#' This example concerns the development of a polymer similar to that used in
#' the plastic covering of mobile phones. The desired profile of the polymer
#' was low warp and high strength. Four constituents (glas, crtp, mica, and
#' amtp) were varied in the polymer formulation by means of a 17 run mixture
#' design. For each new polymer, i.e., each new experiment in the mixture
#' design, 14 responses relating to both warp and strength were measured on the
#' product. The objective of the data analysis was to uncover which combination
#' of factors (the four ingredients) gave polymers with low warp and high
#' strength. The data set contains 10 missing values (NA).
#'
#'
#' @name lowarp
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item num mixture number
#' \item glas glas constituent
#' \item crtp crtp constituent
#' \item mica mica constituent
#' \item amtp amtp constituent
#' \item wrp1 warp response 1
#' \item wrp2 warp response 2
#' \item wrp3 warp response 3
#' \item wrp4 warp response 4
#' \item wrp5 warp response 5
#' \item wrp6 warp response 6
#' \item wrp7 warp response 7
#' \item wrp8 warp response 8
#' \item st1 strength response 1
#' \item st2 strength response 2
#' \item st3 strength response 3
#' \item st4 strength response 4
#' \item st5 strength response 5
#' \item st6 strength response 6
#' }
#' @return Data frame (numeric type only; the first column can be transformed
#' into row names) with 17 rows and 19 columns corresponding to the subject's
#' 'num'ber (column 1), the 4 constituent variables (columns 2-5), the 8 warp
#' responses (columns 6-13), and the 6 strength responses (columns 14-19). For
#' details see the 'Format' section above.
#' @references Eriksson et al. (2006). Multi- and Megarvariate Data Analysis.
#' Umetrics Academy. pp.16, 77, 209.
#' @keywords datasets
NULL

#' 'mark' Dataset
#'
#' Examination marks obtained by French students in Mathematics, Physics,
#' French and English
#'
#'
#' @name mark
#' @docType data
#' @format A data frame with the following parameters:
#' \itemize{
#' \item nom names of the students
#' \item math marks in mathematics
#' \item phys marks in physics
#' \item fran marks in french
#' \item angl marks in english
#' }
#' @return Data frame (numeric type except the first column, which can be
#' transformed into row names) with 9 rows and 5 columns, corresponding to the
#' name of the students (column 1), followed by the marks obtained in Maths,
#' Physics, French and English (columns 2-5). For details see the 'Format'
#' section above.
#' @references Baccini (2010). Statistique Descriptive Multidimensionnelle
#' (pour les nuls).
#' @source 'mark' dataset.
#' @keywords datasets
NULL

#' Analysis of the human adult urinary metabolome variations with age, body
#' mass index and gender
#'
#' Urine samples from 183 human adults were analyzed by liquid chromatography
#' coupled to high-resolution mass spectrometry (LTQ Orbitrap) in the negative
#' ionization mode. A total of 109 metabolites were identified or annotated at
#' the MSI level 1 or 2. After retention time alignment with XCMS, peaks were
#' integrated with Quan Browser. After signal drift and batch effect correction
#' of intensities, each urine profile was normalized to the osmolality of the
#' sample. Finally, the data were log10 transformed.
#'
#' @name sacurine
#' @docType data
#' @format A list with the following elements:
#' \itemize{
#' \item dataMatrix a 183 samples x
#' 109 variables matrix of numeric type corresponding to the intensity profiles
#' (values have been log10-transformed)
#' \item sampleMetadata a 183 x 3 data frame, with the volunteers' age
#' ('age', numeric), body mass index ('bmi',
#' numeric), and gender ('gender', factor)
#' \item variableMetadata a 109 x 3 data frame, with the metabolites'
#' MSI identification level ('msiLevel':
#' either 1 or 2), HMDB ID when available ('hmdb', character), chemical class
#' according to the 'super class' taxonomy of HMDB ('chemicalClass', character)
#' }
#' @return List containing the 'dataMatrix' matrix (numeric) of data (samples
#' as rows, variables as columns), the 'sampleMetadata' data frame of sample
#' metadata, and the variableMetadata data frame of variable metadata. Row
#' names of 'dataMatrix' and 'sampleMetadata' are identical. Column names of
#' 'dataMatrix' are identical to row names of 'variableMetadata'. For details
#' see the 'Format' section above.
#' @references Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015).
#' Analysis of the human adult urinary metabolome variations with age, body
#' mass index and gender by implementing a comprehensive workflow for
#' univariate and OPLS statistical analyses. Journal of Proteome Research, DOI:
#' 10.1021/acs.jproteome.5b00354
#' @keywords datasets
NULL

\name{ropls-package}

\alias{ropls-package}

\alias{ropls}

\docType{package}

\title{
  PCA, PLS(-DA) and OPLS(-DA) for multivariate analysis and feature selection of omics data
}

\description{
Latent variable modeling with Principal Component Analysis (PCA) and
Partial Least Squares (PLS) are powerful methods for visualization,
regression, classification, and feature selection of omics data where
the number of variables exceeds the number of samples and with
multicollinearity among variables. Orthogonal Partial Least Squares
(OPLS) enables to separately model the variation correlated (predictive)
to the factor of interest and the uncorrelated (orthogonal)
variation. While performing similarly to PLS, OPLS facilitates
interpretation. Successful applications of these chemometrics techniques
include spectroscopic data such as Raman spectroscopy, nuclear magnetic
resonance (NMR), mass spectrometry (MS) in metabolomics and proteomics,
but also transcriptomics data. In addition to scores, loadings and
weights plots, the package provides metrics and graphics to determine
the optimal number of components (e.g. with the R2 and Q2 coefficients),
check the validity of the model by permutation testing, detect outliers,
and perform feature selection (e.g. with Variable Importance in
Projection or regression coefficients). The package can be accessed via
a user interface on the Workflow4Metabolomics.org online resource for
computational metabolomics (built upon the Galaxy environment).
}

\details{
  \tabular{ll}{
    Package: \tab ropls\cr
    Type: \tab Package\cr
    Version: \tab 1.3.9\cr
    Date: \tab 2016-01-10\cr
    License: \tab CeCILL\cr
    LazyLoad: \tab yes\cr
  }
}

\author{
  Etienne A. Thevenot.

  Maintainer: Etienne A. Thevenot <etienne.thevenot@cea.fr>

  Please cite: Etienne A. Thevenot, Aurelie Roux, Ying Xu, Eric Ezan,
  and Christophe Junot (2015). Analysis of the human adult urinary metabolome
  variations with age, body mass index and gender by implementing a
  comprehensive workflow for univariate and OPLS statistical
  analyses. Journal of Proteome Research. DOI: 10.1021/acs.jproteome.5b00354.

}

\keyword{ package }

\examples{

#### PCA

data(foods) ## see Eriksson et al. (2001); presence of 3 missing values (NA)
head(foods)
foodMN <- as.matrix(foods[, colnames(foods) != "Country"])
rownames(foodMN) <- foods[, "Country"]
head(foodMN)
foo.pca <- opls(foodMN)

#### PLS with a single response

data(cornell) ## see Tenenhaus, 1998
head(cornell)
cornell.pls <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                    cornell[, "y"])
## Complementary graphics
plot(cornell.pls, typeVc = c("outlier", "predict-train", "xy-score", "xy-weight"))

#### PLS with multiple (quantitative) responses

data(lowarp) ## see Eriksson et al. (2001); presence of NAs
head(lowarp)
lowarp.pls <- opls(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
                   as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
                                      grepl("^st", colnames(lowarp))]))

#### PLS-DA

data(sacurine)
attach(sacurine)
sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])

#### OPLS-DA

sacurine.oplsda <- opls(dataMatrix, sampleMetadata[, "gender"], predI = 1, orthoI = NA)

detach(sacurine)

}

\name{ropls-package}

\alias{ropls-package}

\alias{ropls}

\docType{package}

\title{
  PCA, PLS(-DA) and OPLS(-DA) with R
}

\description{
  Function and datasets for PCA, PLS, and OPLS regression, classification, and cross-validation with the NIPALS algorithm
}

\details{
  \tabular{ll}{
    Package: \tab ropls\cr
    Type: \tab Package\cr
    Version: \tab 0.99.4\cr
    Date: \tab 2015-07-08\cr
    License: \tab CeCILL\cr
    LazyLoad: \tab yes\cr
  }
}

\author{
  Etienne A. Thevenot.

  Maintainer: Etienne A. Thevenot <etienne.thevenot@cea.fr>

  Please cite: Etienne A. Thevenot, Aurelie Roux, Ying Xu, Eric Ezan,
  and Christophe Junot. Analysis of the human adult urinary metabolome
  variations with age, body mass index and gender by implementing a
  comprehensive workflow for univariate and OPLS statistical analyses. submitted.

}

\keyword{ package }

\examples{

#### PCA

data(foods) ## see Eriksson et al. (2001); presence of 3 missing values (NA)
head(foods)
fooMN <- as.matrix(foods[, colnames(foods) != "Country"])
rownames(fooMN) <- foods[, "Country"]
head(fooMN)
fooPcaLs <- roplsF(fooMN)


#### PLS with a single response

data(cornell) ## see Tenenhaus, 1998
head(cornell)
corPlsLs <- roplsF(as.matrix(cornell[, grep("x", colnames(cornell))]),
                   matrix(cornell[, "y"], ncol = 1))
## Complementary graphics
corPlsLs <- roplsF(as.matrix(cornell[, grep("x", colnames(cornell))]),
                   matrix(cornell[, "y"], ncol = 1),
                   plotVc = c("outlier", "predict-train", "xy-score", "xy-weight"))

#### PLS with multiple (quantitative) responses

data(lowarp) ## see Eriksson et al. (2001); presence of NAs
head(lowarp)
lowPlsLs <- roplsF(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
                   as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
                                      grepl("^st", colnames(lowarp))]))


#### PLS-DA

data(sacurine)
attach(sacurine)
sacPlsLs <- roplsF(profileMN, matrix(sampleDF[, "genderVc"], ncol = 1))

#### OPLS-DA

sacPlsLs <- roplsF(profileMN, matrix(sampleDF[, "genderVc"], ncol = 1), predN
= 1, orthoN = NA)

detach(sacurine)

}

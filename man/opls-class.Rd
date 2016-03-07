\name{opls-class}
\docType{class}

\alias{opls-class}

\title{Class "opls"}

\description{The opls object class}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("opls", ...)}
  or by calling the \code{opls} function
}

\section{Slots}{
  \describe{
    \item{\code{typeC}:}{
      character: model type (PCA, PLS, PLS-DA, OPLS, or OPLS-DA)
    }
    \item{\code{descriptionMC}:}{
      character matrix: Description of the data
      set (number of samples, variables, etc.)
    }
    \item{\code{modelDF}:}{
      Data frame with the model overview (number of components, R2X,
      R2X(cum), R2Y, R2Y(cum), Q2, Q2(cum), significance, iterations)
    }
    \item{\code{summaryDF}:}{
      Data frame with the model summary (cumulated R2X, R2Y and Q2); RMSEE
      is the square root of the mean error between the actual and the
      predicted responses
    }
    \item{\code{subsetVi}:}{
      Integer vector: Indices of observations in the training data set
    }
    \item{\code{pcaVarVn}:}{
      PCA: Numerical vector of variances of length: predI
    }
    \item{\code{vipVn}:}{
      PLS(-DA): Numerical vector of Variable Importance in Projection;
      OPLS(-DA): Numerical vector of Variable Importance for Prediction
      (VIP4,p from Galindo-Prieto et al, 2014)
    }
    \item{\code{orthoVipVn}:}{
      OPLS(-DA): Numerical vector of Variable Importance for Orthogonal
      Modeling (VIP4,o from Galindo-Prieto et al, 2014)
    }
    \item{\code{coefficientMN}:}{
      (O)PLS(-DA): Numerical matrix of regression coefficients (B; dimensions:
      ncol(x) x number of responses; B = W*C' and Y = XB + F
    }
    \item{\code{xMeanVn}:}{
      Numerical vector: variable means of the 'x' matrix
    }
    \item{\code{xSdVn}:}{
      Numerical vector: variable standard deviations of the 'x' matrix
    }
    \item{\code{yMeanVn}:}{
      (O)PLS: Numerical vector: variable means of the 'y' response (transformed
      into a dummy matrix in case it is of 'character' mode initially)
    }
    \item{\code{ySdVn}:}{
      (O)PLS: Numerical vector: variable standard deviations of the 'y' response (transformed
      into a dummy matrix in case it is of 'character' mode initially)
    }
    \item{\code{xZeroVarVi}:}{
      Numerical vector: indices of variables with variance < 2.22e-16
      which were excluded from 'x' before building the model
    }
    \item{\code{scoreMN}:}{
      Numerical matrix of x scores (T; dimensions: nrow(x) x predI)
      X = TP' + E; Y = TC' + F
    }
    \item{\code{loadingMN}:}{
      Numerical matrix of x loadings (P; dimensions: ncol(x) x predI)
      X = TP' + E
    }
    \item{\code{weightMN}:}{
      (O)PLS: Numerical matrix of x weights (W; same dimensions as loadingMN)
    }
    \item{\code{orthoScoreMN}:}{
      OPLS: Numerical matrix of orthogonal scores (Tortho; dimensions:
      nrow(x) x number of orthogonal components)
    }
    \item{\code{orthoLoadingMN}:}{
      OPLS: Numerical matrix of orthogonal loadings (Portho; dimensions:
      ncol(x) x number of orthogonal components)
    }
    \item{\code{orthoWeightMN}:}{
      OPLS: Numerical matrix of orthogonal weights (same dimensions as orthoLoadingMN)
    }
    \item{\code{cMN}:}{
      (O)PLS: Numerical matrix of Y weights (C); dimensions: number of
      responses or number of classes in case of qualitative response
      with more than 2 classes x number of predictive components; Y = TC' + F
    }
    \item{\code{coMN}:}{
      (O)PLS: Numerical matrix of Y orthogonal weights; dimensions: number of
      responses or number of classes in case of qualitative response
      with more than 2 classes x number of orthogonal components
    }
    \item{\code{uMN}:}{
      (O)PLS: Numerical matrix of Y scores (U; same dimensions as scoreMN); Y = UC' + G
    }
    \item{\code{weightStarMN}:}{
      Numerical matrix of projections (W*; same dimensions as loadingMN); whereas
      columns of weightMN are derived from successively deflated matrices,
      columns of weightStarMN relate to the original 'x' matrix: T = XW*; W*=W(P'W)inv
    }
    \item{\code{suppLs}:}{
      List of additional objects to be used internally by the 'print',
      'plot', and 'predict' methods
    }
  }
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(object = "opls")}: print method for
      (O)PLS(-DA) models}
    \item{show}{\code{signature(object = "opls")}: show method for
      (O)PLS(-DA) models}
    \item{plot}{\code{signature(object = "opls")}: plot method for
      (O)PLS(-DA) models}
    \item{predict}{\code{signature(object = "opls")}: plot method for
      (O)PLS(-DA) models}
    \item{coefficients}{\code{signature(object = "opls")}: }
    \item{fitted}{\code{signature(object = "opls")}: returns predictions
      of the (O)PLS(-DA) model}
    \item{residuals}{\code{signature(object = "opls")}: returns the
      residuals of the (O)PLS(-DA) regression model}
    \item{tested}{\code{signature(object = "opls")}: returns the
      predictions of the (O)PLS(-DA) trained on the 'subset' observations only}
    \item{getSummaryDF}{\code{signature(object = "opls")}: returns
      summary information about the model}
    \item{getPcaVarVn}{\code{signature(object = "opls")}: returns the variance of the components (score vectors)}
    \item{getScoreMN}{\code{signature(object = "opls")}: returns the
      (orthogonal) score matrix of the PCA/(O)PLS(-DA) model}
    \item{getLoadingMN}{\code{signature(object = "opls")}: returns the
      (orthogonal) loading matrix of the PCA/(O)PLS(-DA) model}
    \item{getWeightMN}{\code{signature(object = "opls")}: returns the
      (orthogonal) weight matrix of the (O)PLS(-DA) model}
    \item{getVipVn}{\code{signature(object = "opls")}: returns the
      (orthogonal) VIP values of the (O)PLS(-DA) model}
    \item{getSubsetVi}{\code{signature(object = "opls")}: extracts the indices of the samples used for building the model (when
  a subset argument has been specified)}
  }
}

\author{
Etienne Thevenot (CEA)
}

\seealso{
  \code{\link{opls}},
}

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


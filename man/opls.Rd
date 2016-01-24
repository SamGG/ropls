\name{opls}

\alias{opls}
\alias{opls.default}

\title{
    PCA, PLS(-DA), and OPLS(-DA)
}

\description{
   PCA, PLS, and OPLS regression, classification, and cross-validation with the NIPALS algorithm
}

\usage{
\method{opls}{default}(x,
y = NULL,
predI = NA,
orthoI = 0,
algoC = c("default", "nipals", "svd")[1],
crossvalI = 7,
log10L = FALSE,
permI = 20,
scaleC = c("center", "pareto", "standard")[3],
subset = NULL,
printL = TRUE,
plotL = TRUE,

.sinkC = NULL,
...)
}

\arguments{
  \item{x}{
    Numerical data frame or matrix (observations x variables); NAs are allowed
  }
  \item{y}{
    Response to be modelled: Either 1) 'NULL' for PCA (default) or
    2) a numerical vector (same length as 'x' row number) for single response
    (O)PLS, or 3) a numerical matrix (same row number as 'x')
    for multiple response PLS or 4) a factor (same length as 'x' row
    number) for (O)PLS-DA. Note that, for convenience, character vectors are also accepted
    for (O)PLS-DA as well as single column numerical (resp. character) matrices for (O)PLS (respectively (O)PLS-DA). NAs are allowed in
    numeric responses.
  }
  \item{predI}{
    Integer: number of components (predictive componenents in case of
    PLS and OPLS) to extract; for OPLS, predI is (automatically) set to
    1; if set to NA [default], autofit is performed: a maximum of 10 components are
    extracted until (i) PCA case: the variance is less than the mean
    variance of all components (note that this rule requires all
    components to be computed and can be quite time-consuming for large
    datasets) or (ii) PLS case: either R2Y of the component is < 0.01 (N4
    rule) or Q2Y is < 0 (for more than 100 observations) or 0.05 otherwise (R1 rule)
  }
  \item{orthoI}{
    Integer: number of orthogonal components (for OPLS only); when set
    to 0 [default], PLS will be performed; otherwise OPLS will be
    peformed; when set to NA, OPLS is performed and the number of
    orthogonal components is automatically computed by using the
    cross-validation (with a maximum of 9 orthogonal components).
  }
  \item{algoC}{
    Default algorithm is 'svd' for PCA (in case of no missing values in
    'x'; 'nipals' otherwise) and 'nipals' for PLS and OPLS; when asking to use 'svd' for PCA on an 'x' matrix containing missing values, NAs are set to half the minimum of non-missing values and a warning is generated
  }
  \item{crossvalI}{
    Integer: number of cross-validation segments (default is 7); The
    number of samples (rows of 'x') must be at least >= crossvalI
  }
  \item{log10L}{
    Should the 'x' matrix be log10 transformed? Zeros are set to 1 prior
    to transformation
  }
  \item{permI}{
    Integer: number of random permutations of response labels to
    estimate R2Y and Q2Y significance by permutation testing [default is
    20 for single response models (without train/test partition), and 0 otherwise]
  }
  \item{scaleC}{
    Character: either mean-centering ('center'),
    mean-centering and pareto scaling ('pareto'), or mean-centering and unit
    variance scaling ('standard') [default]
  }
  \item{subset}{
    Integer vector: indices of the observations to be used for training (in a
    classification scheme); use NULL [default] for no
    partition of the dataset; use 'odd' for a partition of the dataset
    in two equal sizes (with respect to the classes proportions)
  }
  \item{printL}{
    Logical: Should informations regarding the data set and the model be
    printed? [default = TRUE]
  }
  \item{plotL}{
    Logical: Should the 'summary' plot be displayed? [default = TRUE]
  }
  \item{.sinkC}{
    Character: Name of the file for R output diversion [default = NULL:
    no diversion]; Diversion of messages is required for the integration
    into Galaxy
  }
  \item{\dots}{
    Currently not used.
  }
}


\value{
  \item{typeC}{
    Character: model type (PCA, PLS, PLS-DA, OPLS, or OPLS-DA)
  }
  \item{descriptionMC}{
    Character matrix: Description of the data set (number of samples,
    variables, etc.)
  }
  \item{modelDF}{
    Data frame with the model overview (number of components, R2X,
    R2X(cum), R2Y, R2Y(cum), Q2, Q2(cum), significance, iterations)
  }
  \item{summaryDF}{
    Data frame with the model summary (cumulated R2X, R2Y and Q2); RMSEE
    is the square root of the mean error between the actual and the
    predicted responses
  }
  \item{subset}{
    Integer vector: Indices of observations in the training data set
  }
  \item{pcaVarVn}{
    PCA: Numerical vector of variances of length: predI
  }
  \item{vipVn}{
    PLS(-DA): Numerical vector of Variable Importance in Projection;
    OPLS(-DA): Numerical vector of Variable Importance for Prediction
    (VIP4,p from Galindo-Prieto et al, 2014)
  }
  \item{orthoVipVn}{
    OPLS(-DA): Numerical vector of Variable Importance for Orthogonal
    Modeling (VIP4,o from Galindo-Prieto et al, 2014)
  }
  \item{fitted}{
    (O)PLS: Predicted reponse based on the reference
    subset (which is the whole 'x' matrix in case of no partition: subset = NULL)
  }
  \item{tested}{
    (O)PLS: In case of training/test partition indicated via the subset
    argument, predicted response based on the test subset
  }
  \item{coefficients}{
    (O)PLS(-DA): Numerical matrix of regression coefficients (B; dimensions:
  ncol(x) x number of responses; B = W*C' and Y = XB + F
  }
  \item{residuals}{
    (O)PLS: Numerical vector equal to the difference between the true
  and the predicted responses
  }
  \item{xMeanVn}{
     Numerical vector: variable means of the 'x' matrix
  }
  \item{xSdVn}{
    Numerical vector: variable standard deviations of the 'x' matrix
  }
  \item{yMeanVn}{
    (O)PLS: Numerical vector: variable means of the 'y' response (transformed
    into a dummy matrix in case it is of 'character' mode initially)
  }
  \item{ySdVn}{
    (O)PLS: Numerical vector: variable standard deviations of the 'y' response (transformed
    into a dummy matrix in case it is of 'character' mode initially)
  }
  \item{xZeroVarVi}{
    Numerical vector: indices of variables with variance < 2.22e-16
  which were excluded from 'x' before building the model
  }
  \item{scoreMN}{
    Numerical matrix of x scores (T; dimensions: nrow(x) x predI)
    X = TP' + E; Y = TC' + F
  }
  \item{loadingMN}{
    Numerical matrix of x loadings (P; dimensions: ncol(x) x predI)
    X = TP' + E
  }
  \item{weightMN}{
    (O)PLS: Numerical matrix of x weights (W; same dimensions as loadingMN)
  }
    \item{orthoScoreMN}{
    OPLS: Numerical matrix of orthogonal scores (Tortho; dimensions:
  nrow(x) x number of orthogonal components)
  }
  \item{orthoLoadingMN}{
    OPLS: Numerical matrix of orthogonal loadings (Portho; dimensions:
  ncol(x) x number of orthogonal components)
  }
  \item{orthoWeightMN}{
    OPLS: Numerical matrix of orthogonal weights (same dimensions as orthoLoadingMN)
  }
  \item{cMN}{
    (O)PLS: Numerical matrix of Y weights (C; dimensions: number of
  responses or number of classes in case of qualitative response) x
  number of predictive components; Y = TC' + F
  }
  \item{uMN}{
    (O)PLS: Numerical matrix of Y scores (U; same dimensions as scoreMN); Y = UC' + G
  }
  \item{weightStarMN}{
    Numerical matrix of projections (W*; same dimensions as loadingMN); whereas
    columns of weightMN are derived from successively deflated matrices,
    columns of weightStarMN relate to the original 'x' matrix: T = XW*; W*=W(P'W)inv
  }
  \item{suppLs}{
    List of additional objects to be used internally by the 'print',
  'plot', and 'predict' methods
    }
}

\references{
  Eriksson et al. (2006). Multi- and Megarvariate Data
  Analysis. Umetrics Academy.
  Rosipal and Kramer (2006). Overview and recent advances in partial least squares
  Tenenhaus (1990). La regression PLS : theorie et pratique. Technip.
  Wehrens (2011). Chemometrics with R. Springer.
  Wold et al. (2001). PLS-regression: a basic tool of chemometrics
}

\author{
  Etienne Thevenot (CEA)
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


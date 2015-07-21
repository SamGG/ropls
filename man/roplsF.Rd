\name{roplsF}

\alias{roplsF}

\title{
    PCA, PLS(-DA), and OPLS(-DA)
}

\description{
   PCA, PLS, and OPLS regression, classification, and cross-validation with the NIPALS algorithm
}

\usage{
roplsF(xMN,
yMCN = NULL,
predI = NA,
orthoI = 0,
plotVc = c("correlation", "none", "outlier", "overview",
"permutation", "predict-train", "predict-test", "summary", "x-loading",
"x-score", "x-variance", "xy-score", "xy-weight")[8],
algoC = c("default", "nipals", "svd")[1],
crossvalI = 7,
fileFig.pdfC = NULL,
fileInfo.txtC = NULL,
log10L = FALSE,
permI = 0,
scaleC = c("center", "pareto", "standard")[3],
testVi = NULL,
verboseC = c("all", "none", "overview", "summary")[1],
parAsColVcn = NA,
parCexN = 1,
parCompVi = c(1, 2),
parDevNewL = TRUE,
parEllipsesL = NA,
parLabVc = NA,
parTitleL = TRUE,
parTopLoadI = 3)
}

\arguments{
  \item{xMN}{
    Numerical matrix (observations * X variables); NAs are allowed
  }
  \item{yMCN}{
    Numerical matrix or 1-column character matrix
    (observations * Y variables); when set to NULL [default], PCA will be performed;
    otherwise PLS is performed; when yMCN is a 1 column character matrix
    (for classification), it will be transformed into a dummy matrix
    within the roplsF function; NAs are allowed (except for the 'classification' case)
  }
  \item{predI}{
    Integer: number of components (predictive componenents in case of
    PLS and OPLS) to extract; for OPLS, predI is (automatically) set to
    1; if set to NA [default], autofit is performed: components are
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
    cross-validation
  }
  \item{plotVc}{
    Character vector: the following plots are available:
    'correlation': Variable correlations with the components, 'none': No
    plotting, 'outlier':
    Observation diagnostics (score and orthogonal distances), 'overview':
    Model overview showing R2Ycum and Q2cum (or 'Variance explained' for
    PCA), 'permutation': Scatterplot of R2Y and Q2Y actual and simulated
    models after random permutation of response values;
    'predict-train' and 'predict-test': Predicted vs Actual Y for reference
    and test sets (only if Y has a single column), 'summary' [default]: 4-plot
    summary showing overview, permutation, x-loading and x-score together,
    'x-variance': Spread of raw
    variables corresp. with min, median, and max variances,
    'x-loading': X-loadings, 'x-score': X-Scores, 'xy-score': XY-Scores,
    'xy-weight': XY-Weights
  }
  \item{algoC}{
    Default algorithm is 'svd' for PCA (in case of no missing values in
    xMN; 'nipals' otherwise) and 'nipals' for PLS and OPLS; when asking to use 'svd' for PCA on an xMN matrix containing missing values, NAs are set to half the minimum of non-missing values and a warning is generated
  }
  \item{crossvalI}{
    Integer: number of cross-validation segments (default is 7); The
    number of samples (rows of xMN) must be at least >= crossvalI
  }
  \item{fileFig.pdfC}{
    Figure filename (e.g. in case of batch mode) ending with '.pdf'; for multiple graphics,
    set parLayL to TRUE; default is NULL (no saving; displaying instead)
  }
  \item{fileInfo.txtC}{
    Log filename (e.g. in case of batch mode); default is NULL (no log
    file; interactive printing instead)
  }
  \item{log10L}{
    Should the xMN matrix be log10 transformed? Zeros are set to 1 prior
    to transformation
  }
  \item{permI}{
    Integer: number of random permutations of response labels to
    estimate R2Y and Q2Y confidence [default: 0]
  }
  \item{scaleC}{
    Character: either mean-centering ('center'),
    mean-centering and pareto scaling ('pareto'), or mean-centering and unit
    variance scaling ('standard') [default]
  }
  \item{testVi}{
    Integer vector: indices of the observations to be tested (in a
    classification scheme); use NULL [default] for no
    partition of the dataset; use 'odd' for a partition of the dataset
    in two equal sizes (with respect to the classes proportions)
  }
  \item{verboseC}{
    Character: Messages and warnings to be printed: 'all' [default],
    'none', model 'overview' or 'summary'
  }
  \item{parAsColVcn}{
    Optional character or numeric vector to be converted into colors (with the
    .colorF function) for the score plot; default is NA [ie colors will be
    converted from yMCN in case of (O)PLS(-DA) or will be 'black' for PCA]
  }
  \item{parCexN}{
    Numeric: amount by which plotting text should be magnified relative to the default
  }
  \item{parCompVi}{
    Integer vector of length 2: indices of the two components to be displayed
    on the score plot (first two components by default)
  }
  \item{parDevNewL}{
    Should the graphics be displayed in a new window [default]; If
    FALSE, parLayL must be set to FALSE also
  }
  \item{parEllipsesL}{
    Should the Mahalanobis ellipses be drawn? If 'NA' [default],
    ellipses are drawn when either a parAsColVcn is provided (PCA
    case), or when yMCN is a single column matrix of character type
    ((O)PLS-DA cases).
  }
  \item{parLabVc}{
    Optional character vector for the labels of observations on the
    plot; default is NA [ie rownames of xMN, if available, or indices of
    xMN, otherwise, will be used]
  }
  \item{parTitleL}{
    Should the titles of the plots be printed on the graphics (default =
    TRUE); It may be convenient to set this argument to FALSE when
    the user wishes to add specific titles a posteriori
  }
  \item{parTopLoadI}{
    Integer: number of variables most contributing to loadings (colored
    in red on the plot) for which loading values and correlation
    with the components will be printed to facilitate interpretation
  }
}


\value{
  \item{typC}{
    Character: model type (PCA, PLS, PLS-DA, OPLS, or OPLS-DA)
  }
  \item{predI}{
    Integer: number of (predictive) components
  }
  \item{xMeanVn}{
     Numerical vector: variable means of the xMN matrix
  }
  \item{xSdVn}{
    Numerical vector: variable standard deviations of the xMN matrix
  }
  \item{testVi}{
    Integer vector: indices of observations in the test data set
  }
  \item{tMN}{
    Numerical matrix of xMN scores (dimensions: nrow(xMN) x predI)
    X = TP' + E; Y = TC' + F
  }
  \item{pMN}{
    Numerical matrix of xMN loadings (dimensions: ncol(xMN) x predI)
    X = TP' + E
  }
  \item{varVn}{
    PCA: Numerical vector of variances of length: predI
  }
  \item{modelDF}{
    Data frame with the model overview (number of components, R2X,
    R2X(cum), R2Y, R2Y(cum), Q2, Q2(cum), significance, iterations)
  }
  \item{summaryDF}{
    Data frame with the model summary (cumulated R2X, R2Y and Q2); RMSEE
    is the square root of the mean error between the actual and the
    predicted responses (model with all components)
  }
  \item{yMeanVn}{
    (O)PLS: Numerical vector: variable means of the yMCN matrix (transformed
    into a dummy matrix in case it is of 'character' mode initially)
  }
  \item{ySdVn}{
    (O)PLS: Numerical vector: variable standard deviations of the yMCN matrix (transformed
    into a dummy matrix in case it is of 'character' mode initially)
  }
  \item{wMN}{
    (O)PLS: Numerical matrix of xMN-yMCN weights (same dimensions as pMN)
  }
  \item{cMN}{
    (O)PLS: Numerical matrix of yMCN loadings; dimensions: ncol(yMCN or
    corresponding dummy matrix) x predI; Y = TC' + F
  }
  \item{uMN}{
    (O)PLS: Numerical matrix of yMCN scores (same dimensions as tMN); Y = UC' + G
  }
  \item{rMN}{
    Numerical matrix of projections W* (same dimensions as pMN); whereas
    columns of wMN are derived from successively deflated matrices,
    columns of rMN relate to the original X matrix: T = XW*
  }
  \item{bMN}{
    (O)PLS: Numerical matrix of regression coefficients (dimensions: ncol(xMN) x
    ncol(yMCN)); B = W*C' and Y = XB + F
  }
  \item{vipVn}{
    (O)PLS: Numerical vector of Variable Importance in Projection
    (computed with all components)
  }
  \item{orthoI}{
    OPLS: Integer: number of orthogonal components
  }
  \item{tOrthoMN}{
    OPLS: Numerical matrix of orthogonal scores (dimensions: nrow(xMN) x orthoI)
  }
  \item{pOrthoMN}{
    OPLS: Numerical matrix of orthogonal loadings (dimensions: ncol(xMN) x orthoI)
  }
  \item{wOrthoMN}{
    OPLS: Numerical matrix of orthogonal weights (same dimensions as pOrthoMN)
  }
  \item{.char2numF}{
    (O)PLS-DA: Function for character to numeric convertion of the
    response
  }
  \item{xZeroVarVi}{
    Numerical vector: indices of variables with variance < 2.22e-16
  }
  \item{permMN}{
    (O)PLS: Numerical matrix of the R2Y and Q2Y values from the actual
    model and the models after random permutation of response values (if
    permI > 0)
  }
  \item{tCompMN}{
    Score plot matrix
  }
  \item{pCompMN}{
    Loading plot matrix
  }
  \item{topLoadMN}{
    Matrix (numeric; dimensions = topLoadN x 4): For the 'topLoadN'
    variables most contributing to loadings, gives the loadings values
    and the correlation with each of the 2 'parCompN' displayed component.
  }
  \item{yPredMCN}{
    (O)PLS: Predicted reponse based on the reference
    subset (which is the whole xMN matrix in case of no partition: testVi = NA)
  }
  \item{yTestMCN}{
    (O)PLS: In case of training/test partition indicated via the testVi
  argument, predicted response based on the test subset
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
sacPlsLs <- roplsF(dataMatrix, matrix(sampleMetadata[, "genderVc"], ncol = 1))

#### OPLS-DA

sacPlsLs <- roplsF(dataMatrix, matrix(sampleMetadata[, "genderVc"], ncol = 1), predI
= 1, orthoI = NA)

detach(sacurine)


}


#' PCA, PLS(-DA), and OPLS(-DA)
#'
#' PCA, PLS, and OPLS regression, classification, and cross-validation with the
#' NIPALS algorithm
#'
#' @name opls
#' @rdname opls
#' @aliases opls opls,ExpressionSet-method opls,data.frame-method
#' opls,matrix-method
#' @docType methods
#' @param x Numerical data frame or matrix (observations x variables; NAs are
#' allowed); or ExpressionSet object with non empty assayData, for PCA, and
#' phenoData@data, for (O)PLS(-DA), slots
#' @param y Response to be modelled: Either 1) 'NULL' for PCA (default) or 2) a
#' numerical vector (same length as 'x' row number) for single response (O)PLS,
#' or 3) a numerical matrix (same row number as 'x') for multiple response PLS,
#' 4) a factor (same length as 'x' row number) for (O)PLS-DA, or 5) a character
#' indicating the name of the column of the phenoData@data to be used, when x
#' is an ExpressionSet object. Note that, for convenience, character vectors
#' are also accepted for (O)PLS-DA as well as single column numerical (resp.
#' character) matrix for (O)PLS (respectively (O)PLS-DA). NAs are allowed in
#' numeric responses.
#' @param predI Integer: number of components (predictive componenents in case
#' of PLS and OPLS) to extract; for OPLS, predI is (automatically) set to 1; if
#' set to NA [default], autofit is performed: a maximum of 10 components are
#' extracted until (i) PCA case: the variance is less than the mean variance of
#' all components (note that this rule requires all components to be computed
#' and can be quite time-consuming for large datasets) or (ii) PLS case: either
#' R2Y of the component is < 0.01 (N4 rule) or Q2Y is < 0 (for more than 100
#' observations) or 0.05 otherwise (R1 rule)
#' @param orthoI Integer: number of orthogonal components (for OPLS only); when
#' set to 0 [default], PLS will be performed; otherwise OPLS will be peformed;
#' when set to NA, OPLS is performed and the number of orthogonal components is
#' automatically computed by using the cross-validation (with a maximum of 9
#' orthogonal components).
#' @param algoC Default algorithm is 'svd' for PCA (in case of no missing
#' values in 'x'; 'nipals' otherwise) and 'nipals' for PLS and OPLS; when
#' asking to use 'svd' for PCA on an 'x' matrix containing missing values, NAs
#' are set to half the minimum of non-missing values and a warning is generated
#' @param crossvalI Integer: number of cross-validation segments (default is
#' 7); The number of samples (rows of 'x') must be at least >= crossvalI
#' @param log10L Should the 'x' matrix be log10 transformed? Zeros are set to 1
#' prior to transformation
#' @param permI Integer: number of random permutations of response labels to
#' estimate R2Y and Q2Y significance by permutation testing [default is 20 for
#' single response models (without train/test partition), and 0 otherwise]
#' @param scaleC Character: either no centering nor scaling ('none'),
#' mean-centering only ('center'), mean-centering and pareto scaling
#' ('pareto'), or mean-centering and unit variance scaling ('standard')
#' [default]
#' @param subset Integer vector: indices of the observations to be used for
#' training (in a classification scheme); use NULL [default] for no partition
#' of the dataset; use 'odd' for a partition of the dataset in two equal sizes
#' (with respect to the classes proportions)
#' @param printL Logical: Should informations regarding the data set and the
#' model be printed? [default = TRUE]
#' @param plotL Logical: Should the 'summary' plot be displayed? [default =
#' TRUE]
#' @param .sinkC Character: Name of the file for R output diversion [default =
#' NULL: no diversion]; Diversion of messages is required for the integration
#' into Galaxy
#' @param ... Currently not used.
#' @return An S4 object of class 'opls' containing the following slots:
#' \itemize{
#' \item typeC Character: model type (PCA, PLS, PLS-DA, OPLS, or OPLS-DA)
#' \item descriptionMC Character matrix: Description of the data set (number
#' of samples, variables, etc.)
#' \item modelDF Data frame with the model overview (number of components, R2X, R2X(cum), R2Y, R2Y(cum), Q2, Q2(cum),
#' significance, iterations)
#' \item summaryDF Data frame with the model summary (cumulated R2X, R2Y and Q2); RMSEE is the square root of the mean
#' error between the actual and the predicted responses
#' \item subsetVi Integer vector: Indices of observations in the training data set
#' \item pcaVarVn PCA: Numerical vector of variances of length: predI
#' \item vipVn PLS(-DA): Numerical vector of Variable Importance in Projection; OPLS(-DA): Numerical vector of Variable Importance for Prediction (VIP4,p from Galindo-Prieto et al, 2014)
#' \item orthoVipVn OPLS(-DA): Numerical vector of Variable Importance for Orthogonal Modeling (VIP4,o from Galindo-Prieto et al, 2014)
#' \item xMeanVn Numerical vector: variable means of the 'x' matrix
#' \item xSdVn Numerical vector: variable standard deviations of the 'x' matrix
#' \item yMeanVn (O)PLS: Numerical vector: variable means of the 'y' response (transformed into a dummy matrix in case it is of 'character' mode initially)
#' \item ySdVn (O)PLS: Numerical vector: variable standard deviations of the 'y' response (transformed into a dummy matrix in case it is of 'character' mode initially)
#' \item xZeroVarVi Numerical vector: indices of variables with variance < 2.22e-16 which were excluded from 'x' before building the model
#' \item scoreMN Numerical matrix of x scores (T; dimensions: nrow(x) x predI) X = TP' + E; Y = TC' + F
#' \item loadingMN Numerical matrix of x loadings (P; dimensions: ncol(x) x predI) X = TP' + E
#' \item weightMN (O)PLS: Numerical matrix of x weights (W; same dimensions as loadingMN)
#' \item orthoScoreMN OPLS: Numerical matrix of orthogonal scores (Tortho; dimensions: nrow(x) x number of orthogonal components)
#' \item orthoLoadingMN OPLS: Numerical matrix of orthogonal loadings (Portho; dimensions: ncol(x) x number of orthogonal components)
#' \item orthoWeightMN OPLS: Numerical matrix of orthogonal weights (same dimensions as orthoLoadingMN)
#' \item cMN (O)PLS: Numerical matrix of Y weights (C; dimensions: number of responses or number of classes in case of qualitative response) x number of predictive components; Y = TC' + F
#' \item coMN) (O)PLS: Numerical matrix of Y orthogonal weights; dimensions: number of responses or number of classes in case of qualitative response with more than 2 classes x number of orthogonal components
#' \item uMN (O)PLS: Numerical matrix of Y scores (U; same dimensions as scoreMN); Y = UC' + G
#' \item weightStarMN Numerical matrix of projections (W*; same dimensions as loadingMN); whereas columns of weightMN are derived from successively deflated matrices, columns of weightStarMN relate to the original 'x' matrix: T = XW*; W*=W(P'W)inv
#' \item suppLs List of additional objects to be used internally by the 'print', 'plot', and 'predict' methods
#' }
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @references Eriksson et al. (2006). Multi- and Megarvariate Data Analysis.
#' Umetrics Academy.  Rosipal and Kramer (2006). Overview and recent advances
#' in partial least squares Tenenhaus (1990). La regression PLS : theorie et
#' pratique. Technip.  Wehrens (2011). Chemometrics with R. Springer.  Wold et
#' al. (2001). PLS-regression: a basic tool of chemometrics
#' @examples
#'
#' #### PCA
#'
#' data(foods) ## see Eriksson et al. (2001); presence of 3 missing values (NA)
#' head(foods)
#' foodMN <- as.matrix(foods[, colnames(foods) != "Country"])
#' rownames(foodMN) <- foods[, "Country"]
#' head(foodMN)
#' foo.pca <- opls(foodMN)
#'
#' #### PLS with a single response
#'
#' data(cornell) ## see Tenenhaus, 1998
#' head(cornell)
#' cornell.pls <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
#'                     cornell[, "y"])
#'
#' ## Complementary graphics
#'
#' plot(cornell.pls, typeVc = c("outlier", "predict-train", "xy-score", "xy-weight"))
#'
#' #### PLS with multiple (quantitative) responses
#'
#' data(lowarp) ## see Eriksson et al. (2001); presence of NAs
#' head(lowarp)
#' lowarp.pls <- opls(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
#'                    as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
#'                                       grepl("^st", colnames(lowarp))]))
#'
#' #### PLS-DA
#'
#' data(sacurine)
#' attach(sacurine)
#' sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])
#'
#' #### OPLS-DA
#'
#' sacurine.oplsda <- opls(dataMatrix, sampleMetadata[, "gender"], predI = 1, orthoI = NA)
#'
#' detach(sacurine)
#'
#' @export
setGeneric("opls",
           function(x, ...) standardGeneric("opls"))


#' Tested method for (O)PLS models
#'
#' Returns predictions of the (O)PLS(-DA) model on the out of the box samples
#' (when a 'subset' of samples has been selected when training the model)
#'
#' @aliases tested tested,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Predictions (either a vector, factor, or matrix depending on the y
#' response used for training the model)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' testedorMN <- dataMatrix
#' responseFc <- sampleMetadata[, "gender"]
#'
#' sacurine.plsda <- opls(testedorMN,
#'                        responseFc,
#'                        subset = "odd")
#'
#' trainVi <- getSubsetVi(sacurine.plsda)
#'
#' table(responseFc[trainVi], fitted(sacurine.plsda))
#'
#' detach(sacurine)
#'
#' @rdname tested
#' @export
setGeneric("tested",
           function(object, ...) standardGeneric("tested"))


#' getSummaryDF method for PCA/(O)PLS models
#'
#' Summary of model metrics
#'
#' @aliases getSummaryDF getSummaryDF,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Data frame
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getSummaryDF(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getSummaryDF
#' @export
setGeneric("getSummaryDF",
           function(object, ...) {standardGeneric("getSummaryDF")})


#' getPcaVarVn method for PCA models
#'
#' Variance of the components (score vectors)
#'
#' @aliases getPcaVarVn getPcaVarVn,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Numeric vector with the same length as the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.pca <- opls(dataMatrix)
#'
#' getPcaVarVn(sacurine.pca)
#'
#' detach(sacurine)
#'
#' @rdname getPcaVarVn
#' @export
setGeneric("getPcaVarVn",
           function(object, ...) {standardGeneric("getPcaVarVn")})


#' getScoreMN method for PCA/(O)PLS(-DA) models
#'
#' (Orthogonal) scores of the (O)PLS(-DA) model
#'
#' @aliases getScoreMN getScoreMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal score matrix be returned
#' (default is FALSE and the predictive score matrix is returned)
#' @param ... Currently not used.
#' @return Numeric matrix with a number of rows equal to the number of samples
#' and a number of columns equal to the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getScoreMN(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getScoreMN
#' @export
setGeneric("getScoreMN",
           function(object, ...) {standardGeneric("getScoreMN")})


#' getLoadingMN method for PCA/(O)PLS(-DA) models
#'
#' (Orthogonal) loadings of the PCA/(O)PLS(-DA) model
#'
#' @aliases getLoadingMN getLoadingMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal loading matrix be returned
#' @param ... Currently not used.
#' (default is FALSE and the predictive loading matrix is returned)
#' @return Numeric matrix with a number of rows equal to the number of
#' variables and a number of columns equal to the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getLoadingMN(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getLoadingMN
#' @export
setGeneric("getLoadingMN",
           function(object, ...) {standardGeneric("getLoadingMN")})


#' getWeightMN method for (O)PLS(-DA) models
#'
#' (Orthogonal) weights of the (O)PLS(-DA) model
#'
#'
#' @aliases getWeightMN getWeightMN,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal weight matrix be returned
#' @param ... Currently not used.
#' (default is FALSE and the predictive weight matrix is returned)
#' @return Numeric matrix with a number of rows equal to the number of
#' variables and a number of columns equal to the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getWeightMN(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getWeightMN
#' @export
setGeneric("getWeightMN",
           function(object, ...) {standardGeneric("getWeightMN")})


#' getVipVn method for (O)PLS(-DA) models
#'
#' (Orthogonal) VIP of the (O)PLS(-DA) model
#'
#'
#' @aliases getVipVn getVipVn,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param orthoL Logical: Should the orthogonal VIP be returned (default is
#' FALSE and the predictive VIP is returned)
#' @param ... Currently not used.
#' @return Numeric vector with a length equal to the number of variables and a
#' number of columns equal to the number of components
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @references Galindo-Prieto B., Eriksson L. and Trygg J. (2014). Variable
#' influence on projection (VIP) for orthogonal projections to latent
#' structures (OPLS). Journal of Chemometrics 28, 623-632.
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' sacurine.plsda <- opls(dataMatrix,
#'                        sampleMetadata[, "gender"])
#'
#' getVipVn(sacurine.plsda)
#'
#' detach(sacurine)
#'
#' @rdname getVipVn
#' @export
setGeneric("getVipVn",
           function(object, ...) {standardGeneric("getVipVn")})


#' getSubsetVi method for (O)PLS(-DA) models
#'
#' Extracts the indices of the samples used for building the model (when a
#' subset argument has been specified)
#'
#' @aliases getSubsetVi getSubsetVi,opls-method
#' @param object An S4 object of class \code{opls}, created by \code{opls}
#' function.
#' @param ... Currently not used.
#' @return Integer vector with the indices of the samples used for training
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'
#' data(sacurine)
#' attach(sacurine)
#'
#' predictorMN <- dataMatrix
#' responseFc <- sampleMetadata[, "gender"]
#'
#' sacurine.plsda <- opls(predictorMN,
#'                        responseFc,
#'                        subset = "odd")
#'
#' trainVi <- getSubsetVi(sacurine.plsda)
#'
#' table(responseFc[trainVi], fitted(sacurine.plsda))
#'
#' detach(sacurine)
#'
#' @rdname getSubsetVi
#' @export
setGeneric("getSubsetVi",
           function(object, ...) {standardGeneric("getSubsetVi")})


#' Checking the consistency of an ExpressionSet instance with W4M format
#'
#' @param eset An S4 object of class ExpressionSet.
#' @param ... Currently not used.
#' @return Invisible TRUE logical in case of success (otherwise generates an
#' error)
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#' sacSet <- fromW4M(file.path(path.package("ropls"), "extdata"))
#' print(checkW4M(sacSet))
#' @rdname checkW4M
#' @export
setGeneric("checkW4M", function(eset, ...) {standardGeneric("checkW4M")})


#' Exporting ExpressionSet instance into 3 tabulated files.
#'
#' The 3 .tsv files are written with the indicated \code{file} prefix, and
#' '_dataMatrix.tsv', '_sampleMetadata.tsv', and '_variableMetadata.tsv'
#' suffices, respectively. Note that the \code{dataMatrix} is transposed before
#' export (e.g., the samples are written column wise in the 'dataMatrix.tsv'
#' exported file).
#'
#' @param eset An S4 object of class \code{ExpressionSet}
#' function.
#' @param filePrefixC Character: common prefix (including repository full path)
#' of the three file names: for example, the 'c:/mydata/setname' value will
#' result in writting the 'c:/mydata/setname_dataMatrix.tsv',
#' 'c:/mydata/setname_sampleMetadata.tsv', and
#' 'c:/mydata/setname_variableMetadata.tsv' files.
#' @param verboseL Logical: should comments be printed?
#' @param ... Currently not used.
#' @return No object returned.
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples
#'  sacSet <- fromW4M(file.path(path.package("ropls"), "extdata"))
#'  toW4M(sacSet)
#' @rdname toW4M
#' @export
setGeneric("toW4M", function(eset, ...) standardGeneric("toW4M"))

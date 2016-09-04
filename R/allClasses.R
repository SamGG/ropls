#' Class "opls"
#'
#' An S4 class to store PCA and (O)PLS(-DA) models: Objects can be created by calls of the form
#' \code{new("opls", ...)} or by calling the \code{opls} function
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#'
#' @slot typeC character: model type (PCA, PLS, PLS-DA, OPLS, or OPLS-DA)
#' @slot descriptionMC character matrix: Description of the data set (number of samples, variables, etc.)
#' @slot modelDF data frame with the model overview (number of components, R2X, R2X(cum), R2Y, R2Y(cum), Q2, Q2(cum), significance, iterations)
#' @slot summaryDF data frame with the model summary (cumulated R2X, R2Y and Q2); RMSEE is the square root of the mean error between the actual and the predicted responses
#' @slot subsetVi Integer vector: Indices of observations in the training data set
#' @slot pcaVarVn PCA: Numerical vector of variances of length: predI
#' @slot vipVn PLS(-DA): Numerical vector of Variable Importance in Projection; OPLS(-DA): Numerical vector of Variable Importance for Prediction (VIP4,p from Galindo-Prieto et al, 2014)
#' @slot orthoVipVn OPLS(-DA): Numerical vector of Variable Importance for Orthogonal Modeling (VIP4,o from Galindo-Prieto et al, 2014)
#' @slot coefficientMN (O)PLS(-DA): Numerical matrix of regression coefficients (B; dimensions: ncol(x) x number of responses; B = W*C' and Y = XB + F
#' @slot xMeanVn Numerical vector: variable means of the 'x' matrix
#' @slot xSdVn Numerical vector: variable standard deviations of the 'x' matrix
#' @slot yMeanVn (O)PLS: Numerical vector: variable means of the 'y' response (transformed into a dummy matrix in case it is of 'character' mode initially)
#' @slot ySdVn (O)PLS: Numerical vector: variable standard deviations of the 'y' response (transformed into a dummy matrix in case it is of 'character' mode initially)
#' @slot xZeroVarVi Numerical vector: indices of variables with variance < 2.22e-16 which were excluded from 'x' before building the model
#' @slot scoreMN Numerical matrix of x scores (T; dimensions: nrow(x) x predI) X = TP' + E; Y = TC' + F
#' @slot loadingMN Numerical matrix of x loadings (P; dimensions: ncol(x) x predI) X = TP' + E
#' @slot weightMN (O)PLS: Numerical matrix of x weights (W; same dimensions as loadingMN)
#' @slot orthoScoreMN OPLS: Numerical matrix of orthogonal scores (Tortho; dimensions: nrow(x) x number of orthogonal components)
#' @slot orthoLoadingMN OPLS: Numerical matrix of orthogonal loadings (Portho; dimensions: ncol(x) x number of orthogonal components)
#' @slot orthoWeightMN OPLS: Numerical matrix of orthogonal weights (same dimensions as orthoLoadingMN)
#' @slot cMN (O)PLS: Numerical matrix of Y weights (C); dimensions: number of responses or number of classes in case of qualitative response with more than 2 classes x number of predictive components; Y = TC' + F
#' @slot coMN (O)PLS: Numerical matrix of Y orthogonal weights; dimensions: number of responses or number of classes in case of qualitative response with more than 2 classes x number of orthogonal components
#' @slot uMN (O)PLS: Numerical matrix of Y scores (U; same dimensions as scoreMN); Y = UC' + G
#' @slot weightStarMN Numerical matrix of projections (W*; same dimensions as loadingMN); whereas columns of weightMN are derived from successively deflated matrices, columns of weightStarMN relate to the original 'x' matrix: T = XW*; W*=W(P'W)inv
#' @slot suppLs List of additional objects to be used internally by the 'print', 'plot', and 'predict' methods
#' @name opls-class
#' @rdname opls-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("opls", ...)} or by calling the \code{opls} function
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @seealso \code{\link{opls}}
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
#' @exportClass opls
setClass(Class = "opls",
         representation = representation(typeC = "character",
             descriptionMC = "matrix",
             modelDF = "data.frame",
             summaryDF = "data.frame",
             subsetVi = "numeric",
             pcaVarVn = "numeric",
             vipVn = "numeric",
             orthoVipVn = "numeric",
             coefficientMN = "matrix",
             xMeanVn = "numeric",
             xSdVn = "numeric",
             yMeanVn = "numeric",
             ySdVn = "numeric",
             xZeroVarVi = "numeric",
             scoreMN = "matrix",
             loadingMN = "matrix",
             weightMN = "matrix",
             orthoScoreMN = "matrix",
             orthoLoadingMN = "matrix",
             orthoWeightMN = "matrix",
             cMN = "matrix",
             uMN = "matrix",
             weightStarMN = "matrix",
             coMN = "matrix",
             suppLs = "list"))

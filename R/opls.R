opls <-
function (x, ...)
    UseMethod ("opls")

opls.default <- function(x,
                         y = NULL,
                         predI = NA,
                         orthoI = 0,

                         algoC = c("default", "nipals", "svd")[1],
                         crossvalI = 7,
                         log10L = FALSE,
                         permI = 20,
                         scaleC = c("none", "center", "pareto", "standard")[4],
                         subset = NULL,

                         printL = TRUE,
                         plotL = TRUE,

                         .sinkC = NULL,
                         ...) {

    if(!is.null(.sinkC)) ##  Diversion of messages is required for the integration into Galaxy
        sink(.sinkC, append = TRUE)

    ## Checking arguments
    ##-------------------

    ## x -> xMN

    if(is.data.frame(x)) {
        if(!all(sapply(x, data.class) == "numeric")) {
            stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
        } else
            x <- as.matrix(x)
    } else if(is.matrix(x)) {
        if(mode(x) != "numeric")
            stop("'x' matrix must be of 'numeric' mode", call. = FALSE)
    } else
        stop("'x' must be either a data.frame or a matrix", call. = FALSE)

    if(any(apply(x, 2, function(colVn) all(is.na(colVn)))))
        stop("'x' contains columns with 'NA' only", call. = FALSE)

    xMN <- x

    ## y -> yMCN

    yLevelVc <- NULL

    if(!is.null(y)) {

        if(is.vector(y)) {

            if(!(mode(y) %in% c("character", "numeric")))
                stop("'y' vector must be of 'character' or 'numeric' type", call. = FALSE)

            if(length(y) != nrow(xMN))
                stop("'y' vector length must be equal to the number of rows of 'x'", call. = FALSE)

            yMCN <- matrix(y, ncol = 1)

        } else if(is.factor(y)) {

            if(length(y) != nrow(xMN))
                stop("'y' factor length must be equal to the number of rows of 'x'", call. = FALSE)

            yLevelVc <- levels(y)

            yMCN <- matrix(as.character(y), ncol = 1)

        } else if(is.matrix(y)) {

            if(!(mode(y) %in% c("character", "numeric")))
                stop("'y' matrix must be of 'character' or 'numeric' type", call. = FALSE)

            if(nrow(y) != nrow(xMN))
                stop("'x' and 'y' must have the same number of rows", call. = FALSE)

            if(ncol(y) > 1 && mode(y) != "numeric")
                stop("Multiple response 'y' matrices must be of numeric type", call. = FALSE)

            yMCN <- y

        } else
            stop("'y' must be either a vector, a factor, or a matrix", call. = FALSE)

    } else
        yMCN <- NULL

    ## NA in Y only possible for multi-response regression (i.e., Y is a numeric matrix)

    if(!is.null(yMCN) &&
       ncol(yMCN) == 1 &&
       any(is.na(drop(yMCN))))
        stop("In case of single response modeling, 'y' must not contain missing values", call. = FALSE)

    if(!is.logical(log10L))
        stop("'log10L' must be a logical", call. = FALSE)

    if(permI < 0 || (permI - floor(permI)) > 1e-10)
        stop("'permI' must be an integer", call. = FALSE)

    if(permI > 0 && (is.null(yMCN) || ncol(yMCN) > 1)) {
        ## warning("Permutation testing available for single response (O)PLS(-DA) models only", call. = FALSE)
        permI <- 0
    }

    if(permI > 0 && !is.null(subset)) {
        permI <- 0
        warning("'permI' set to 0 because train/test partition is selected", call. = FALSE)
    }

    if(!(algoC %in% c('default', 'nipals', 'svd')))
        stop("'algoC' must be either 'default', 'nipals', or 'svd'", call. = FALSE)

    if(algoC == "default")
        algoC <- ifelse(is.null(yMCN) && !any(is.na(c(xMN))), "svd", "nipals")

    if(!is.null(yMCN) && algoC != "nipals")
        stop("'nipals' algorithm must be used for (O)PLS(-DA)", call. = FALSE)

    if((is.na(orthoI) || orthoI > 0) && is.null(yMCN))
        stop("'y' response cannot be NULL for OPLS(-DA) modeling", call. = FALSE)

    if(!is.null(yMCN)) {
        if(is.na(orthoI) || orthoI > 0) {
            if(ncol(yMCN) > 1) {
                stop("OPLS regression only available for a single 'y' response", call. = FALSE)
            } else if(mode(yMCN) == "character" && length(unique(drop(yMCN))) > 2)
                stop("OPLS-DA only available for binary classification (use PLS-DA for multiple classes)", call. = FALSE)
        }
    }

    if(is.na(orthoI) || orthoI > 0)
        if(is.na(predI) || predI > 1) {
            predI <- 1
            warning("OPLS: number of predictive components ('predI' argument) set to 1", call. = FALSE)
        }

    if(!is.na(predI) && !is.na(orthoI) && ((predI + orthoI) > min(dim(xMN))))
        stop("The sum of 'predI' (", predI, ") and 'orthoI' (", orthoI, ") exceeds the minimum dimension of the 'x' data matrix (", min(dim(xMN)), ")" , call. = FALSE)

    if(!(length(scaleC) == 1 && scaleC %in% c('none', 'center', 'pareto', 'standard')))
        stop("'scaleC' must be either 'none', 'center', 'pareto', or 'standard'", call. = FALSE)

    if(!is.null(subset) && (is.null(yMCN) || ncol(yMCN) > 1))
        stop("train/test partition with 'subset' only available for (O)PLS(-DA) models of a single 'y' response", call. = FALSE)

    if(!is.null(subset) &&
       !(mode(subset) == 'character' && subset == 'odd') &&
       !all(subset %in% 1:nrow(xMN)))
        stop("'subset' must be either set to 'odd' or an integer vector of 'x' row numbers", call. = FALSE)

    if(crossvalI > nrow(xMN))
        stop("'crossvalI' must be less than the row number of 'x'", call. = FALSE)


    ## Constants
    ##----------

    epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16


    ## Character to numeric convertion function (for classification)
    ##--------------------------------------------------------------

    if(!is.null(yMCN) && mode(yMCN) == "character") {

        if(!is.null(yLevelVc)) {
            claVc <- yLevelVc
        } else
            claVc <- sort(unique(drop(yMCN)))

        if(length(claVc) == 2) {
            ## binary response kept as a single vector for OPLS-DA computations
            .char2numF <- function(inpMCN,
                                   c2nL = TRUE) {

                if(c2nL) {

                    outMCN <- inpMCN == claVc[2]
                    mode(outMCN) <- "numeric"

                } else {

                    outMCN <- matrix(claVc[as.numeric(inpMCN > 0.5) + 1],
                                     ncol = 1,
                                     dimnames = dimnames(inpMCN))

                }

                return(outMCN)

            }

        } else
            .char2numF <- function(inpMCN,
                                   c2nL = TRUE) {

                if(c2nL) {

                    outMCN  <- t(sapply(drop(inpMCN),
                                        function(claC) as.numeric(claVc == claC)))
                    colnames(outMCN) <- claVc


                } else {

                    outMCN <- t(t(apply(inpMCN, 1,
                                        function(rowVn) claVc[which(rowVn == max(rowVn))[1]])))
                    colnames(outMCN) <- "y1"

                }

                return(outMCN)

            }

    } else
        .char2numF <- NULL


    ##------------------------------------
    ##   Computations
    ##------------------------------------


    ## rownames and colnames

    if(is.null(rownames(xMN)))
        if(!is.null(yMCN) && !is.null(rownames(yMCN))) {
            rownames(xMN) <- rownames(yMCN)
        } else
            rownames(xMN) <- paste0("s", 1:nrow(xMN))
    if(is.null(colnames(xMN)))
        colnames(xMN) <- paste0("x", 1:ncol(xMN))

    if(!is.null(yMCN)) {
        if(is.null(rownames(yMCN)))
            rownames(yMCN) <- rownames(xMN)
        if(is.null(colnames(yMCN)))
            colnames(yMCN) <- paste0("y", 1:ncol(yMCN))
    }

    ## Log10 transformation
    ##---------------------

    if(log10L)
        xMN <- .log10F(xMN)


    ## Test indices
    ##-------------

    obsIniVi <- 1:nrow(xMN)

    if(!is.null(subset)) {
        if(length(subset) == 1 && subset == "odd") {

            if(mode(yMCN) == "numeric")
                subset <- seq(1, nrow(xMN), by = 2)
            else {
                subset <- integer()
                for(claC in unique(drop(yMCN)))
                    subset <- c(subset,
                                which(drop(yMCN) == claC)[seq(1, sum(drop(yMCN) == claC), by = 2)])
                subset <- sort(subset)
            }
        }
        if(crossvalI > length(subset))
            stop("'crossvalI' must be less than the number of samples in the subset", call. = FALSE)
    }


    ## Filtering out zero variance variables
    ##--------------------------------------

    xVarIndLs <- list()
    xVarIndLs[[1]] <- 1:nrow(xMN)

    if(!is.null(subset)) {
        xVarIndLs[[1]] <- subset
    } ## else if(!is.null(yMCN) && ncol(yMCN) == 1 && nrow(xMN) >= 2 * crossvalI)
      ##   for(cvkN in 1:crossvalI)
      ##       xVarIndLs <- c(xVarIndLs, list(setdiff(1:nrow(xMN), cvaOutLs[[cvkN]])))

    xVarVarLs <- lapply(xVarIndLs,
                        function(xVarVi) {
                            apply(xMN[xVarVi, , drop = FALSE],
                                  2,
                                  function(colVn) var(colVn, na.rm = TRUE))
                        })

    xZeroVarVi <- integer()
    for(k in 1:length(xVarVarLs))
        xZeroVarVi <- union(xZeroVarVi, which(xVarVarLs[[k]] < epsN))

    if(length(xZeroVarVi) > 0) {
        names(xZeroVarVi) <- colnames(xMN)[xZeroVarVi]
        xMN <- xMN[, -xZeroVarVi, drop = FALSE]
        warning("The variance of the ",
                length(xZeroVarVi),
                " following variables is less than ",
                signif(epsN, 2),
                " in the full or partial (cross-validation) dataset: these variables will be removed:\n",
                paste(names(xZeroVarVi), collapse = ", "),
                call. = FALSE)
    }


    ## Core
    ##-----

    opLs <- .coreF(xMN = xMN,
                   yMCN = yMCN,
                   orthoI = orthoI,
                   predI = predI,
                   scaleC = scaleC,
                   algoC = algoC,
                   crossvalI = crossvalI,
                   subset = subset,
                   .char2numF = .char2numF)


    if(is.null(opLs[["suppLs"]][["yMCN"]])) {
        opLs[["typeC"]] <- "PCA"
    } else {
        if(ncol(opLs[["suppLs"]][["yMCN"]]) > 1 || mode(opLs[["suppLs"]][["yMCN"]]) == "numeric")
            opLs[["typeC"]] <- "PLS"
        else
            opLs[["typeC"]] <- "PLS-DA"
    }
    if(opLs[["summaryDF"]][, "ort"] > 0)
        opLs[["typeC"]] <- paste("O", opLs[["typeC"]], sep = "")


    opLs[["xZeroVarVi"]] <- xZeroVarVi
    opLs[["suppLs"]][["yLevelVc"]] <- yLevelVc

    ## Permutation testing (Szymanska et al, 2012)

    if(permI > 0) {

        modSumVc <- colnames(opLs[["summaryDF"]])

        permMN <- matrix(0,
                         nrow = 1 + permI,
                         ncol = length(modSumVc),
                         dimnames = list(NULL, modSumVc))

        perSimVn <- numeric(1 + permI)
        perSimVn[1] <- 1


        permMN[1, ] <- as.matrix(opLs[["summaryDF"]])

        for(k in 1:permI) {

            yVcn <- drop(opLs[["suppLs"]][["yMCN"]])
            if(is.null(opLs[["subset"]])) {
                yPerVcn <- sample(yVcn)
            } else {
                yPerVcn <- numeric(nrow(xMN))
                refVi <- opLs[["subset"]]
                tesVi <- setdiff(1:nrow(xMN), refVi)
                yPerVcn[refVi] <- sample(yVcn[refVi])
                yPerVcn[tesVi] <- yVcn[tesVi]
            }
            yPerMCN <- matrix(yPerVcn, ncol = 1)

            perLs <- .coreF(xMN = xMN,
                            yMCN = yPerMCN,
                            orthoI = opLs[["summaryDF"]][, "ort"],
                            predI = opLs[["summaryDF"]][, "pre"],
                            scaleC = scaleC,
                            algoC = algoC,
                            crossvalI = crossvalI,
                            subset = opLs[["subset"]],
                            .char2numF = .char2numF)

            permMN[1 + k, ] <- as.matrix(perLs[["summaryDF"]])

            perSimVn[1 + k] <- .similarityF(opLs[["suppLs"]][["yMCN"]], yPerMCN,
                                            .char2numF = .char2numF,
                                            charL = mode(opLs[["suppLs"]][["yMCN"]]) == "character")

        }

        permMN <- cbind(permMN, sim = perSimVn)

        perPvaVn <- c(pR2Y = (1 + length(which(permMN[-1, "R2Y(cum)"] >= permMN[1, "R2Y(cum)"]))) / (nrow(permMN) - 1),
                      pQ2 = (1 + length(which(permMN[-1, "Q2(cum)"] >= permMN[1, "Q2(cum)"]))) / (nrow(permMN) - 1))
        opLs[["summaryDF"]][, "pR2Y"] <- perPvaVn["pR2Y"]
        opLs[["summaryDF"]][, "pQ2"] <- perPvaVn["pQ2"]

        opLs[["suppLs"]][["permMN"]] <- permMN

    }

    ##------------------------------------
    ##   Numerical results
    ##------------------------------------

    opLs[["descriptionMC"]] <- rbind(samples = ifelse(is.null(subset),
                                         nrow(xMN),
                                         length(subset)),
                                     X_variables = ncol(xMN),
                                     near_zero_excluded_X_variables = length(opLs[["xZeroVarVi"]]))

    totN <- length(c(xMN))
    nasN <- sum(is.na(c(xMN)))

    if(!is.null(opLs[["suppLs"]][["yMCN"]])) {

        opLs[["descriptionMC"]] <- rbind(opLs[["descriptionMC"]],
                                         Y_variables = ncol(opLs[["suppLs"]][["yMCN"]]))
        totN <- totN + length(c(opLs[["suppLs"]][["yMCN"]]))
        nasN <- nasN + sum(is.na(c(opLs[["suppLs"]][["yMCN"]])))

    }

    opLs[["descriptionMC"]] <- rbind(opLs[["descriptionMC"]],
                                     missing_values = paste0(nasN, " (", round(nasN / totN * 100), "%)"))

    ## Raw summary
    ##------------

    opLs[["topLoadI"]] <- 3

    if(ncol(xMN) > opLs[["topLoadI"]]) {
        xVarVn <- apply(xMN, 2, var)
        names(xVarVn) <- 1:length(xVarVn)
        xVarVn <- sort(xVarVn)
        xVarSorVin <- as.numeric(names(xVarVn[seq(1, length(xVarVn), length = opLs[["topLoadI"]])]))
        opLs[["suppLs"]][["xSubIncVarMN"]] <- xMN[, xVarSorVin, drop = FALSE]
    } else
        opLs[["suppLs"]][["xSubIncVarMN"]] <- xMN

    if(ncol(xMN) <= 100) {

        xCorMN <- cor(xMN, use = "pairwise.complete.obs")
        xCorMN[lower.tri(xCorMN, diag = TRUE)] <- 0

        if(ncol(xMN) > opLs[["topLoadI"]]) {

            xCorNexDF <- which(abs(xCorMN) >= sort(abs(xCorMN), decreasing = TRUE)[opLs[["topLoadI"]] + 1],
                               arr.ind = TRUE)

            xCorDisMN <- matrix(0,
                                nrow = nrow(xCorNexDF),
                                ncol = nrow(xCorNexDF),
                                dimnames = list(colnames(xMN)[xCorNexDF[, "row"]],
                                    colnames(xMN)[xCorNexDF[, "col"]]))

            for(k in 1:nrow(xCorDisMN))
                xCorDisMN[k, k] <- xCorMN[xCorNexDF[k, "row"], xCorNexDF[k, "col"]]

        } else
            xCorDisMN <- xCorMN

        opLs[["suppLs"]][["xCorMN"]] <- xCorDisMN

        rm(xCorDisMN)

    }

    ## fitted
    ##-------

    if(!is.null(opLs[["suppLs"]][["yPreMN"]])) {

        if(mode(opLs[["suppLs"]][["yMCN"]]) == "character") {

            yPredMCN <- opLs[["suppLs"]][[".char2numF"]](opLs[["suppLs"]][["yPreMN"]],
                                              c2nL = FALSE)
            if(is.vector(y)) {
                fitted <- c(yPredMCN)
                names(fitted) <- rownames(yPredMCN)
            } else if(is.factor(y)) {
                fitted <- c(yPredMCN)
                names(fitted) <- rownames(yPredMCN)
                fitted <- factor(fitted, levels = yLevelVc)
            } else if(is.matrix(y)) {
                fitted <- yPredMCN
            } else
                stop() ## this case should not happen

        } else {

            yPredMCN <- opLs[["suppLs"]][["yPreMN"]]

            if(is.vector(y)) {
                fitted <- c(yPredMCN)
                names(fitted) <- rownames(yPredMCN)
            } else if(is.matrix(y)) {
                fitted <- yPredMCN
            } else
                stop() ## this case should not happen

        }

        opLs[["fitted"]] <- fitted

    }

    ## tested (when subset is non null)
    ##---------------------------------

    if(!is.null(opLs[["suppLs"]][["yTesMN"]])) {

        if(mode(opLs[["suppLs"]][["yMCN"]]) == "character") {
            opLs[["yTestMCN"]] <- opLs[["suppLs"]][[".char2numF"]](opLs[["suppLs"]][["yTesMN"]],
                                                                   c2nL = FALSE)
            if(is.vector(y)) {
                tested <- c(opLs[["yTestMCN"]])
                names(tested) <- rownames(opLs[["yTestMCN"]])
            } else if(is.factor(y)) {
                tested <- c(opLs[["yTestMCN"]])
                names(tested) <- rownames(opLs[["yTestMCN"]])
                tested <- factor(tested, levels = yLevelVc)
            } else if(is.matrix(y)) {
                tested <- opLs[["yTestMCN"]]
            } else
                stop() ## this case should not happen
        } else {

            opLs[["yTestMCN"]] <- opLs[["suppLs"]][["yTesMN"]]

            if(is.vector(y)) {
                tested <- c(opLs[["yTestMCN"]])
                names(tested) <- rownames(opLs[["yTestMCN"]])
            } else if(is.matrix(y)) {
                tested <- opLs[["yTestMCN"]]
            } else
                stop() ## this case should not happen

        }

        opLs[["tested"]] <- tested

    }

    ## residuals (O)PLS(-DA)
    ##----------------------

    if(opLs[["typeC"]] != "PCA") {

        if(is.null(opLs[["subset"]])) {
            trainVi <- 1:length(opLs[["fitted"]])
        } else
            trainVi <- opLs[["subset"]]

        if(mode(opLs[["suppLs"]][["yMCN"]]) == "numeric") {
            opLs[["residuals"]] <- y[trainVi] - opLs[["fitted"]]
        } else
            opLs[["residuals"]] <- as.numeric(as.character(c(y)[trainVi]) != as.character(c(opLs[["fitted"]])))

    }

    ## Defining the 'opls' class
    ##--------------------------

    class(opLs) <- "opls"

    ## Printing
    ##---------

    if(printL) {
        print.opls(opLs)
        warnings()
    }

    ## Plotting
    ##---------

    if(plotL)
        plot(opLs, typeVc = "summary")

    ## Closing connection
    ##-------------------

    if(!is.null(.sinkC)) ## Used in the Galaxy module
        sink()

    ## Returning
    ##----------

    return(invisible(opLs))


} ## end of opls.default



## Core algorithms for PCA, PLS(-DA), and OPLS(-DA)
.coreF <- function(xMN,
                   yMCN,
                   orthoI,
                   predI,
                   scaleC,
                   algoC,
                   crossvalI,
                   subset,
                   .char2numF = .char2numF){

    epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16

    ##------------------------------------
    ##   Initialization
    ##------------------------------------

    varVn <- NULL      ## PCA only
    yMeanVn <- NULL    ## (O)PLS only
    ySdVn <- NULL      ## (O)PLS only
    wMN <- NULL        ## (O)PLS only
    cMN <- NULL        ## (O)PLS only
    uMN <- NULL        ## (O)PLS only
    rMN <- NULL        ## (O)PLS only
    bMN <- NULL        ## (O)PLS only
    vipVn <- NULL      ## (O)PLS only
    yPreMN <- NULL     ## (O)PLS only
    yTesMN <- NULL     ## (O)PLS only
    toMN <- NULL       ## OPLS only
    poMN <- NULL       ## OPLS only
    woMN <- NULL       ## OPLS only
    coMN <- NULL       ## OPLS only
    orthoVipVn <- NULL ## OPLS only


    ## Missing values
    ##---------------


    naxVi <- which(is.na(c(xMN)))
    naxL <- length(naxVi) > 0
    nayVi <- integer()
    nayL <- FALSE

    if(algoC == "svd" && length(which(is.na(c(xMN)))) > 0) {
        minN <- min(c(xMN[!is.na(xMN)])) / 2
        xMN[is.na(xMN)] <- minN
        warning("Missing values set to ", round(minN, 1), " (half minimum value) for 'svd' algorithm to be used", call. = FALSE)
    }

    if(!is.null(yMCN)) {
        nayVi <- which(is.na(c(yMCN)))
        nayL <- length(nayVi) > 0
        if(nayL && ncol(yMCN) == 1)
            stop("Current implementation does not handle missing values in single response models") ## TO DO
    }

    ## yMCN 'character' to 'numeric' conversion

    yMN <- yMCN

    if(!is.null(yMCN)) {

        if(mode(yMCN) == "character")
            yMN <- .char2numF(yMCN)

        ## training and a test partition

        if(!is.null(subset)) {

            xTesMN <- xMN[-subset, , drop = FALSE]
            xMN <- xMN[subset, , drop = FALSE]
            yMN <- yMN[subset, , drop = FALSE]

        }

    }


    ## Observation names
    ##------------------

    if(!is.null(rownames(xMN))) {
        obsNamVc <- rownames(xMN)
    } else
        obsNamVc <- as.character(1:nrow(xMN))


    ## Autofit
    ##---------------

    autNcoL <- autNcpL <- FALSE
    autMaxN <- min(c(10, dim(xMN)))

    if(is.na(orthoI)) {
        if(autMaxN == 1) {
            orthoI <- 0
            predI <- 1
            warning("The data contain a single variable (or sample): A PLS model with a single component will be built", call. = FALSE)
        } else {
            orthoI <- autMaxN - 1
            predI <- 1
            autNcoL <- TRUE
        }
    }
    if(is.na(predI)) {
        if(orthoI > 0) {
            if(autMaxN == 1) {
                orthoI <- 0
                warning("The data contain a single variable (or sample): A PLS model with a single component will be built", call. = FALSE)
            } else
                warning("OPLS(-DA): The number of predictive component is set to 1 for a single response model", call. = FALSE)
            predI <- 1
            if((predI + orthoI) > min(dim(xMN)))
                stop("The sum of 'predI' (", predI, ") and 'orthoI' (", orthoI, ") exceeds the minimum dimension of the 'x' data matrix (", min(dim(xMN)), ")" , call. = FALSE)
        } else {
            predI <- autMaxN
            autNcpL <- TRUE
        }
    }


    ##------------------------------------
    ##   Preprocessing
    ##------------------------------------


    ## X variable variances
    ##---------------------

    xVarVn <- apply(xMN, 2, function(colVn) var(colVn, na.rm = TRUE))


    ## X-Scaling
    ##---------------

    xMeanVn <- apply(xMN, 2, function(colVn) mean(colVn, na.rm = TRUE))

    switch(scaleC,
           none = {
               xMeanVn <- rep(0, ncol(xMN))
               xSdVn <- rep(1, times = ncol(xMN))
           },
           center = {
               xSdVn <- rep(1, times = ncol(xMN))
           },
           pareto = {
               xSdVn <- apply(xMN, 2, function(colVn) sqrt(sd(colVn, na.rm = TRUE)))
           },
           standard = {
               xSdVn <- apply(xMN, 2, function(colVn) sd(colVn, na.rm = TRUE))
           })

    xMN <- scale(xMN, center = xMeanVn, scale = xSdVn)


    if(!is.null(colnames(xMN))) {
        xvaNamVc <- colnames(xMN)
    } else
        xvaNamVc <- paste("x", 1:ncol(xMN), sep = "")

    preNamVc <- paste("p", 1:predI, sep = "")

    pMN <- matrix(0,
                  nrow = ncol(xMN),
                  ncol = predI,
                  dimnames = list(xvaNamVc, preNamVc))

    tMN <- uMN <- matrix(0,
                         nrow = nrow(xMN),
                         ncol = predI,
                         dimnames = list(obsNamVc, preNamVc))

    ssxTotN <- sum(xMN^2, na.rm = TRUE)


    if(is.null(yMCN)) {


        ##------------------------------------
        ##   PCA
        ##------------------------------------


        varVn <- numeric(predI)
        names(varVn) <- preNamVc
        vSumVn <- sum(apply(xMN, 2, function(y) var(y, na.rm = TRUE))) ## xMN is centered

        modelDF <- as.data.frame(matrix(0,
                                      nrow = predI,
                                      ncol = 3,
                                      dimnames = list(preNamVc, c("R2X", "R2X(cum)", "Iter."))))

        switch(algoC,

               nipals = {

                   ## NIPALS
                   ##-------

                   xOldMN <- xMN

                   for(hN in 1:predI) {

                       iteN <- 1
                       tOldVn <- xOldMN[, 1]
                       pOldVn <- rep(0, ncol(xMN))

                       repeat {

                           if(naxL) {
                               pNewVn <- numeric(length(pOldVn))
                               for(j in 1:length(pNewVn)) {
                                   comVl <- complete.cases(xOldMN[, j]) &
                                       complete.cases(tOldVn)
                                   pNewVn[j] <- crossprod(xOldMN[comVl, j], tOldVn[comVl]) / drop(crossprod(tOldVn[comVl]))
                               }
                           } else {
                               pNewVn <- crossprod(xOldMN, tOldVn) / drop(crossprod(tOldVn))
                           }

                           pNewVn <- pNewVn / sqrt(drop(crossprod(pNewVn)))

                           if(naxL) {
                               tNewVn <- numeric(length(tOldVn))
                               for(i in 1:length(tNewVn)) {
                                   comVl <- complete.cases(xOldMN[i, ])
                                   tNewVn[i] <- crossprod(xOldMN[i, comVl], pNewVn[comVl])
                               }
                           } else {
                               tNewVn <- xOldMN %*% pNewVn
                           }


                           if(sqrt(drop(crossprod(pNewVn - pOldVn))) < 1e-6 || iteN > 100) {

                               break

                           } else {

                               tOldVn <- tNewVn
                               pOldVn <- pNewVn
                               iteN <- iteN + 1

                           }

                       }

                       tMN[, hN] <- tNewVn
                       pMN[, hN] <- pNewVn
                       varVn[hN] <- 1 / (min(nrow(xMN), ncol(xMN)) - 1) * drop(crossprod(tNewVn))
                       xOldMN <- xOldMN - tcrossprod(tNewVn, pNewVn)

                       modelDF[hN, "R2X"] <- sum(tcrossprod(tMN[, hN], pMN[, hN])^2) / ssxTotN
                       modelDF[hN, "Iter."] <- iteN

                   } ## for(hN in 1:predI) {

               }, ## nipals

               svd = {

                   ## SVD algorithm
                   ## PCA (svd, Wehrens11, p48)
                   ##-----------------------------

                   pcaSvdLs <- svd(tcrossprod(xMN))

                   tMN <- pcaSvdLs[["u"]] %*% diag(sqrt(pcaSvdLs[["d"]]))

                   pMN <- t(solve(tMN, xMN))

                   varVn <- pcaSvdLs[["d"]] / (nrow(xMN) - 1)
                   ##  length   class    mode typeof  size
                   ##      50 numeric numeric double 2'304
                   ##  249.014 202.008 ... 14.658 0
                   ## Names:  t1 t2 ... t49 t50

                   rm(pcaSvdLs)

                   tMN <- tMN[, 1:predI, drop = FALSE]
                   pMN <- pMN[, 1:predI, drop = FALSE]
                   varVn <- varVn[1:predI]
                   rownames(tMN) <- obsNamVc
                   rownames(pMN) <- xvaNamVc
                   names(varVn) <- colnames(pMN) <- colnames(tMN) <- preNamVc

                   modelDF[, "R2X"] <- round(varVn / vSumVn, 3)

               }) ## svd

        modelDF[, "R2X(cum)"] <- cumsum(modelDF[, "R2X"])

        if(autNcpL) {
            vSelVl <- cumsum(varVn) / vSumVn > 0.5
            vSelVi <- which(vSelVl)
            if(sum(vSelVl) == 0) {
                warning("The maximum number of components for the automated mode (", autMaxN, ") has been reached whereas the cumulative variance ", round(tail(modelDF[, "R2X(cum)"], 1) * 100), "% is still less than 50%.", call. = FALSE)
            } else
                predI <- max(which(vSelVl)[1])

            tMN <- tMN[, 1:predI, drop = FALSE]
            pMN <- pMN[, 1:predI, drop = FALSE]
            varVn <- varVn[1:predI]
            modelDF <- modelDF[1:predI, , drop = FALSE]
        }

        summaryDF <- modelDF[predI, c("R2X(cum)"), drop = FALSE]

    } else { ## if(is.null(yMCN))


        ## Y-Scaling
        ##---------------

        yMeanVn <- apply(yMN, 2, function(colVn) mean(colVn, na.rm = TRUE))

        if(mode(yMCN) == "character") {

            ySdVn <- apply(yMN, 2, function(colVn) sd(colVn, na.rm = TRUE))

        } else {

            switch(scaleC,
                   none = {
                       yMeanVn <- rep(0, times = ncol(yMN))
                       ySdVn <- rep(1, times = ncol(yMN))
                   },
                   center = {
                       ySdVn <- rep(1, times = ncol(yMN))
                   },
                   pareto = {
                       ySdVn <- apply(yMN, 2, function(colVn) sqrt(sd(colVn, na.rm = TRUE)))
                   },
                   standard = {
                       ySdVn <- apply(yMN, 2, function(colVn) sd(colVn, na.rm = TRUE))
                   })

        }

        yMN <- scale(yMN, center = yMeanVn, scale = ySdVn)


        if(!is.null(colnames(yMN))) {
            yvaNamVc <- colnames(yMN)
        } else
            yvaNamVc <- paste("y", 1:ncol(yMN), sep = "")


        wMN <- pMN
        uMN <- tMN

        cMN <- matrix(0,
                      nrow = ncol(yMN),
                      ncol = predI,
                      dimnames = list(yvaNamVc, preNamVc))


        ## Cross-validation variables

        cvfNamVc <- paste("cv", 1:crossvalI, sep = "")
        cvfOutLs <- split(1:nrow(xMN), rep(1:crossvalI, length = nrow(xMN)))

        prkVn <- numeric(crossvalI) ## PRESS for each cv fold

        ## rules

        ru1ThrN <- ifelse(orthoI == 0,
                          ifelse(nrow(xMN) > 100, yes = 0, no = 0.05), ## PLS
                          0.01) ## OPLS

        ## SSY total

        ssyTotN <- rs0N <- sum(yMN^2, na.rm = TRUE)


        hN <- 1


        if(orthoI == 0) {

            ##------------------------------------
            ##   PLS
            ##------------------------------------


            xnMN <- xMN
            ynMN <- yMN
            ## 'n' stands for NIPALS (OPLS followed by PLS); original matrices are needed for cross-validation

            modelDF <- as.data.frame(matrix(NA,
                                          nrow = predI,
                                          ncol = 8,
                                          dimnames = list(preNamVc, c("R2X", "R2X(cum)", "R2Y", "R2Y(cum)", "Q2", "Q2(cum)", "Signif.", "Iter."))))
            for(j in 1:ncol(modelDF))
                mode(modelDF[, j]) <- ifelse(colnames(modelDF)[j] == "Signif.", "character", "numeric")

            rssN <- rs0N

            while(hN < (predI + 1)) {

                iteN <- 1

                uVn <- ynMN[, 1, drop = FALSE]
                tOldVn <- matrix(0, nrow = nrow(xMN))

                repeat {

                    if(naxL || nayL) {
                        wVn <- numeric(ncol(xnMN))
                        for(j in 1:ncol(xnMN)) {
                            comVl <- complete.cases(xnMN[, j]) &
                                complete.cases(uVn)
                            wVn[j] <- crossprod(xnMN[comVl, j], uVn[comVl]) / drop(crossprod(uVn[comVl]))
                        }
                    } else
                        wVn <- crossprod(xnMN, uVn) / drop(crossprod(uVn))

                    wVn <- wVn / sqrt(drop(crossprod(wVn)))

                    if(naxL) {
                        tVn <- numeric(nrow(xnMN))
                        for(i in 1:nrow(xnMN)) {
                            comVl <- complete.cases(xnMN[i, ])
                            tVn[i] <- crossprod(xnMN[i, comVl], wVn[comVl]) / drop(crossprod(wVn[comVl]))
                        }
                    } else
                        tVn <- xnMN %*% wVn

                    if(nayL) {
                        cVn <- numeric(ncol(ynMN))
                        for(j in 1:ncol(ynMN)) {
                            comVl <- complete.cases(ynMN[, j])
                            cVn[j] <- crossprod(ynMN[comVl, j], tVn[comVl]) / drop(crossprod(tVn[comVl]))
                        }
                    } else
                        cVn <- crossprod(ynMN, tVn) / drop(crossprod(tVn))


                    if(ncol(ynMN) == 1 ||
                       drop(sqrt(crossprod((tOldVn - tVn) / tVn))) < 1e-6) {

                        break

                    } else {

                        if(nayL) {
                            uVn <- numeric(nrow(xnMN))
                            for(i in 1:nrow(xnMN)) {
                                comVl <- complete.cases(ynMN[i, ])
                                uVn[i] <- crossprod(ynMN[i, comVl], cVn[comVl]) / drop(crossprod(cVn[comVl]))
                            }
                        } else
                            uVn <- ynMN %*% cVn / drop(crossprod(cVn))

                        tOldVn <- tVn

                        iteN <- iteN + 1

                    }

                }

                if(naxL) {
                    pVn <- numeric(ncol(xnMN))
                    for(j in 1:ncol(xnMN)) {
                        comVl <- complete.cases(xnMN[, j])
                        pVn[j] <- crossprod(xnMN[comVl, j], tVn[comVl]) / drop(crossprod(tVn[comVl]))
                    }
                } else
                    pVn <- crossprod(xnMN, tVn) / drop(crossprod(tVn))

                wMN[, hN] <- wVn
                tMN[, hN] <- tVn
                pMN[, hN] <- pVn

                cMN[, hN] <- cVn
                uMN[, hN] <- uVn

                if(naxL)
                    modelDF[hN, "R2X"] <- sum((tcrossprod(tMN[, hN], pMN[, hN])[!is.na(xMN)])^2) / ssxTotN
                else
                    modelDF[hN, "R2X"] <- sum(tcrossprod(tMN[, hN], pMN[, hN])^2) / ssxTotN

                if(nayL)
                    modelDF[hN, "R2Y"] <- sum((tcrossprod(tMN[, hN], cMN[, hN])[!is.na(yMN)])^2) / ssyTotN
                else
                    modelDF[hN, "R2Y"] <- sum(tcrossprod(tMN[, hN], cMN[, hN])^2) / ssyTotN
                modelDF[hN, "Iter."] <- iteN

                ## cross-validation (PRESS computation)

                for(k in 1:crossvalI) {

                    ckxMN <- xnMN[-cvfOutLs[[k]], , drop = FALSE]
                    ckyMN <- ynMN[-cvfOutLs[[k]], , drop = FALSE]

                    ckuVn <- ckyMN[, 1, drop = FALSE]
                    cktOldVn <- 0

                    nkxL <- any(is.na(c(ckxMN)))
                    nkyL <- any(is.na(c(ckyMN)))
                    nkuL <- any(is.na(ckuVn))


                    repeat {

                        if(nkxL || nkyL) {
                            ckwVn <- numeric(ncol(ckxMN))
                            for(j in 1:ncol(ckxMN)) {
                                comVl <- complete.cases(ckxMN[, j]) &
                                    complete.cases(ckuVn)
                                ## ckwVn[j] <- crossprod(ckxMN[comVl, j], ckuVn[comVl])
                                ckwVn[j] <- crossprod(ckxMN[comVl, j], ckuVn[comVl]) / drop(crossprod(ckuVn[comVl]))
                            }
                        } else
                            ckwVn <- drop(crossprod(ckxMN, ckuVn)) / drop(crossprod(ckuVn))

                        ckwVn <- ckwVn / sqrt(drop(crossprod(ckwVn)))

                        if(nkxL) {
                            cktVn <- numeric(nrow(ckxMN))
                            for(i in 1:nrow(ckxMN)) {
                                comVl <- complete.cases(ckxMN[i, ])
                                cktVn[i] <- crossprod(ckxMN[i, comVl], ckwVn[comVl]) / drop(crossprod(ckwVn[comVl]))
                                ## cktVn[i] <- crossprod(ckxMN[i, comVl], ckwVn[comVl])
                            }
                        } else
                            cktVn <- ckxMN %*% ckwVn

                        if(nkyL) {
                            ckcVn <- numeric(ncol(ckyMN))
                            for(j in 1:ncol(ckyMN)) {
                                comVl <- complete.cases(ckyMN[, j])
                                ckcVn[j] <- crossprod(ckyMN[comVl, j], cktVn[comVl]) / drop(crossprod(cktVn[comVl]))
                            }
                        } else
                            ckcVn <- crossprod(ckyMN, cktVn) / drop(crossprod(cktVn))

                        if(ncol(ckyMN) == 1 ||
                           drop(sqrt(crossprod((cktOldVn - cktVn) / cktVn))) < 1e-6) {

                            break

                        } else {

                            if(nkyL) {
                                ckuVn <- numeric(nrow(ckxMN))
                                for(i in 1:nrow(ckxMN)) {
                                    comVl <- complete.cases(ckyMN[i, ])
                                    ckuVn[i] <- crossprod(ckyMN[i, comVl], ckcVn[comVl]) / drop(crossprod(ckcVn[comVl]))
                                }
                            } else
                                ckuVn <- ckyMN %*% ckcVn / drop(crossprod(ckcVn))

                            cktOldVn <- cktVn

                        }

                    }

                    if(any(is.na(xnMN[cvfOutLs[[k]], ]))) {
                        prxVn <- numeric(length(cvfOutLs[[k]]))
                        for(r in 1:length(prxVn)) {
                            comVl <- complete.cases(xnMN[cvfOutLs[[k]][r], ])
                            ## prxVn[r] <- crossprod(xnMN[cvfOutLs[[k]][r], comVl], ckwVn[comVl])
                            prxVn[r] <- crossprod(xnMN[cvfOutLs[[k]][r], comVl], ckwVn[comVl]) / drop(crossprod(ckwVn[comVl]))
                        }
                        prkVn[k] <- sum((ynMN[cvfOutLs[[k]], , drop = FALSE] - prxVn %*% t(ckcVn))^2, na.rm = TRUE)
                    } else
                        prkVn[k] <- sum((ynMN[cvfOutLs[[k]], , drop = FALSE] - xnMN[cvfOutLs[[k]], , drop = FALSE] %*% ckwVn %*% t(ckcVn))^2, na.rm = TRUE)

                } ## for(k in 1:crossvalI) {

                prsN <- sum(prkVn)

                modelDF[hN, "Q2"] <- 1 - prsN / rssN

                if(modelDF[hN, "R2Y"] < 0.01) {
                    modelDF[hN, "Signif."] <- "N4"
                } else if(modelDF[hN, "Q2"] < ru1ThrN) {
                    modelDF[hN, "Signif."] <- "NS"
                } else
                    modelDF[hN, "Signif."] <- "R1"

                if(autNcpL && modelDF[hN, "Signif."] != "R1" && hN >= 1)
                    break

                rssN <- sum((ynMN - tcrossprod(tVn, cVn))^2, na.rm = TRUE)

                xnMN <- xnMN - tcrossprod(tVn, pVn)
                ynMN <- ynMN - tcrossprod(tVn, cVn)

                hN <- hN + 1

            } ## for(hN in 1:predI) {

            rm(ckxMN)
            rm(ckyMN)

            modelDF[, "R2X(cum)"] <- cumsum(modelDF[, "R2X"])
            modelDF[, "R2Y(cum)"] <- cumsum(modelDF[, "R2Y"])
            modelDF[, "Q2(cum)"] <- 1 - cumprod(1 - modelDF[, "Q2"])

            if(autNcpL) {

                hN <- hN - 1

                if(hN == 0)
                    stop("No model was built because the first predictive component was already not significant;\nSelect a number of predictive components of 1 if you want the algorithm to compute a model despite this.", call. = FALSE)

                if(hN == autMaxN)
                    warning("The maximum number of components in the automated mode (", autMaxN, ") has been reached whereas R2Y (", round(modelDF[hN, 'R2Y'] * 100), "%) is still above 1% and Q2Y (", round(modelDF[hN, 'Q2'] * 100), "%) is still above ", round(ru1ThrN * 100), "%.", call. = FALSE)

                wMN <- wMN[, 1:hN, drop = FALSE]
                tMN <- tMN[, 1:hN, drop = FALSE]
                pMN <- pMN[, 1:hN, drop = FALSE]
                cMN <- cMN[, 1:hN, drop = FALSE]
                uMN <- uMN[, 1:hN, drop = FALSE]

                preNamVc <- preNamVc[1:hN]

                predI <- hN

                modelDF <- modelDF[1:hN, , drop = FALSE]

            }

            summaryDF <- modelDF[predI, c("R2X(cum)", "R2Y(cum)", "Q2(cum)")]


            ## WeightStar matrix (W*)

            if(predI == 1) {

                rMN <- wMN

            } else {

                pwMN <- crossprod(pMN, wMN)
                rMN <- wMN %*% solve(pwMN)
                colnames(rMN) <- preNamVc

            }

            rm(xnMN)
            rm(ynMN)


            ## Regression coefficients

            bMN <- tcrossprod(rMN, cMN)


            ## Predicted values

            yPreScaMN <- tcrossprod(tMN, cMN)
            if(nayL && ncol(yMN) == 1)
                yPreScaMN <- yPreScaMN[!is.na(yMN), , drop = FALSE]

            yPreMN <- scale(scale(yPreScaMN,
                                   FALSE,
                                   1 / ySdVn),
                             -yMeanVn,
                             FALSE)
            attr(yPreMN, "scaled:center") <- NULL
            attr(yPreMN, "scaled:scale") <- NULL


            if(!is.null(subset)) {
                yActMCN <- yMCN[subset, , drop = FALSE]
            } else
                yActMCN <- yMCN

            if(mode(yMCN) == "character") {
                yActMN <- .char2numF(yActMCN)
                                     ## , c2nLs = c2nLs)
            } else
                yActMN <- yActMCN


            summaryDF[, "RMSEE"] <- sqrt(.errorF(yActMN, yPreMN)^2 * nrow(yActMN) / (nrow(yActMN) - (1 + predI))) ## for SIMCA compatibility


            if(!is.null(subset)) { ## tRaining/tEst partition

                xteMN <- scale(xTesMN, xMeanVn, xSdVn)

                if(naxL) {
                    yTesScaMN <- matrix(0, nrow = nrow(xteMN), ncol = ncol(bMN))
                    for(j in 1:ncol(yTesScaMN))
                        for(i in 1:nrow(yTesScaMN)) {
                            comVl <- complete.cases(xteMN[i, ])
                            yTesScaMN[i, j] <- crossprod(xteMN[i, comVl], bMN[comVl, j])
                        }
                } else
                    yTesScaMN <- xteMN %*% bMN

                if(nayL)
                    yTesScaMN <- yTesScaMN[setdiff(1:row(yTesScaMN), union(subset, nayVi)), , drop = FALSE]
                ## yTesScaMN <- yTesScaMN[!is.na(yMCN[testVi, ]), , drop = FALSE]

                yTesMN <- scale(scale(yTesScaMN,
                                       FALSE,
                                       1 / ySdVn),
                                 -yMeanVn,
                                 FALSE) ## predicted values
                attr(yTesMN, "scaled:center") <- NULL
                attr(yTesMN, "scaled:scale") <- NULL

                if(mode(yMCN) == "character") {
                    yTestMCN <- .char2numF(yTesMN,
                                           c2nL = FALSE)
                } else
                    yTestMCN <- yTesMN

                yTesActMCN <- yMCN[setdiff(1:nrow(yMCN), subset), , drop = FALSE] ## actual values
                if(mode(yMCN) == "character") {
                    yTesActMN <- .char2numF(yTesActMCN)
                } else
                    yTesActMN <- yTesActMCN

                summaryDF[, "RMSEP"] <- .errorF(c(yTesMN), c(yTesActMN))

            } else
                yTestMCN <- NULL


        } else { ## orthoI > 0


            ##------------------------------------
            ##   OPLS
            ##------------------------------------


            ## Trygg and Wold (2002).
            ## Orthogonal projections to latent structures (O-PLS).
            ## Journal of Chemometrics. 16:119-128.

            orthoNamVc <- paste("o", 1:orthoI, sep = "")

            toMN <- matrix(0,
                           nrow = nrow(xMN),
                           ncol = orthoI,
                           dimnames = list(obsNamVc, orthoNamVc))
            woMN <- poMN <- matrix(0,
                                   nrow = ncol(xMN),
                                   ncol = orthoI,
                                   dimnames = list(xvaNamVc, orthoNamVc))
            coMN <- matrix(0,
                           nrow = ncol(yMN),
                           ncol = orthoI,
                           dimnames = list(yvaNamVc, orthoNamVc))

            modelDF <- as.data.frame(matrix(NA,
                                            nrow = 1 + orthoI + 1,
                                            ncol = 7,
                                            dimnames = list(c("p1", orthoNamVc, "sum"), c("R2X", "R2X(cum)", "R2Y", "R2Y(cum)", "Q2", "Q2(cum)", "Signif."))))
            for(j in 1:ncol(modelDF))
                mode(modelDF[, j]) <- ifelse(colnames(modelDF)[j] == "Signif.", "character", "numeric")

            xcvTraLs <- lapply(cvfOutLs,
                               function(obsVi)
                               xMN[-obsVi, , drop = FALSE])

            xcvTesLs <- lapply(cvfOutLs,
                               function(obsVi)
                               xMN[obsVi, , drop = FALSE])

            ycvTraLs <- lapply(cvfOutLs,
                               function(obsVi)
                               yMN[-obsVi, , drop = FALSE])

            ycvTesLs <- lapply(cvfOutLs,
                               function(obsVi)
                               yMN[obsVi, , drop = FALSE])

            ## full dataset added as crossvalI + 1 item

            xcvTraLs <- c(xcvTraLs,
                          list(xMN))

            ycvTraLs <- c(ycvTraLs,
                          list(yMN))


            breL <- FALSE

            for(noN in 1:(orthoI + 1)) {

                if(breL)
                    break

                for(cvN in 1:length(xcvTraLs)) {
                    ## cvN < length(xcvTraLs): cross-validation
                    ## cvN == length(xcvTraLs): full dataset

                    xcvTraMN <- xcvTraLs[[cvN]]
                    ycvTraMN <- ycvTraLs[[cvN]]

                    if(ncol(ycvTraMN) > 1) {

                        ## step -|1 [case vector y (p121) | matrix Y (p127)]

                        if(naxL || nayL)
                            wwMN <- apply(ycvTraMN,
                                          2,
                                          function(colVn) {
                                              wwjVn <- numeric(ncol(xcvTraMN))
                                              for(j in 1:ncol(xcvTraMN)) {
                                                  comVl <- complete.cases(xcvTraMN[, j]) & complete.cases(colVn)
                                                  wwjVn[j] <- crossprod(xcvTraMN[comVl,j], colVn[comVl]) / drop(crossprod(colVn[comVl]))
                                              }
                                              wwjVn
                                          })
                        else
                            wwMN <- apply(ycvTraMN,
                                          2,
                                          function(colVn)
                                          crossprod(xcvTraMN, colVn) / drop(crossprod(colVn)))

                        ## step -|2

                        wwSvdLs <- svd(wwMN)
                        wwNcpVin <- which(wwSvdLs[["d"]]^2 > epsN * sum(wwSvdLs[["d"]]^2))

                        twMN <- wwSvdLs[["u"]][, wwNcpVin, drop = FALSE] %*% diag(wwSvdLs[["d"]][wwNcpVin], nrow = length(wwNcpVin))

                    }


                    ## step -|4

                    uOldVn <- ycvTraMN[, 1, drop = FALSE]

                    repeat {

                        ## step 1|5

                        if(naxL || nayL) {
                            wVn <- numeric(ncol(xcvTraMN))
                            for(j in 1:ncol(xcvTraMN)) {
                                comVl <- complete.cases(xcvTraMN[, j]) &
                                    complete.cases(uOldVn)
                                wVn[j] <- crossprod(xcvTraMN[comVl, j], uOldVn[comVl]) / drop(crossprod(uOldVn[comVl]))
                            }
                        } else
                            wVn <- crossprod(xcvTraMN, uOldVn) / drop(crossprod(uOldVn))

                        ## step 2|6

                        wVn <- wVn / sqrt(drop(crossprod(wVn)))

                        ## step 3|7

                        if(naxL) {
                            tVn <- numeric(nrow(xcvTraMN))
                            for(i in 1:nrow(xcvTraMN)) {
                                comVl <- complete.cases(xcvTraMN[i, ])
                                tVn[i] <- crossprod(xcvTraMN[i, comVl], wVn[comVl]) / drop(crossprod(wVn[comVl]))
                            }
                        } else
                            tVn <- xcvTraMN %*% wVn

                        ## step 4|8

                        if(nayL) {
                            cVn <- numeric(ncol(ycvTraMN))
                            for(j in 1:ncol(ycvTraMN)) {
                                comVl <- complete.cases(ycvTraMN[, j])
                                cVn[j] <- crossprod(ycvTraMN[comVl, j], tVn[comVl]) / drop(crossprod(tVn[comVl]))
                            }
                        } else
                            cVn <- crossprod(ycvTraMN, tVn) / drop(crossprod(tVn))

                        ## step 5|9

                        if(nayL) {
                            uVn <- numeric(nrow(xcvTraMN))
                            for(i in 1:nrow(xcvTraMN)) {
                                comVl <- complete.cases(ycvTraMN[i, ])
                                uVn[i] <- crossprod(ycvTraMN[i, comVl], cVn[comVl]) / drop(crossprod(cVn[comVl]))
                            }
                        } else
                            uVn <- ycvTraMN %*% cVn / drop(crossprod(cVn))

                        if(nayL) {
                            comVl <- complete.cases(uOldVn)
                            dscN <- drop(sqrt(crossprod((uVn[comVl] - uOldVn[comVl] / uVn[comVl]))))
                        } else
                            dscN <- drop(sqrt(crossprod((uVn - uOldVn) / uVn)))

                        if(ncol(ycvTraMN) == 1 ||
                           dscN < 1e-10) {

                            break

                        } else {

                            uOldVn <- uVn

                        }

                    } ## end of repeat

                    ## step 6|

                    if(naxL) {
                        pVn <- numeric(ncol(xcvTraMN))
                        for(j in 1:ncol(xcvTraMN)) {
                            comVl <- complete.cases(xcvTraMN[, j])
                            pVn[j] <- crossprod(xcvTraMN[comVl, j], tVn[comVl]) / drop(crossprod(tVn[comVl]))
                        }
                    } else
                        pVn <- crossprod(xcvTraMN, tVn) / drop(crossprod(tVn))

                    ## step 7|

                    if(ncol(ycvTraMN) > 1)
                        for(j in 1:ncol(twMN))
                            woVn <- pVn - drop(crossprod(twMN[, j, drop = FALSE], pVn)) / drop(crossprod(twMN[, j, drop = FALSE])) * twMN[, j, drop = FALSE]
                    else
                        woVn <- pVn - drop(crossprod(wVn, pVn)) / drop(crossprod(wVn)) * wVn

                    ## step 8|

                    woVn <- woVn / sqrt(drop(crossprod(woVn)))

                    ## step 9|

                    if(naxL) {
                        toVn <- numeric(nrow(xcvTraMN))
                        for(i in 1:nrow(xcvTraMN)) {
                            comVl <- complete.cases(xcvTraMN[i, ])
                            toVn[i] <- crossprod(xcvTraMN[i, comVl], woVn[comVl]) / drop(crossprod(woVn[comVl]))
                        }
                    } else
                        toVn <- xcvTraMN %*% woVn

                    if(nayL) {
                        coVn <- numeric(ncol(ycvTraMN))
                        for(j in 1:ncol(ycvTraMN)) {
                            comVl <- complete.cases(ycvTraMN[, j])
                            coVn[j] <- crossprod(ycvTraMN[comVl, j], toVn[comVl]) / drop(crossprod(toVn[comVl]))
                        }
                    } else
                        coVn <- crossprod(ycvTraMN, toVn) / drop(crossprod(toVn))

                    ## step 10|

                    if(naxL) {
                        poVn <- numeric(ncol(xcvTraMN))
                        for(j in 1:ncol(xcvTraMN)) {
                            comVl <- complete.cases(xcvTraMN[, j])
                            poVn[j] <- crossprod(xcvTraMN[comVl, j], toVn[comVl]) / drop(crossprod(toVn[comVl]))
                        }
                    } else
                        poVn <- crossprod(xcvTraMN, toVn) / drop(crossprod(toVn))

                    ## step 12|

                    if(cvN <= crossvalI) { ## cross-validation

                        xcvTesMN <- xcvTesLs[[cvN]]
                        ycvTesMN <- ycvTesLs[[cvN]]

                        if(any(is.na(xcvTesMN))) {
                            prxVn <- numeric(nrow(xcvTesMN))
                            for(r in 1:length(prxVn)) {
                                comVl <- complete.cases(xcvTesMN[r, ])
                                prxVn[r] <- crossprod(xcvTesMN[r, comVl], wVn[comVl]) / drop(crossprod(wVn[comVl]))
                            }
                            prkVn[cvN] <- sum((ycvTesMN - prxVn %*% t(cVn))^2, na.rm = TRUE)
                        } else
                            prkVn[cvN] <- sum((ycvTesMN - xcvTesMN %*% wVn %*% t(cVn))^2, na.rm = TRUE)

                        if(naxL) {
                            toTesVn <- numeric(nrow(xcvTesMN))
                            for(i in 1:nrow(xcvTesMN)) {
                                comVl <- complete.cases(xcvTesMN[i, ])
                                toTesVn[i] <- crossprod(xcvTesMN[i, comVl], woVn[comVl]) / drop(crossprod(woVn[comVl]))
                            }
                        } else
                            toTesVn <- xcvTesMN %*% woVn

                        xcvTesLs[[cvN]] <- xcvTesMN - tcrossprod(toTesVn, poVn)

                        if(cvN == crossvalI) {
                            q2N <- 1 - sum(prkVn) / rs0N
                            if(noN == 1) {
                                modelDF["p1", "Q2(cum)"] <- modelDF["p1", "Q2"] <- q2N
                            } else {
                                modelDF[noN, "Q2(cum)"] <- q2N - modelDF["p1", "Q2"]
                                modelDF[noN, "Q2"] <- q2N - sum(modelDF[1:(noN - 1), "Q2"], na.rm = TRUE)
                            }
                        }

                    } else { ## cvN == crossvalI + 1 (full matrix)

                        ## R2Xp computed later on (since they are updated)

                        ## R2Yp

                        if(nayL) {
                            r2yN <- sum((tcrossprod(tVn, cVn)[!is.na(yMN)])^2) / ssyTotN
                        } else
                            r2yN <- sum(tcrossprod(tVn, cVn)^2) / ssyTotN

                        if(noN == 1) {
                            modelDF["p1", "R2Y(cum)"] <- modelDF["p1", "R2Y"] <- r2yN
                        } else {
                            modelDF[noN, "R2Y(cum)"] <- r2yN - modelDF["p1", "R2Y"]
                            modelDF[noN, "R2Y"] <- r2yN - sum(modelDF[1:(noN - 1), "R2Y"], na.rm = TRUE)
                        }

                        if(noN <= orthoI) {

                            ## R2Xoi (R2Yoi is 0)

                            if(naxL) {
                                modelDF[paste0("o", noN), "R2X"] <- sum((tcrossprod(toVn, poVn)[!is.na(xMN)])^2) / ssxTotN
                            } else
                                modelDF[paste0("o", noN), "R2X"] <- sum(tcrossprod(toVn, poVn)^2) / ssxTotN

                            poMN[, noN] <- poVn
                            toMN[, noN] <- toVn
                            woMN[, noN] <- woVn
                            coMN[, noN] <- coVn

                        }

                        if(modelDF[noN, "R2Y"] < 0.01) {
                            modelDF[noN, "Signif."] <- "N4"
                        } else if(modelDF[noN, "Q2"] < ru1ThrN) {
                            modelDF[noN, "Signif."] <- "NS"
                        } else
                            modelDF[noN, "Signif."] <- "R1"

                        if(autNcoL && modelDF[noN, "Signif."] != "R1" && noN > 2) {
                            breL <- TRUE
                            break
                        } else {
                            cMN[, 1] <- cVn
                            pMN[, 1] <- pVn
                            tMN[, 1] <- tVn
                            uMN[, 1] <- uVn
                            wMN[, 1] <- wVn
                        }
                    }

                    if(breL)
                        break

                    ## step 11|

                    if(noN < orthoI + 1)
                        xcvTraLs[[cvN]] <- xcvTraMN - tcrossprod(toVn, poVn)

                } ## for(cvN in 1:length(xcvTraLs)) {

            } ## for(noN in 1:(orthoI + 1)) {

            rm(xcvTraLs)
            rm(xcvTesLs)
            rm(ycvTraLs)

            ## R2X

            if(naxL) {
                modelDF["p1", "R2X(cum)"] <- modelDF["p1", "R2X"] <- sum((tcrossprod(tMN, pMN)[!is.na(xMN)])^2) / ssxTotN
            } else
                modelDF["p1", "R2X(cum)"] <- modelDF["p1", "R2X"] <- sum(tcrossprod(tMN, pMN)^2) / ssxTotN

            modelDF[1:(1 + orthoI), "R2X(cum)"] <- cumsum(modelDF[1:(1 + orthoI), "R2X"])

            if(autNcoL) {

                if(all(modelDF[, "Signif."] == "R1", na.rm = TRUE)) {
                    orthoI <- noN - 1
                } else
                    orthoI <- noN - 3

                if(orthoI == 0)
                    stop("No model was built because the first orthogonal component was already not significant;\nSelect a number of orthogonal components of 1 if you want the algorithm to compute a model despite this.", call. = FALSE)

                if(orthoI == autMaxN - 1)
                    warning("The maximum number of orthogonal components in the automated mode (", autMaxN - 1, ") has been reached whereas R2Y (", round(modelDF[1 + orthoI, 'R2Y'] * 100), "%) is above 1% and Q2Y (", round(modelDF[1 + orthoI, 'Q2'] * 100), "%) is still above ", round(ru1ThrN * 100), "%.", call. = FALSE)

                poMN <- poMN[, 1:orthoI, drop = FALSE]
                toMN <- toMN[, 1:orthoI, drop = FALSE]
                woMN <- woMN[, 1:orthoI, drop = FALSE]
                coMN <- coMN[, 1:orthoI, drop = FALSE]

                orthoNamVc <- orthoNamVc[1:orthoI]
                modelDF <- modelDF[c(1:(orthoI + 1), nrow(modelDF)), ]

            }

            ## R2X

            modelDF["sum", "R2X(cum)"] <- modelDF[1 + orthoI, "R2X(cum)"]

            ## R2Y

            modelDF["sum", "R2Y(cum)"] <- sum(modelDF[, "R2Y"], na.rm = TRUE)

            ## Q2

            modelDF["sum", "Q2(cum)"] <- sum(modelDF[, "Q2"], na.rm = TRUE)

            summaryDF <- modelDF["sum", c("R2X(cum)", "R2Y(cum)", "Q2(cum)")]

            ## WeightStar matrix (W*)

            rMN <- wMN ## only 1 predictive component for OPLS


            ## Regression coefficients

            bMN <- tcrossprod(rMN, cMN)


            ## Predicted values

            yPreScaMN <- tcrossprod(tMN, cMN)

            yPreMN <- scale(scale(yPreScaMN,
                                    FALSE,
                                    1 / ySdVn),
                              -yMeanVn,
                              FALSE)
            attr(yPreMN, "scaled:center") <- NULL
            attr(yPreMN, "scaled:scale") <- NULL

            if(!is.null(subset)) {
                yActMCN <- yMCN[subset, , drop = FALSE]
            } else
                yActMCN <- yMCN

            if(mode(yMCN) == "character") {
                yActMN <- .char2numF(yActMCN)
            } else
                yActMN <- yActMCN

            summaryDF[, "RMSEE"] <- sqrt(.errorF(yActMN, yPreMN)^2 * nrow(yActMN) / (nrow(yActMN) - (1 + predI + orthoI)))


            if(!is.null(subset)) { ## tRaining/tEst partition

                xteMN <- scale(xTesMN, xMeanVn, xSdVn)

                for(noN in 1:orthoI) {
                    if(naxL) {
                        xtoMN <- matrix(0, nrow = nrow(xteMN), ncol = 1)
                        for(i in 1:nrow(xtoMN)) {
                            comVl <- complete.cases(xteMN[i, ])
                            xtoMN[i, ] <- crossprod(xteMN[i, comVl], woMN[comVl, noN]) / drop(crossprod(woMN[comVl, noN]))
                        }
                    } else
                        xtoMN <- xteMN %*% woMN[, noN]

                    xteMN <- xteMN - tcrossprod(xtoMN, poMN[, noN])
                }

                if(naxL) {
                    yTesScaMN <- matrix(0, nrow = nrow(xteMN), ncol = ncol(bMN), dimnames = list(rownames(xteMN), colnames(bMN)))
                    for(j in 1:ncol(yTesScaMN))
                        for(i in 1:nrow(yTesScaMN)) {
                            comVl <- complete.cases(xteMN[i, ])
                            yTesScaMN[i, j] <- crossprod(xteMN[i, comVl], bMN[comVl, j])
                        }
                } else
                    yTesScaMN <- xteMN %*% bMN

                if(nayL)
                    yTesScaMN <- yTesScaMN[!is.na(yMCN[setdiff(1:nrow(yMCN), subset), ]), , drop = FALSE]

                yTesMN <- scale(scale(yTesScaMN,
                                      FALSE,
                                      1 / ySdVn),
                                -yMeanVn,
                                FALSE)
                attr(yTesMN, "scaled:center") <- NULL
                attr(yTesMN, "scaled:scale") <- NULL

                if(mode(yMCN) == "character") {
                    yTestMCN <- .char2numF(yTesMN,
                                           c2nL = FALSE)
                } else
                    yTestMCN <- yTesMN

                yTesActMCN <- yMCN[setdiff(1:nrow(yMCN), subset), , drop = FALSE] ## actual values
                if(mode(yMCN) == "character") {
                    yTesActMN <- .char2numF(yTesActMCN)
                } else
                    yTesActMN <- yTesActMCN

                summaryDF[, "RMSEP"] <- .errorF(c(yTesMN), c(yTesActMN))

            } else
                yTestMCN <- NULL


        } ## end of OPLS


        ## VIP (specific implementation required for OPLS(-DA))

        if(orthoI == 0) { ## sum(vipVn^2) == nrow(wMN) [number of features]

            ssyVn <-  sapply(1:ncol(tMN),
                             function(j) sum(drop(tcrossprod(tMN[, j], cMN[, j])^2)))

            vipVn <- sqrt(nrow(wMN) * rowSums(sweep(wMN^2,
                                                    2,
                                                    ssyVn,
                                                    "*")) / sum(ssyVn))

        } else {

            sxpVn <- sapply(1:ncol(tMN),
                            function(h)
                            sum(drop(tcrossprod(tMN[, h], pMN[, h])^2)))
            sxpCumN <- sum(sxpVn)
            sxoVn <- sapply(1:ncol(toMN),
                            function(h)
                            sum(drop(tcrossprod(toMN[, h], poMN[, h])^2)))
            sxoCumN <- sum(sxoVn)
            ssxCumN <- sxpCumN + sxoCumN

            sypVn <- sapply(1:ncol(tMN),
                            function(h)
                            sum(drop(tcrossprod(tMN[, h], cMN[, h])^2)))
            sypCumN <- sum(sypVn)
            syoVn <- sapply(1:ncol(toMN),
                            function(h)
                            sum(drop(tcrossprod(toMN[, h], coMN[, h])^2)))
            syoCumN <- sum(syoVn)
            ssyCumN <- sypCumN + syoCumN

            ## VIP4,p [sum(vipVn^2) == nrow(wMN) instead of nrow(wMN) / 2 in the formula (but not in the figure) of the paper]

            kpN <- nrow(wMN) / (sxpCumN / ssxCumN + sypCumN / ssyCumN)

            pNorMN <- sweep(pMN, 2, sqrt(colSums(pMN^2)), "/") ## normalized loadings

            vipVn <- sqrt(kpN * (rowSums(sweep(pNorMN^2, 2, sxpVn, "*")) / ssxCumN + rowSums(sweep(pNorMN^2, 2, sypVn, "*")) / ssyCumN))

            ## VIP4,o [sum(orthoVipVn^2) == nrow(wMN) instead of nrow(wMN) / 2 in the formula (but not in the figure) of the paper]

            koN <- nrow(wMN) / (sxoCumN / ssxCumN + syoCumN / ssyCumN)

            poNorMN <- sweep(poMN, 2, sqrt(colSums(poMN^2)), "/")

            orthoVipVn <- sqrt(koN * (rowSums(sweep(poNorMN^2, 2, sxoVn, "*")) / ssxCumN + rowSums(sweep(poNorMN^2, 2, syoVn, "*")) / ssyCumN))


        }

    }

    summaryDF[, "pre"] <- predI
    summaryDF[, "ort"] <- orthoI
    rownames(summaryDF) <- "Total"

    sigNamVc <- c("R2X", "R2X(cum)", "R2Y", "R2Y(cum)", "Q2", "Q2(cum)", "RMSEE", "RMSEP")
    for(namC in intersect(colnames(modelDF), sigNamVc))
        modelDF[, namC] <- signif(modelDF[, namC], 3)
    for(namC in intersect(colnames(summaryDF), sigNamVc))
        summaryDF[, namC] <- signif(summaryDF[, namC], 3)

    ##------------------------------------
    ##   Returning
    ##------------------------------------


    retLs <- list(typeC = NULL,
                  descriptionMC = NULL,
                  modelDF = modelDF,
                  summaryDF = summaryDF,
                  subset = subset,

                  pcaVarVn = varVn,
                  vipVn = vipVn,
                  orthoVipVn = orthoVipVn,
                  fitted = NULL,
                  tested = NULL,
                  coefficients = bMN,
                  residuals = NULL,

                  xMeanVn = xMeanVn,
                  xSdVn = xSdVn,
                  yMeanVn = yMeanVn,
                  ySdVn = ySdVn,
                  xZeroVarVi = NULL,

                  scoreMN = tMN,
                  loadingMN = pMN,
                  weightMN = wMN,
                  orthoScoreMN = toMN,
                  orthoLoadingMN = poMN,
                  orthoWeightMN = woMN,
                  cMN = cMN,
                  uMN = uMN,
                  weightStarMN = rMN,
                  coMN = coMN,

                  suppLs = list(.char2numF = .char2numF,
                      yLevelVc = NULL,
                      algoC = algoC,
                      naxL = naxL,
                      nayL = nayL,
                      nayVi = nayVi,
                      permMN = NULL,
                      scaleC = scaleC,
                      topLoadI = NULL,
                      yMCN = yMCN,
                      xSubIncVarMN = NULL,
                      xCorMN = NULL,
                      xModelMN = xMN,
                      yModelMN = yMN,
                      yPreMN = yPreMN,
                      yTesMN = yTesMN))

} ## .coreF


.errorF <- function(x, y)
    sqrt(mean(drop((x - y)^2), na.rm = TRUE))


.log10F <- function(inpMN) {

    if(length(which(inpMN < 0)) > 0)
        stop("Negative values in the table to be log10 transformed", call. = FALSE)

    zerMN <- inpMN == 0

    inpMN[zerMN] <- 1

    return(log10(inpMN))

} ## .log10F

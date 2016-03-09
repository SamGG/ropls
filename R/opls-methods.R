setGeneric("opls", function(x, ...) standardGeneric("opls"))

setMethod("opls", signature(x = "data.frame"),
          function(x, ...) {
              if(!all(sapply(x, data.class) == "numeric")) {
                  stop("'x' data frame must contain columns of 'numeric' vectors only", call. = FALSE)
              } else
                  x <- as.matrix(x)
              opl <- opls(x, ...)
              opl
          })

setMethod("opls", signature(x = "matrix"),
          function(x,
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
        subsetL <- TRUE
        if(length(subset) == 1 && subset == "odd") {

            if(mode(yMCN) == "numeric")
                subsetVi <- seq(1, nrow(xMN), by = 2)
            else {
                subsetVi <- integer()
                for(claC in unique(drop(yMCN)))
                    subsetVi <- c(subsetVi,
                                  which(drop(yMCN) == claC)[seq(1, sum(drop(yMCN) == claC), by = 2)])
                subsetVi <- sort(subsetVi)
            }
        } else
            subsetVi <- subset
        if(crossvalI > length(subsetVi))
            stop("'crossvalI' must be less than the number of samples in the subset", call. = FALSE)
    } else {
        subsetL <- FALSE
        subsetVi <- numeric()
    }


    ## Filtering out zero variance variables
    ##--------------------------------------

    xVarIndLs <- list()
    xVarIndLs[[1]] <- 1:nrow(xMN)

    if(subsetL) {
        xVarIndLs[[1]] <- subsetVi
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

    opl <- .coreF(xMN = xMN,
                  yMCN = yMCN,
                  orthoI = orthoI,
                  predI = predI,
                  scaleC = scaleC,
                  algoC = algoC,
                  crossvalI = crossvalI,
                  subsetL = subsetL,
                  subsetVi = subsetVi,
                  .char2numF = .char2numF)

    opl@suppLs[["y"]] <- y

    if(is.null(opl@suppLs[["yMCN"]])) {
        opl@typeC <- "PCA"
    } else {
        if(ncol(opl@suppLs[["yMCN"]]) > 1 || mode(opl@suppLs[["yMCN"]]) == "numeric")
            opl@typeC <- "PLS"
        else
            opl@typeC <- "PLS-DA"
    }
    if(opl@summaryDF[, "ort"] > 0)
        opl@typeC <- paste("O", opl@typeC, sep = "")

    opl@xZeroVarVi <- xZeroVarVi
    ## opl@suppLs[["yLevelVc"]] <- yLevelVc


    ## Permutation testing (Szymanska et al, 2012)

    if(permI > 0) {

        modSumVc <- colnames(opl@summaryDF)

        permMN <- matrix(0,
                         nrow = 1 + permI,
                         ncol = length(modSumVc),
                         dimnames = list(NULL, modSumVc))

        perSimVn <- numeric(1 + permI)
        perSimVn[1] <- 1


        permMN[1, ] <- as.matrix(opl@summaryDF)

        for(k in 1:permI) {

            yVcn <- drop(opl@suppLs[["yMCN"]])
            if(!subsetL) {
                yPerVcn <- sample(yVcn)
            } else {
                yPerVcn <- numeric(nrow(xMN))
                refVi <- opl@subsetVi
                tesVi <- setdiff(1:nrow(xMN), refVi)
                yPerVcn[refVi] <- sample(yVcn[refVi])
                yPerVcn[tesVi] <- yVcn[tesVi]
            }
            yPerMCN <- matrix(yPerVcn, ncol = 1)

            perOpl <- .coreF(xMN = xMN,
                             yMCN = yPerMCN,
                             orthoI = opl@summaryDF[, "ort"],
                             predI = opl@summaryDF[, "pre"],
                             scaleC = scaleC,
                             algoC = algoC,
                             crossvalI = crossvalI,
                             subsetL = subsetL,
                             subsetVi = opl@subsetVi,
                             .char2numF = .char2numF)

            permMN[1 + k, ] <- as.matrix(perOpl@summaryDF)

            perSimVn[1 + k] <- .similarityF(opl@suppLs[["yMCN"]], yPerMCN,
                                            .char2numF = .char2numF,
                                            charL = mode(opl@suppLs[["yMCN"]]) == "character")

        }

        permMN <- cbind(permMN, sim = perSimVn)

        perPvaVn <- c(pR2Y = (1 + length(which(permMN[-1, "R2Y(cum)"] >= permMN[1, "R2Y(cum)"]))) / (nrow(permMN) - 1),
                      pQ2 = (1 + length(which(permMN[-1, "Q2(cum)"] >= permMN[1, "Q2(cum)"]))) / (nrow(permMN) - 1))
        opl@summaryDF[, "pR2Y"] <- perPvaVn["pR2Y"]
        opl@summaryDF[, "pQ2"] <- perPvaVn["pQ2"]

        opl@suppLs[["permMN"]] <- permMN

    }

    ##------------------------------------
    ##   Numerical results
    ##------------------------------------

    opl@descriptionMC <- rbind(samples = ifelse(!subsetL,
                                   nrow(xMN),
                                   length(subsetVi)),
                               X_variables = ncol(xMN),
                               near_zero_excluded_X_variables = length(opl@xZeroVarVi))

    totN <- length(c(xMN))
    nasN <- sum(is.na(c(xMN)))

    if(!is.null(opl@suppLs[["yMCN"]])) {

        opl@descriptionMC <- rbind(opl@descriptionMC,
                                         Y_variables = ncol(opl@suppLs[["yMCN"]]))
        totN <- totN + length(c(opl@suppLs[["yMCN"]]))
        nasN <- nasN + sum(is.na(c(opl@suppLs[["yMCN"]])))

    }

    opl@descriptionMC <- rbind(opl@descriptionMC,
                               missing_values = paste0(nasN, " (", round(nasN / totN * 100), "%)"))


    ## Raw summary
    ##------------

    opl@suppLs[["topLoadI"]] <- 3

    if(ncol(xMN) > opl@suppLs[["topLoadI"]]) {
        xVarVn <- apply(xMN, 2, var)
        names(xVarVn) <- 1:length(xVarVn)
        xVarVn <- sort(xVarVn)
        xVarSorVin <- as.numeric(names(xVarVn[seq(1, length(xVarVn), length = opl@suppLs[["topLoadI"]])]))
        opl@suppLs[["xSubIncVarMN"]] <- xMN[, xVarSorVin, drop = FALSE]
    } else
        opl@suppLs[["xSubIncVarMN"]] <- xMN

    if(ncol(xMN) <= 100) {

        xCorMN <- cor(xMN, use = "pairwise.complete.obs")
        xCorMN[lower.tri(xCorMN, diag = TRUE)] <- 0

        if(ncol(xMN) > opl@suppLs[["topLoadI"]]) {

            xCorNexDF <- which(abs(xCorMN) >= sort(abs(xCorMN), decreasing = TRUE)[opl@suppLs[["topLoadI"]] + 1],
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

        opl@suppLs[["xCorMN"]] <- xCorDisMN

        rm(xCorDisMN)

    }

    ## Printing
    ##---------

    if(printL) {
        show(opl)
        warnings()
    }

    ## Plotting
    ##---------

    if(plotL)
        plot(opl, typeVc = "summary")

    ## Closing connection
    ##-------------------

    if(!is.null(.sinkC)) ## Used in the Galaxy module
        sink()

    ## Returning
    ##----------

    return(invisible(opl))

}) ## end of opls.default


## Core algorithms for PCA, PLS(-DA), and OPLS(-DA)
.coreF <- function(xMN,
                   yMCN,
                   orthoI,
                   predI,
                   scaleC,
                   algoC,
                   crossvalI,
                   subsetL,
                   subsetVi,
                   .char2numF = .char2numF){

    epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16

    ##------------------------------------
    ##   Initialization
    ##------------------------------------

    varVn <- numeric()      ## PCA only
    yMeanVn <- numeric()    ## (O)PLS only
    ySdVn <- numeric()      ## (O)PLS only
    wMN <- matrix()         ## (O)PLS only
    cMN <- matrix()         ## (O)PLS only
    uMN <- matrix()         ## (O)PLS only
    rMN <- matrix()         ## (O)PLS only
    bMN <- matrix()         ## (O)PLS only
    vipVn <- numeric()      ## (O)PLS only
    yPreMN <- matrix()      ## (O)PLS only
    yTesMN <- matrix()      ## (O)PLS only
    toMN <- matrix()        ## OPLS only
    poMN <- matrix()        ## OPLS only
    woMN <- matrix()        ## OPLS only
    coMN <- matrix()        ## OPLS only
    orthoVipVn <- numeric() ## OPLS only
    mode(wMN) <- mode(cMN) <- mode(uMN) <- mode(rMN) <- mode(bMN) <- mode(yPreMN) <- mode(yTesMN) <- "numeric"
    mode(toMN) <- mode(poMN) <- mode(woMN) <- mode(coMN) <- "numeric"


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

        if(subsetL) {

            xTesMN <- xMN[-subsetVi, , drop = FALSE]
            xMN <- xMN[subsetVi, , drop = FALSE]
            yMN <- yMN[subsetVi, , drop = FALSE]

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


            if(subsetL) {
                yActMCN <- yMCN[subsetVi, , drop = FALSE]
            } else
                yActMCN <- yMCN

            if(mode(yMCN) == "character") {
                yActMN <- .char2numF(yActMCN)
                                     ## , c2nLs = c2nLs)
            } else
                yActMN <- yActMCN


            summaryDF[, "RMSEE"] <- sqrt(.errorF(yActMN, yPreMN)^2 * nrow(yActMN) / (nrow(yActMN) - (1 + predI))) ## for SIMCA compatibility


            if(subsetL) { ## tRaining/tEst partition

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
                    yTesScaMN <- yTesScaMN[setdiff(1:row(yTesScaMN), union(subsetVi, nayVi)), , drop = FALSE]
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

                yTesActMCN <- yMCN[setdiff(1:nrow(yMCN), subsetVi), , drop = FALSE] ## actual values
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

            if(subsetL) {
                yActMCN <- yMCN[subsetVi, , drop = FALSE]
            } else
                yActMCN <- yMCN

            if(mode(yMCN) == "character") {
                yActMN <- .char2numF(yActMCN)
            } else
                yActMN <- yActMCN

            summaryDF[, "RMSEE"] <- sqrt(.errorF(yActMN, yPreMN)^2 * nrow(yActMN) / (nrow(yActMN) - (1 + predI + orthoI)))


            if(subsetL) { ## tRaining/tEst partition

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
                    yTesScaMN <- yTesScaMN[!is.na(yMCN[setdiff(1:nrow(yMCN), subsetVi), ]), , drop = FALSE]

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

                yTesActMCN <- yMCN[setdiff(1:nrow(yMCN), subsetVi), , drop = FALSE] ## actual values
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

    opl <- new("opls")
    opl@typeC <- character()
    opl@descriptionMC <- matrix()
    opl@modelDF <- modelDF
    opl@summaryDF <- summaryDF
    opl@subsetVi <- subsetVi

    opl@pcaVarVn <- varVn
    opl@vipVn <- vipVn
    opl@orthoVipVn <- orthoVipVn
    opl@coefficientMN <- bMN

    opl@xMeanVn <- xMeanVn
    opl@xSdVn <- xSdVn
    opl@yMeanVn <- yMeanVn
    opl@ySdVn <- ySdVn
    opl@xZeroVarVi <- numeric()

    opl@scoreMN <- tMN
    opl@loadingMN <- pMN
    opl@weightMN <- wMN
    opl@orthoScoreMN <- toMN
    opl@orthoLoadingMN <- poMN
    opl@orthoWeightMN <- woMN
    opl@cMN <- cMN
    opl@uMN <- uMN
    opl@weightStarMN <- rMN
    opl@coMN <- coMN

    opl@suppLs <- list(.char2numF = .char2numF,
                       ## yLevelVc = NULL,
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
                       y = NULL,
                       xModelMN = xMN,
                       yModelMN = yMN,
                       yPreMN = yPreMN,
                       yTesMN = yTesMN)

    return(opl)


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


setMethod("show", "opls",
          function(object) {

    cat(object@typeC, "\n", sep = "")

    cat(object@descriptionMC["samples", ],
        " samples x ",
        object@descriptionMC["X_variables", ],
        " variables",
        ifelse(grepl("PLS", object@typeC),
               paste0(" and ", ncol(object@suppLs[["yMCN"]]), " response", ifelse(ncol(object@suppLs[["yMCN"]]) > 1, "s", "")),
               ""), "\n", sep = "")

    cat(object@suppLs[["scaleC"]], " scaling of predictors",
        ifelse(object@typeC == "PCA",
               "",
               paste0(" and ",
                      ifelse(mode(object@suppLs[["yMCN"]]) == "character" && object@suppLs[["scaleC"]] != "standard",
                             "standard scaling of ",
                             ""),
                      "response(s)")), "\n", sep = "")

    if(substr(object@descriptionMC["missing_values", ], 1, 1) != "0")
        cat(object@descriptionMC["missing_values", ], " NAs\n", sep = "")

    if(substr(object@descriptionMC["near_zero_excluded_X_variables", ], 1, 1) != "0")
        cat(object@descriptionMC["near_zero_excluded_X_variables", ],
            " excluded variables (near zero variance)\n", sep = "")

    optDigN <- options()[["digits"]]
    options(digits = 3)
    print(object@summaryDF)
    options(digits = optDigN)

}) ## show


setMethod("print", "opls",
          function(x, ...) {

    cat("\n1) Data set:\n", sep = "")

    cat(x@descriptionMC["samples", ],
            " samples x ",
            x@descriptionMC["X_variables", ],
            " variables",
            ifelse(grepl("PLS", x@typeC),
                   paste0(" and ", ncol(x@suppLs[["yMCN"]]), " response", ifelse(ncol(x@suppLs[["yMCN"]]) > 1, "s", "")),
                   ""),
        "\n", sep = "")

    cat(x@descriptionMC["missing_values", ], " NAs\n", sep = "")

    cat(x@descriptionMC["near_zero_excluded_X_variables", ], " excluded variables (near zero variance)\n", sep = "")

    cat(x@suppLs[["scaleC"]], " x", ifelse(x@typeC != "PCA", " and y", ""), " scaling\n", sep = "")

    cat("Summary of the ", x@suppLs[["topLoadI"]], " increasing variance spaced raw variables:\n", sep = "")
    print(summary(x@suppLs[["xSubIncVarMN"]]))

    if(!is.null(x@suppLs[["xCorMN"]])) {
        cat("Correlations between the X-variables:\n")
        print(signif(x@suppLs[["xCorMN"]], 2))
        cat("\n", sep = "")
    }

    cat("\n2) Model: ", x@typeC, "\n", sep = "")

    cat("Correlations between variables and first 2 components:\n", sep = "")

    if(x@summaryDF[, "pre"] + x@summaryDF[, "ort"] < 2) {

        warning("A single component model has been selected by cross-validation", call. = FALSE)

        tCompMN <- x@scoreMN
        pCompMN <- x@loadingMN

    } else {

        if(x@summaryDF[, "ort"] > 0) {
            tCompMN <- cbind(x@scoreMN[, 1], x@orthoScoreMN[, 1])
            pCompMN <- cbind(x@loadingMN[, 1], x@orthoLoadingMN[, 1])
            colnames(pCompMN) <- colnames(tCompMN) <- c("h1", "o1")
        } else {
            tCompMN <- x@scoreMN[, 1:2]
            pCompMN <- x@loadingMN[, 1:2]
        }

    }

    cxtCompMN <- cor(x@suppLs[["xModelMN"]], tCompMN,
                     use = "pairwise.complete.obs")

    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]])) {

        pexVi <- integer(x@suppLs[["topLoadI"]] * ncol(pCompMN) * 2) ## 'ex'treme values

        for(k in 1:ncol(pCompMN)) {

            pkVn <-  pCompMN[, k]

            pexVi[1:(2 * x@suppLs[["topLoadI"]]) + 2 * x@suppLs[["topLoadI"]] * (k - 1)] <- c(order(pkVn)[1:x@suppLs[["topLoadI"]]],
                                                                                rev(order(pkVn, decreasing = TRUE)[1:x@suppLs[["topLoadI"]]]))

        }

    } else
        pexVi <- 1:ncol(x@suppLs[["xModelMN"]])

    pxtCompMN <- cbind(pCompMN,
                       cxtCompMN)

    if(ncol(pCompMN) == 1) {
        colnames(pxtCompMN)[2] <- paste0("cor_", colnames(pxtCompMN)[2])
    } else
        colnames(pxtCompMN)[3:4] <- paste0("cor_", colnames(pxtCompMN)[3:4])

    topLoadMN <- pxtCompMN

    topLoadMN <- topLoadMN[pexVi, , drop = FALSE]

    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]]) &&
       ncol(pCompMN) > 1) {

        topLoadMN[(2 * x@suppLs[["topLoadI"]] + 1):(4 * x@suppLs[["topLoadI"]]), c(1, 3)] <- NA
        topLoadMN[1:(2 * x@suppLs[["topLoadI"]]), c(2, 4)] <- NA

    }
    print(signif(topLoadMN, 2))

    message("")

    optDigN <- options()[["digits"]]
    options(digits = 3)
    print(x@modelDF)
    options(digits = optDigN)

}) ## print


setMethod("plot", signature(x = "opls"),
          function(x,
                   y,
                   typeVc = c("correlation",
                       "outlier",
                       "overview",
                       "permutation",
                       "predict-train",
                       "predict-test",
                       "summary",
                       "x-loading",
                       "x-score",
                       "x-variance",
                       "xy-score",
                       "xy-weight")[7],
                   parAsColFcVn = NA,
                   parCexN = 0.8,
                   parCompVi = c(1, 2),
                   parDevNewL = TRUE,
                   parEllipsesL = NA,
                   parLabVc = NA,
                   parTitleL = TRUE,
                   file.pdfC = NULL,

                   .sinkC = NULL,
                   ...) {


    if(!is.null(.sinkC)) ##  Diversion of messages is required for the integration into Galaxy
        sink(.sinkC, append = TRUE)


    if("summary" %in% typeVc) {
        if(!is.null(x@suppLs[["permMN"]]))
            typeVc <- c("permutation",
                        "overview",
                        "outlier",
                        "x-score")
        else
            typeVc <- c("overview",
                        "outlier",
                        "x-score",
                        "x-loading")
    }


    ## Checking arguments
    ##-------------------

    if(!all(typeVc %in% c('correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight')))
        stop("'typeVc' elements must be either 'correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight'", call. = FALSE)

    if('predict-test' %in% typeVc && length(x@subsetVi) == 0)
        stop("For the 'predict-test' graphic to be generated, 'subset' must not be kept to NULL", call. = FALSE)

    if(!any(is.na(parLabVc))) {
        if(length(x@subsetVi) > 0 && length(parLabVc) != nrow(x@suppLs[["yMCN"]])) {
            stop("When 'subset' is not NULL, 'parLabVc' vector length must be equal to the number of train + test samples (here: ", nrow(x@suppLs[["yMCN"]]), ").", call. = FALSE)
        } else if(length(parLabVc) != nrow(x@scoreMN))
            stop("'parLabVc' vector length must be equal to the number of 'x' rows")
        if(mode(parLabVc) != "character")
            stop("'parLabVc' must be of 'character' type")
    }

    if(!any(is.na(parAsColFcVn))) {
        if(length(x@subsetVi) > 0 && length(parAsColFcVn) != nrow(x@suppLs[["yMCN"]])) {
            stop("When 'subset' is not NULL, 'parAsColFcVn' vector length must be equal to the number of train + test samples (here: ", nrow(x@suppLs[["yMCN"]]), ").", call. = FALSE)
        } else if(length(parAsColFcVn) != nrow(x@scoreMN))
            stop("'parAsColFcVn' vector length must be equal to the number of 'x' rows")
        if(!(mode(parAsColFcVn) %in% c("character", "numeric")))
            stop("'parAsColFcVn' must be of 'character' or 'numeric' type")
        if(is.character(parAsColFcVn)) {
            parAsColFcVn <- factor(parAsColFcVn)
            warning("Character 'parAsColFcVn' set to a factor", call. = FALSE)
        }
    }

    if(is.null(x@suppLs[["permMN"]]) && 'permutation' %in% typeVc)
        stop("'permI' must be > 0 for 'permutation' graphic to be plotted", call. = FALSE)

    if(x@summaryDF[, "ort"] > 0)
        if(parCompVi[1] != 1) {
            parCompVi[1] <- 1
            warning("OPLS: first component to display ('parCompVi' first value) set to 1", call. = FALSE)
        }

    if("xy-weight" %in% typeVc &&
       substr(x@typeC, 1, 3) != "PLS")
       ## (is.null(yMCN) || is.na(x@summaryDF[, "ort"]) || x@summaryDF[, "ort"] > 0))
        stop("'xy-weight graphic can be displayed only for PLS(-DA) models", call. = FALSE)

    if(any(grepl('predict', typeVc)))
       if(is.null(x@suppLs[["yMCN"]]) ||
          ncol(x@suppLs[["yMCN"]]) > 1 ||
          (mode(x@suppLs[["yMCN"]]) == "character" && length(unique(drop(x@suppLs[["yMCN"]]))) > 2))
    ## if(any(grepl('predict', typeVc)) && is.matrix(x@fitted"]]) && ncol(x@fitted"]]) > 1)
        ## if(any(grepl('predict', typeVc)) && (is.null(yMCN) || ncol(yMCN) != 1))
           stop("'predict' graphics available for single response regression or binary classification only", call. = FALSE)

    if(is.na(parEllipsesL)) {
        if((all(is.na(parAsColFcVn)) && grepl("-DA$", x@typeC)) ||
           (!all(is.na(parAsColFcVn)) && is.factor(parAsColFcVn))) {
            parEllipsesL <- TRUE
        } else
            parEllipsesL <- FALSE
        ## if((x@typeC == "PCA" && !all(is.na(parAsColFcVn)) && is.factor(parAsColFcVn)) || ## PCA case
        ##    grepl("-DA$", x@typeC)) { ## (O)PLS-DA cases
        ##     parEllipsesL <- TRUE
        ## } else
        ##     parEllipsesL <- FALSE
    } else if(parEllipsesL && !grepl("-DA$", x@typeC) && (all(is.na(parAsColFcVn)) || !is.factor(parAsColFcVn)))
        stop("Ellipses can be plotted for PCA (or PLS regression) only if the 'parAsColFcVn' is a factor",
             call. = FALSE)

    if(x@summaryDF[, "pre"] + x@summaryDF[, "ort"] < 2) {

        if(!all(typeVc %in% c("permutation", "overview"))) {
            warning("Single component model: only 'overview' and 'permutation' (in case of single response (O)PLS(-DA)) plots available", call. = FALSE)
            typeVc <- "overview"
            if(!is.null(x@suppLs[["permMN"]]))
                typeVc <- c("permutation", typeVc)
        }

        tCompMN <- x@scoreMN
        pCompMN <- x@loadingMN

    } else {

        if(x@summaryDF[, "ort"] > 0) {
            if(parCompVi[2] > x@summaryDF[, "ort"] + 1)
                stop("Selected orthogonal component for plotting (ordinate) exceeds the total number of orthogonal components of the model", call. = FALSE)
            tCompMN <- cbind(x@scoreMN[, 1], x@orthoScoreMN[, parCompVi[2] - 1])
            pCompMN <- cbind(x@loadingMN[, 1], x@orthoLoadingMN[, parCompVi[2] - 1])
            colnames(pCompMN) <- colnames(tCompMN) <- c("h1", paste("o", parCompVi[2] - 1, sep = ""))
        } else {
            if(max(parCompVi) > x@summaryDF[, "pre"])
                stop("Selected component for plotting as ordinate exceeds the total number of predictive components of the model", call. = FALSE)
            tCompMN <- x@scoreMN[, parCompVi, drop = FALSE]
            pCompMN <- x@loadingMN[, parCompVi, drop = FALSE]
        }

    }

    ## if(ncol(tCompMN) > 1) {

    ##     mahInvCovMN <- solve(cov(tCompMN))

    ##     pcaResMN <- cbind(sdsVn = apply(tCompMN,
    ##                           1,
    ##                           function(x) sqrt(t(as.matrix(x)) %*% mahInvCovMN %*% as.matrix(x))),
    ##                       odsVn = apply(x@suppLs[["xModelMN"]] - tcrossprod(tCompMN, pCompMN),
    ##                           1,
    ##                           function(x) sqrt(drop(crossprod(x[complete.cases(x)])))))

    ## } else
    ##     pcaResMN <- NULL

    cxtCompMN <- cor(x@suppLs[["xModelMN"]], tCompMN,
                     use = "pairwise.complete.obs")

    if(!is.null(x@suppLs[["yModelMN"]]))
        cytCompMN <- cor(x@suppLs[["yModelMN"]], tCompMN, use = "pairwise.complete.obs")


    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]])) {

        pexVi <- integer(x@suppLs[["topLoadI"]] * ncol(pCompMN) * 2) ## 'ex'treme values

        for(k in 1:ncol(pCompMN)) {

            pkVn <-  pCompMN[, k]

            pexVi[1:(2 * x@suppLs[["topLoadI"]]) + 2 * x@suppLs[["topLoadI"]] * (k - 1)] <- c(order(pkVn)[1:x@suppLs[["topLoadI"]]],
                                                                         rev(order(pkVn, decreasing = TRUE)[1:x@suppLs[["topLoadI"]]]))

        }

    } else
        pexVi <- 1:ncol(x@suppLs[["xModelMN"]])

    pxtCompMN <- cbind(pCompMN,
                       cxtCompMN)

    if(ncol(pCompMN) == 1) {
       colnames(pxtCompMN)[2] <- paste0("cor_", colnames(pxtCompMN)[2])
    } else
        colnames(pxtCompMN)[3:4] <- paste0("cor_", colnames(pxtCompMN)[3:4])

    topLoadMN <- pxtCompMN

    topLoadMN <- topLoadMN[pexVi, , drop = FALSE]

    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]]) &&
       ncol(pCompMN) > 1) {

        topLoadMN[(2 * x@suppLs[["topLoadI"]] + 1):(4 * x@suppLs[["topLoadI"]]), c(1, 3)] <- NA
        topLoadMN[1:(2 * x@suppLs[["topLoadI"]]), c(2, 4)] <- NA

    }


    ## Observation and variable names and colors
    ##------------------------------------------

    ## obsLabVc

    if(!any(is.na(parLabVc))) {
        obsLabVc <- parLabVc
    } else if(!is.null(x@suppLs[["yMCN"]]) && ncol(x@suppLs[["yMCN"]]) == 1) { ## (O)PLS of single response
        obsLabVc <- rownames(x@suppLs[["yMCN"]])
    } else { ## PCA
        if(!is.null(rownames(tCompMN))) {
            obsLabVc <- rownames(tCompMN)
        } else
            obsLabVc <- as.character(1:nrow(tCompMN))
    }

    if(length(x@subsetVi) > 0) {
        ## (O)PLS(-DA) models of a single 'y' response
        tesLabVc <- obsLabVc[-x@subsetVi]
        obsLabVc <- obsLabVc[x@subsetVi]
    } else
        tesLabVc <- ""

    ## obsColVc

    if(!any(is.na(parAsColFcVn))) {
        obsColVc <- .colorF(as.vector(parAsColFcVn))[["colVc"]]
        obsLegVc <- as.vector(parAsColFcVn)
    } else if(!is.null(x@suppLs[["yMCN"]]) && ncol(x@suppLs[["yMCN"]]) == 1) { ## (O)PLS of single response
        obsColVc <- .colorF(c(x@suppLs[["yMCN"]]))[["colVc"]]
        obsLegVc <- c(x@suppLs[["yMCN"]])
    } else { ## PCA
        obsColVc <- rep("black", nrow(tCompMN))
        obsLegVc <- NULL
    }

    if(length(x@subsetVi) > 0) {
        ## (O)PLS(-DA) models of a single 'y' response
        tesColVc <- obsColVc[-x@subsetVi]
        obsColVc <- obsColVc[x@subsetVi]
        if(!is.null(obsLegVc)) {
            tesLegVc <- obsLegVc[-x@subsetVi]
            obsLegVc <- obsLegVc[x@subsetVi]
        }
    }


    ## Layout
    ##-------

    if(!parDevNewL && length(typeVc) != 1)
        stop("'typeVc' must be of length 1 when 'parDevNewL' is set to FALSE", call. = FALSE)

    if(parDevNewL) {
        layRowN <- ceiling(sqrt(length(typeVc)))
        if(is.null(file.pdfC))
            dev.new()
        else
            pdf(file.pdfC)
        layout(matrix(1:layRowN^2, byrow = TRUE, nrow = layRowN))
    }

    layL <- !parDevNewL || length(typeVc) > 1


    ## Par
    ##----

    if(layL) {
        marVn <- c(4.6, 4.1, 2.6, 1.6)
    } else
        marVn <- c(5.1, 4.1, 4.1, 2.1)

    par(font=2, font.axis=2, font.lab=2, lwd=2,
        mar=marVn,
        pch=18)


    ## Graph
    ##------

    for(ploC in typeVc)
        .plotF(ploC,
               opl = x,
               obsColVc = obsColVc,
               obsLabVc = obsLabVc,
               obsLegVc = obsLegVc,
               layL = layL,
               parCexN = parCexN,
               parEllipsesL = parEllipsesL,
               parTitleL = parTitleL,
               parCompVi = parCompVi,
               typeVc = typeVc,
               tCompMN = tCompMN,
               pCompMN = pCompMN,
               cxtCompMN = cxtCompMN,
               cytCompMN = cytCompMN,
               ## pcaResMN = pcaResMN,
               topLoadMN = topLoadMN,
               pexVi = pexVi,
               tesColVc = tesColVc,
               tesLabVc = tesLabVc,
               tesLegVc = tesLegVc)

    if(layL)
        par(font=1, font.axis=1, font.lab=1, lwd=1,
            mar=c(5.1, 4.1, 4.1, 2.1),
            pch=1)

    if(!is.null(file.pdfC))
        dev.off()


    ## Closing connection
    ##-------------------

    if(!is.null(.sinkC)) ## Used in the Galaxy module
        sink()


}) ## plot


.plotF <- function(ploC,
                   opl,
                   obsColVc,
                   obsLabVc,
                   obsLegVc,
                   layL,
                   parCexN,
                   parEllipsesL,
                   parTitleL,
                   parCompVi,
                   typeVc,
                   tCompMN,
                   pCompMN,
                   cxtCompMN,
                   cytCompMN,
                   topLoadMN,
                   pexVi,
                   tesColVc,
                   tesLabVc,
                   tesLegVc) {

    ploPclF <- function() {

        xLimVn <- NULL
        yLimVn <- NULL

        ploColVc <- "black"

        if(ploC == "correlation") {

            maiC <- "Variable correlations"

            xLabC <- paste("with t",
                           parCompVi[1],
                           sep = "")

            yLabC <- paste("with t",
                           parCompVi[2],
                           sep = "")

            if(opl@summaryDF[, "ort"] > 0)
                yLabC <- paste("with tOrtho",
                               parCompVi[2] - 1,
                               sep = "")

            yLimVn <- xLimVn <- c(-1, 1)

            ploMN <- cxtCompMN

            if(opl@typeC != "PCA")
                ploMN <- rbind(ploMN,
                               cytCompMN)

        } else if(substr(ploC, 1, 7) == "predict") {

            maiC <- paste0("Predicted vs Actual",
                           paste0(" (",
                                 unlist(strsplit(ploC, "-"))[2],
                                 ")"))
            xLabC <- "predicted"
            yLabC <- "actual"

            if(grepl("train", ploC)) {
                ploNamVc <- obsLabVc
                ploColVc <- obsColVc
            } else {
                ploNamVc <- tesLabVc
                ploColVc <- tesColVc
            }

            ypMN <- eval(parse(text = paste("opl@suppLs[['y", switch(unlist(strsplit(ploC, "-"))[2], train = "Pre", test = "Tes"), "MN']]", sep = ""))) ## predicted

            if(length(opl@subsetVi) == 0) {
                yaMCN <- opl@suppLs[["yMCN"]] ## actual
            } else {
                if(grepl("train", ploC))
                    yaMCN <- opl@suppLs[["yMCN"]][opl@subsetVi, , drop = FALSE]
                else
                    yaMCN <- opl@suppLs[["yMCN"]][-opl@subsetVi, , drop = FALSE]
            }

            if(mode(opl@suppLs[["yMCN"]]) == "character") { ## binary only
                ypMN <- ypMN[, 1, drop = FALSE]
                yaMN <- opl@suppLs[[".char2numF"]](yaMCN)[, 1, drop = FALSE]
            } else
                yaMN <- yaMCN

            ploMN <- cbind(ypMN,
                           yaMN) ## to be modified (when ncol(yPreMCN) > 1)

        } else if(ploC == "x-loading") {

            maiC <- "Loadings"

            xLabC <- paste0("p",
                            parCompVi[1],
                            " (",
                            round(opl@modelDF[parCompVi[1], "R2X"] * 100),
                            "%)")

            yLabC <- paste0("p",
                            parCompVi[2],
                            " (",
                            round(opl@modelDF[parCompVi[2], "R2X"] * 100),
                            "%)")

            ploMN <- pCompMN

            if(!is.null(opl@suppLs[["yMCN"]]) && opl@summaryDF[, "ort"] > 0)
                yLabC <- paste0("pOrtho",
                                parCompVi[2] - 1,
                                " (",
                                round(opl@modelDF[parCompVi[2] - 1, "R2X"] * 100),
                                "%)")


        } else if(ploC == "x-score") {

            maiC <- paste0("Scores (", opl@typeC, ")")

            xLabC <- paste0("t",
                            parCompVi[1],
                            " (",
                            round(opl@modelDF[parCompVi[1], "R2X"] * 100),
                            "%)")

            yLabC <- paste0("t",
                            parCompVi[2],
                            " (",
                            round(opl@modelDF[parCompVi[2], "R2X"] * 100),
                           "%)")

            ploMN <- tCompMN

            if(grepl("^OPLS", opl@typeC))
                yLabC <- paste0("to", parCompVi[2] - 1)

            xLimVn <- c(-1, 1) * max(sqrt(var(ploMN[, 1]) * hotFisN), max(abs(ploMN[, 1])))
            yLimVn <- c(-1, 1) *max(sqrt(var(ploMN[, 2]) * hotFisN), max(abs(ploMN[, 2])))

            ploColVc <- obsColVc

        } else if(ploC == "xy-score") {

            maiC <- "XY-Scores"
            xLabC <- paste("t", parCompVi[1], sep = "")
            yLabC <- paste("u/c", parCompVi[1], sep = "")

            ploMN <- cbind(opl@scoreMN[, parCompVi[1]], opl@uMN[, parCompVi[1]] / opl@cMN[parCompVi[1]])

            ploColVc <- obsColVc

        } else if(ploC == "xy-weight") {

            maiC <- "Weights"
            xLabC <- paste0("w*c", parCompVi[1])
            yLabC <- paste0("w*c", parCompVi[2])

            ploMN <- rbind(opl@weightStarMN[, parCompVi],
                           opl@cMN[, parCompVi])

            pchVn <- rep(17, times = nrow(ploMN))
            ploColVc <- rep("grey", times = nrow(ploMN))

            pchVn[(nrow(opl@weightStarMN) + 1):nrow(ploMN)] <- 15
            ploColVc[(nrow(opl@weightStarMN) + 1):nrow(ploMN)] <- "black"

        }


        if(is.null(xLimVn))
            xLimVn <- range(ploMN[, 1])
        if(is.null(yLimVn))
            yLimVn <- range(ploMN[, 2])

        plot(ploMN,
             main=ifelse(parTitleL, maiC, ""),
             type = "n",
             xlab = xLabC,
             ylab = yLabC,
             xlim = xLimVn,
             ylim = yLimVn)

        abline(v = axTicks(1),
               col = "grey")

        abline(h = axTicks(2),
               col = "grey")

        abline(v = 0)
        abline(h = 0)

        if(ploC == "correlation") {

            lines(cos(radVn),
                  sin(radVn))

            corPexVi <- pexVi
            corPchVn <- rep(18, nrow(cxtCompMN))
            corNamVc <- rownames(cxtCompMN)
            ## corPchVn <- rep(18, ncol(opl@suppLs[["xModelMN"]]))
            ## corNamVc <- colnames(opl@suppLs[["xModelMN"]])


            if(opl@typeC != "PCA") {
                corPexVi <- c(corPexVi, (nrow(cxtCompMN) + 1):nrow(ploMN))
                corPchVn <- c(corPchVn, rep(15, nrow(cytCompMN)))
                corNamVc <- c(corNamVc, rownames(cytCompMN))
            }

            points(ploMN,
                   pch = corPchVn)

            points(ploMN[corPexVi, ],
                   pch = corPchVn[corPexVi],
                   col = "red")

            text(ploMN[corPexVi, ],
                 cex = parCexN,
                 col = "red",
                 labels = corNamVc[corPexVi],
                 pos = 3)

            if(opl@typeC != "PCA" && length(typeVc) == 1)
                legend("topleft",
                       pch = c(18, 15),
                       legend = c("X vars", "Y vars"))

        } else if(substr(ploC, 1, 7) == "predict") {

            abline(0, 1)

            text(ploMN[, 1:2],
                 cex = parCexN,
                 col = ploColVc,
                 labels = ploNamVc)

            if(!is.null(obsLegVc))
                if(ploC == "predict-train") {
                    .legendF(obsLegVc,
                             ploMN)

                } else
                    .legendF(tesLegVc,
                             ploMN)

        } else if(ploC == "x-loading") {

            points(ploMN,
                   col = "grey",
                   pch = 18)

            points(ploMN[pexVi, ],
                   pch = 18,
                   col = "black")

            ## pexLabVc <- colnames(opl@suppLs[["xModelMN"]])[pexVi]
            pexLabVc <- rownames(opl@loadingMN)[pexVi]
            pexLabVc[duplicated(pexLabVc)] <- ""

            text(ploMN[pexVi, ],
                 cex = parCexN,
                 col = "black",
                 labels = pexLabVc,
                 pos = rep(c(4, 2, 3, 1), each = opl@suppLs[["topLoadI"]]))

        } else if(ploC == "x-score") {

            lines(sqrt(var(ploMN[, 1]) * hotFisN) * cos(radVn),
                  sqrt(var(ploMN[, 2]) * hotFisN) * sin(radVn))
            ## Tenenhaus98, p87

            if(!is.null(obsLegVc))
                .legendF(obsLegVc,
                         ploMN)

            text(ploMN,
                 cex = parCexN,
                 col = ploColVc,
                 labels = obsLabVc)

            pu1N <- par("usr")[1]
            pu2N <- par("usr")[2]

            cexRqcN <- ifelse(layL, 0.7, 1)

            mtext(paste("R2X", round(opl@summaryDF[, "R2X(cum)"], 3), sep = "\n"),
                  at = pu1N * ifelse(layL, 1.35, 1.1),
                  cex = cexRqcN,
                  font = 1,
                  line = 3,
                  side = 1)


            if(parEllipsesL) {
                par(lwd = 2)
                for(colC in unique(ploColVc))
                    .ellipseF(ploMN[ploColVc == colC, , drop = FALSE],
                              colC = colC)
            }

            if(!is.null(opl@suppLs[["yMCN"]])) {

                mtext(paste("R2Y", round(opl@summaryDF[, "R2Y(cum)"], 3), sep = "\n"),
                      at = pu1N * ifelse(layL, 1, 0.8),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("Q2Y", round(opl@summaryDF[, "Q2(cum)"], 3), sep = "\n"),
                      at = pu1N * ifelse(layL, 0.6, 0.4),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("RMSEE", round(opl@summaryDF[, "RMSEE"], 3), sep = "\n"),
                      at =  -pu1N * ifelse(layL, 0.6, 0.4),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("pre", opl@summaryDF[, "pre"], sep = "\n"),
                      at = -pu1N * ifelse(layL, 0.92, 0.7),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                if(opl@summaryDF[, "ort"] > 0)
                    mtext(paste("ort", opl@summaryDF[, "ort"], sep = "\n"),
                          at = -pu1N * ifelse(layL, 1.1, 0.9),
                          cex = cexRqcN,
                          font = 1,
                          line = 3,
                          side = 1)

            }

        } else if(ploC == "xy-score") {

            abline(0, 1)

            if(!is.null(obsLegVc))
                .legendF(obsLegVc,
                         ploMN)

            text(ploMN,
                 cex = parCexN,
                 col = ploColVc,
                 labels = obsLabVc)

        } else if(ploC == "xy-weight") {

            text(ploMN[, 1:2],
                 cex = parCexN,
                 col = ploColVc,
                 labels = c(rownames(opl@weightStarMN), rownames(opl@cMN)))

            if(!layL)
                legend("topleft",
                       col = c("grey", "black"),
                       legend = c("X", "Y"),
                       text.col = c("grey", "black"))

        }

    } ## end of ploPclF()


    if(is.null(tCompMN) && ploC %in% c("correlation",
                                      "outlier",
                                      "x-loading",
                                      "x-score",
                                      "xy-weight"))
        warning("No ", ploC, " plotting", call. = FALSE)

    ## Hotteling's T2 (Tenenhaus98, p86)
    ##----------------------------------

    if(!is.null(tCompMN))
        hotFisN <- (nrow(tCompMN) - 1) * 2 * (nrow(tCompMN)^2 - 1) / (nrow(tCompMN) * nrow(tCompMN) * (nrow(tCompMN) - 2)) * qf(0.95, 2, nrow(tCompMN) - 2)


    radVn <- seq(0, 2 * pi, length.out = 100)


    if(ploC == "outlier") {

        ## Observation diagnostics
        ## see Hubert2005 p66
        ##------------------------

        mahInvCovMN <- solve(cov(tCompMN))

        pcaResMN <- cbind(sdsVn = apply(tCompMN,
                              1,
                              function(x) sqrt(t(as.matrix(x)) %*% mahInvCovMN %*% as.matrix(x))),
                          odsVn = apply(opl@suppLs[["xModelMN"]] - tcrossprod(tCompMN, pCompMN),
                              1,
                              function(x) sqrt(drop(crossprod(x[complete.cases(x)])))))

        pcaResThrVn <- c(sqrt(qchisq(0.975, 2)),
                         (mean(pcaResMN[, 2]^(2/3)) + sd(pcaResMN[, 2]^(2/3)) * qnorm(0.975))^(3/2))

        pcaResExtVi <- union(which(pcaResMN[, 1] > pcaResThrVn[1]),
                             which(pcaResMN[, 2] > pcaResThrVn[2]))

        plot(pcaResMN,
             main = "Observation diagnostics",
             type = "n",
             xlab = "Score distance (SD)",
             xlim = c(0, max(pcaResMN[, 1]) * 1.1),
             xpd = TRUE,
             ylab = "Orthgonal distance (OD)",
             ylim = c(0, max(pcaResMN[, 2]) * 1.1))
        abline(v = pcaResThrVn[1],
               lty = "dashed")
        abline(h = pcaResThrVn[2],
               lty = "dashed")

        if(length(pcaResExtVi)) {

            points(pcaResMN[-pcaResExtVi, , drop = FALSE],
                   col = obsColVc[-pcaResExtVi],
                   pch = 18)
            text(pcaResMN[pcaResExtVi, , drop = FALSE],
                 cex = parCexN,
                 col = obsColVc[pcaResExtVi],
                 labels = obsLabVc[pcaResExtVi])

        } else
            points(pcaResMN,
                   col = obsColVc,
                   pch = 18)

    } ## outlier


    if(ploC == "overview") {

        if(opl@typeC == "PCA") {

            barplot(opl@modelDF[, "R2X"] * 100,
                    main = "Variance explained",
                    names.arg = rownames(opl@modelDF),
                    xlab = "PC",
                    ylab = "% of total variance")


        } else {

            if(opl@summaryDF[, "ort"] == 0) {
                modBarDF <- opl@modelDF
            } else
                modBarDF <- opl@modelDF[!(rownames(opl@modelDF) %in% c("rot", "sum")), ]

            barplot(rbind(modBarDF[, "R2Y(cum)"],
                          modBarDF[, "Q2(cum)"]),
                    beside = TRUE,
                    main = "Model overview",
                    names.arg = rownames(modBarDF),
                    xlab = "")

            axis(2, lwd=2, lwd.ticks=1)

            abline(h = 0.5,
                   col = "gray")

            barplot(rbind(modBarDF[, "R2Y(cum)"],
                          modBarDF[, "Q2(cum)"]),
                    add = TRUE,
                    beside = TRUE,
                    col = c("grey", "black"))

            text(1.5,
                 0,
                 adj = c(0, 0.5),
                 col = "white",
                 srt = 90,
                 labels = " R2Y") ## R2Ycum

            text(2.5,
                 0,
                 adj = c(0, 0.5),
                 col = "white",
                 labels = " Q2Y", ## Q2cum
                 srt = 90)

        }

    } ## overview


    if(ploC == "permutation") {

        perDiaVc <- c("R2Y(cum)", "Q2(cum)")

        plot(c(min(opl@suppLs[["permMN"]][, "sim"]), 1),
             c(min(opl@suppLs[["permMN"]][, c("R2Y(cum)", "Q2(cum)")]), 1),
             main = paste0("pR2Y = ",
                 opl@summaryDF[, "pR2Y"],
                 ", pQ2 = ",
                 opl@summaryDF[, "pQ2"]),
             type = "n",
             xlab = expression(Similarity(bold(y), bold(y[perm]))),
             ylab = "")

        points(opl@suppLs[["permMN"]][, "sim"], opl@suppLs[["permMN"]][, "Q2(cum)"],
               col = "black",
               pch = 18)
        abline(h = opl@suppLs[["permMN"]][1, "Q2(cum)"],
               col = "black")
        points(opl@suppLs[["permMN"]][, "sim"], opl@suppLs[["permMN"]][, "R2Y(cum)"],
               col = "grey",
               pch = 18)
        abline(h = opl@suppLs[["permMN"]][1, "R2Y(cum)"],
               col = "grey")
        .legendF(c("R2Y", "Q2Y"),
                 "bottomright",
                 colVc = c("grey", "black"))

    } ## permutation


    if(ploC == "x-variance") {

        par(las=2)

        boxplot(opl@suppLs[["xSubIncVarMN"]],
                main = "X variances (min, med, max)",
                names=rep("", 3),
                xaxt="n",
                yaxt="n")

        axis(at=axTicks(2),
             side=2)

        mtext(substr(colnames(opl@suppLs[["xSubIncVarMN"]]), 1, 9),
              cex=ifelse(layL, 0.7, 0.8),
              at=1:3,
              line=0.2,
              side=1)

        par(las=0)

    } ## "x-variance"


    if(ploC %in% c("correlation",
                   "predict-test",
                   "predict-train",
                   "x-loading",
                   "x-score",
                   "xy-score",
                   "xy-weight"))
        ploPclF()


} ## .plotF


## Transforms a character or numeric vector into colors
.colorF <- function(namVcn) {

    ## 16 color palette without 'gray'
    palVc <- c("blue", "red", "green3", "cyan", "magenta", "#FF7F00", "#6A3D9A", "#B15928", "aquamarine4", "yellow4", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#FFFF99")

    if(is.null(namVcn) || all(is.na(namVcn))) {

        if(!is.null(namVcn)) {

            dev.new()

            palNamVc <- paste0(1:length(palVc),
                               "_",
                               palVc)

            pie(rep(1, length(palVc)),
                col = palVc,
                labels = palNamVc)

            print(matrix(palNamVc, ncol = 1))

        }

        return(palVc)

    } else {

        if(is.character(namVcn)) {

            namFcn <- factor(namVcn)

            if(length(levels(namFcn)) <= length(palVc)) {
                scaVc <- palVc[1:length(levels(namFcn))]
            } else
                scaVc <- c(palVc,
                           rep("gray",
                               length(levels(namFcn)) - length(palVc)))

            names(scaVc) <- levels(namFcn)

            colVc <- scaVc[unlist(sapply(namVcn,
                                         function(scaleC) {
                                             if(is.na(scaleC))
                                                 return(NA)
                                             else
                                                 which(levels(namFcn) == scaleC)
                                         }))]

        } else if(is.numeric(namVcn)) {

            scaVc <- rev(rainbow(100, end = 4/6))
            if(length(namVcn) > 1) {
                colVc <- scaVc[round((namVcn - min(namVcn, na.rm = TRUE)) / diff(range(namVcn, na.rm = TRUE)) * 99) + 1]
            } else
                colVc <- rep("black", length(namVcn))

        } else
            stop("'namVcn' argument must be a vector of either character or numeric mode", call. = FALSE)

        colVc[is.na(colVc)] <- "grey"
        names(colVc) <- namVcn

    }

    return(list(colVc = colVc,
                scaVc = scaVc))

}  ## end of .colorF()


## Draws Mahalanobis ellipse
.ellipseF <- function(xMN,
                      colC = NULL,
                      sxyMN = NULL) {
    ## Adapted from the 'drawMahal' function of the 'chemometrics' package
    ## by P. Filzmoser and K. Varmuza

    if(ncol(xMN) != 2)
        stop("Matrix must have two columns", call. = FALSE)

    radVn <- seq(0, 2 * pi, length.out = 100)

    csqN <- qchisq(0.95, 2) ## ncol(xMN) == 2

    xMeaVn <- colMeans(xMN)
    ##        t1        t2
    ## 1.1771851 0.5661031

    xCovMN <- sxyMN

    if(is.null(xCovMN))
        xCovMN <- cov(xMN)
    ##            t1         t2
    ## t1  1.8079514 -0.9768156
    ## t2 -0.9768156  1.0201432

    xCovSvdLs <- svd(xCovMN, nv = 0)
    ## $ d: num [1:2] 2.467 0.361
    ## $ u: num [1:2, 1:2] -0.829 0.559 0.559 0.829
    ## $ v: NULL

    if(!is.null(colC)) {

        mahMN <- matrix(1, nrow = length(radVn), ncol = 1) %*% xMeaVn + cbind(cos(radVn), sin(radVn)) %*% diag(sqrt(xCovSvdLs[["d"]] * csqN)) %*% t(xCovSvdLs[["u"]])
        lines(mahMN,
              col = colC)

    } else {

        zerVarVin <- which(xCovSvdLs[["d"]] < .Machine$double.eps)

        if(length(zerVarVin))
            stop("Covariance matrix cannot be inverted because of ", length(zerVarVin), " zero eigen values\n", call. = FALSE)
        else
            sxyInvMN <- xCovSvdLs[["u"]] %*% diag(1 / xCovSvdLs[["d"]]) %*% t(xCovSvdLs[["u"]])

        invisible(sxyInvMN)

    }

} ## end of .ellipseF()


## Plots the figure legend
.legendF <- function(namOrLegVcn,
                     locCMN = "topright",
                     txtCexN = 0.7,
                     colVc = NULL) {
    ## Note:
    ##  locCMN: either a character indicating the corner of the plot where the legend is to be plotted or the numeric matrix of point coordinates for the legLocF function below to find the corner where there is most space

    ## Determining the location (corner) for the legend

    legLocF <- function(thrN = 0.2) {

        lefN <- par("usr")[1] + thrN * diff(par("usr")[1:2])
        rigN <- par("usr")[2] - thrN * diff(par("usr")[1:2])
        topN <- par("usr")[4] - thrN * diff(par("usr")[3:4])
        botN <- par("usr")[3] + thrN * diff(par("usr")[3:4])

        locVl <- c(all(ploMN[, 1] > lefN |
                       ploMN[, 2] < topN),
                   all(ploMN[, 1] < rigN |
                       ploMN[, 2] < topN),
                   all(ploMN[, 1] > lefN |
                       ploMN[, 2] > botN),
                   all(ploMN[, 1] < rigN |
                       ploMN[, 2] > botN))
        names(locVl) <- c("topleft", "topright", "bottomleft", "bottomright")

        return(locVl)

    }

    stopifnot(is.character(locCMN) || is.matrix(locCMN))

    if(is.matrix(locCMN)) {

        ploMN <- locCMN

        thrVn <- seq(0, 0.25, by = 0.05)
        locSumVn <- sapply(thrVn, function(thrN) sum(legLocF(thrN = thrN)))

        if(sum(locSumVn) > 0)
            locC <- names(which(legLocF(thrVn[max(which(locSumVn > 0))]))[1])
        else
            locC <- "topleft"

    } else
        locC <- locCMN


    ## Determining the color scale

    if(!is.null(colVc)) {
        scaVc <- colVc
        names(scaVc) <- namOrLegVcn
    } else
        scaVc <- .colorF(namOrLegVcn)[["scaVc"]]

    legTypC <- ifelse(is.character(namOrLegVcn), "cha", "num")

    ## Plotting the legend

    dpx <- diff(par("usr")[1:2])
    dpy <- diff(par("usr")[3:4])
    pu1 <- par("usr")[1]
    pu2 <- par("usr")[2]
    pu3 <- par("usr")[3]
    pu4 <- par("usr")[4]

    if(locC == "topright") {

        xLefN <- pu2 - 0.05 * dpx
        xRigN <- pu2 - 0.02 * dpx
        yBotN <- pu4 - 0.22 * dpy
        yTopN <- pu4 - 0.02 * dpy

    } else if(locC == "topleft") {

        xLefN <- pu1 + 0.02 * dpx
        xRigN <- pu1 + 0.05 * dpx
        yBotN <- pu4 - 0.22 * dpy
        yTopN <- pu4 - 0.02 * dpy

    } else if(locC == "bottomleft") {

        xLefN <- pu1 + 0.02 * dpx
        xRigN <- pu1 + 0.05 * dpx
        yBotN <- pu3 + 0.02 * dpy
        yTopN <- pu3 + 0.22 * dpy

    } else if(locC == "bottomright") {

        xLefN <- pu2 - 0.05 * dpx
        xRigN <- pu2 - 0.02 * dpx
        yBotN <- pu3 + 0.02 * dpy
        yTopN <- pu3 + 0.22 * dpy

    }

    yVn <- seq(yBotN,
               yTopN,
               length = length(scaVc) + 1)

    yBotVn <- yVn[1:(length(yVn) - 1)]
    yTopVn <- yVn[2:length(yVn)]

    rect(xleft = xLefN,
         ybottom = yBotVn,
         xright = xRigN,
         ytop = yTopVn,
         col = scaVc,
         border = NA)

    xLegN <- ifelse(grepl("left", locC), xRigN, xLefN)
    xAdjN <- ifelse(grepl("left", locC), 0, 1)

    if(legTypC == "cha") {
        text(xLegN,
             seq(yBotN + (yTopN - yBotN) / (2 * length(scaVc)),
                 yTopN - (yTopN - yBotN) / (2 * length(scaVc)),
                 length = length(scaVc)),
             adj = c(xAdjN, 0.5),
             cex = txtCexN,
             col = scaVc,
             labels = names(scaVc))
    } else
        text(xLegN,
             seq(yBotN,
                 yTopN,
                 length = 5),
             adj = c(xAdjN, 0.5),
             cex = txtCexN,
             labels = signif(seq(min(namOrLegVcn, na.rm = TRUE), max(namOrLegVcn, na.rm = TRUE), length = 5), 2))

} ## .legendF


.similarityF <- function(x, y,
                         .char2numF,
                         charL = FALSE) {

    if(charL) {
        return(sum(x == y) / length(x))
    } else
        return(cor(x, y, use = "pairwise.complete.obs"))

} ## .similarityF


setMethod("fitted", "opls",
          function(object, ...) {

              if(!is.null(object@suppLs[["yPreMN"]])) {

                  if(mode(object@suppLs[["yMCN"]]) == "character") {

                      yPredMCN <- object@suppLs[[".char2numF"]](object@suppLs[["yPreMN"]],
                                                                c2nL = FALSE)

                      if(is.vector(object@suppLs[["y"]])) {
                          fit <- c(yPredMCN)
                          names(fit) <- rownames(yPredMCN)
                      } else if(is.factor(object@suppLs[["y"]])) {
                          fit <- c(yPredMCN)
                          names(fit) <- rownames(yPredMCN)
                          fit <- factor(fit, levels = levels(object@suppLs[["y"]]))
                      } else if(is.matrix(object@suppLs[["y"]])) {
                          fit <- yPredMCN
                      } else
                          stop() ## this case should not happen

                  } else {

                      yPredMCN <- object@suppLs[["yPreMN"]]

                      if(is.vector(object@suppLs[["y"]])) {
                          fit <- c(yPredMCN)
                          names(fit) <- rownames(yPredMCN)
                      } else if(is.matrix(object@suppLs[["y"]])) {
                          fit <- yPredMCN
                      } else
                          stop() ## this case should not happen

                  }

                  return(fit)

              } else
                  return(NULL)

          }) ## fitted

setGeneric(name = "tested",
           def = function(object, ...) standardGeneric("tested"))

setMethod("tested", "opls",
          function(object) {

              if(!is.null(object@suppLs[["yTesMN"]])) {

                  if(mode(object@suppLs[["yMCN"]]) == "character") {

                      yTestMCN <- object@suppLs[[".char2numF"]](object@suppLs[["yTesMN"]],
                                                                             c2nL = FALSE)
                      if(is.vector(object@suppLs[["y"]])) {
                          test <- c(yTestMCN)
                          names(test) <- rownames(yTestMCN)
                      } else if(is.factor(object@suppLs[["y"]])) {
                          test <- c(yTestMCN)
                          names(test) <- rownames(yTestMCN)
                          test <- factor(test, levels = levels(object@suppLs[["y"]]))
                      } else if(is.matrix(object@suppLs[["y"]])) {
                          test <- yTestMCN
                      } else
                          stop() ## this case should not happen

                  } else {

                      yTestMCN <- object@suppLs[["yTesMN"]]

                      if(is.vector(object@suppLs[["y"]])) {
                          test <- c(yTestMCN)
                          names(test) <- rownames(yTestMCN)
                      } else if(is.matrix(object@suppLs[["y"]])) {
                          test <- yTestMCN
                      } else
                          stop() ## this case should not happen

                  }

                  return(test)

              } else
                  stop("Test results only available for (O)PLS(-DA) models", call. = FALSE)


          })


setMethod("coef", "opls",
          function(object, ...) {
              return(object@coefficientMN)
          }) ## coef


setMethod("residuals", "opls",
          function(object, ...) {

    if(grepl("PLS", object@typeC)) {

        fit <- fitted(object)

        if(length(object@subsetVi) == 0) {
            trainVi <- 1:length(fit)
        } else
            trainVi <- object@subsetVi

        if(mode(object@suppLs[["yMCN"]]) == "numeric") {
            y <- object@suppLs[["y"]]
            if(is.matrix(y))
                y <- y[trainVi, , drop = FALSE]
            else
                y <- y[trainVi]
            return(y - fit)
        } else
            return(as.numeric(as.character(c(object@suppLs[["y"]])[trainVi]) != as.character(c(fit))))

    } else
        stop("'residuals' defined for (O)PLS(-DA) regression models only", call. = FALSE)

}) ## residuals


setMethod("predict", "opls",
          function(object, newdata, ...) {

              if(object@typeC == "PCA")
                  stop("Predictions currently available for (O)PLS(-DA) models only (not PCA)",
                       call. = FALSE)

              if(missing(newdata)) {

                  return(fitted(object))

              } else {

                  if(is.data.frame(newdata)) {
                      if(!all(sapply(newdata, data.class) == "numeric")) {
                          stop("'newdata' data frame must contain numeric columns only", call. = FALSE)
                      } else
                          newdata <- as.matrix(newdata)
                  } else if(is.matrix(newdata)) {
                      if(mode(newdata) != "numeric")
                          stop("'newdata' matrix must be of 'numeric' mode", call. = FALSE)
                  } else
                      stop("'newdata' must be either a data.frame or a matrix", call. = FALSE)

                  if(ncol(newdata) != as.numeric(object@descriptionMC["X_variables", ])) {
                      if(length(object@xZeroVarVi) == 0) {
                          stop("'newdata' number of variables is ",
                               ncol(newdata),
                               " whereas the number of variables used for model training was ",
                               as.numeric(object@descriptionMC["X_variables", ]),
                               ".",
                               call. = FALSE)
                      } else if(ncol(newdata) - as.numeric(object@descriptionMC["X_variables", ]) ==
                                as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ])) {
                          warning(as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ]),
                                  " near zero variance variables excluded during the model training will be removed from 'newdata'.",
                                  call. = FALSE)
                          newdata <- newdata[, -object@xZeroVarVi, drop = FALSE]
                      } else {
                          stop("'newdata' number of variables (",
                               ncol(newdata),
                               ") does not correspond to the number of initial variables (",
                               as.numeric(object@descriptionMC["X_variables", ]),
                               ") minus the number of near zero variance variables excluded during the training (",
                               as.numeric(object@descriptionMC["near_zero_excluded_X_variables", ]),
                               ").",
                               call. = FALSE)
                      }
                  }

                  xteMN <- scale(newdata, object@xMeanVn, object@xSdVn)

                  if(object@summaryDF[, "ort"] > 0) {

                      for(noN in 1:object@summaryDF[, "ort"]) {
                          if(object@suppLs[["naxL"]]) {
                              xtoMN <- matrix(0, nrow = nrow(xteMN), ncol = 1)
                              for(i in 1:nrow(xtoMN)) {
                                  comVl <- complete.cases(xteMN[i, ])
                                  xtoMN[i, ] <- crossprod(xteMN[i, comVl], object@orthoWeightMN[comVl, noN]) / drop(crossprod(object@orthoWeightMN[comVl, noN]))
                              }
                          } else
                              xtoMN <- xteMN %*% object@orthoWeightMN[, noN]

                          xteMN <- xteMN - tcrossprod(xtoMN, object@orthoLoadingMN[, noN])
                      }

                  }

                  if(object@suppLs[["naxL"]]) {
                      yTesScaMN <- matrix(0, nrow = nrow(xteMN), ncol = ncol(object@coefficientMN),
                                          dimnames = list(rownames(xteMN), colnames(object@coefficientMN)))
                      for(j in 1:ncol(yTesScaMN))
                          for(i in 1:nrow(yTesScaMN)) {
                              comVl <- complete.cases(xteMN[i, ])
                              yTesScaMN[i, j] <- crossprod(xteMN[i, comVl], object@coefficientMN[comVl, j])
                          }
                  } else
                      yTesScaMN <- xteMN %*% object@coefficientMN

                  ## if(object@suppLs[["nayL"]])
                  ##     yTesScaMN <- yTesScaMN[!is.na(yMCN[testVi, ]), , drop = FALSE]

                  yTesMN <- scale(scale(yTesScaMN,
                                        FALSE,
                                        1 / object@ySdVn),
                                  -object@yMeanVn,
                                  FALSE)
                  attr(yTesMN, "scaled:center") <- NULL
                  attr(yTesMN, "scaled:scale") <- NULL

                  if(is.factor(fitted(object))) {

                      yTestMCN <- object@suppLs[[".char2numF"]](yTesMN,
                                                                     c2nL = FALSE)
                      predMCNFcVcn <- as.character(yTestMCN)
                      names(predMCNFcVcn) <- rownames(newdata)
                      predMCNFcVcn <- factor(predMCNFcVcn, levels = levels(object@suppLs[["y"]]))

                  } else if(is.vector(fitted(object))) {

                      if(is.character(fitted(object))) {

                          yTestMCN <- object@suppLs[[".char2numF"]](yTesMN,
                                                                         c2nL = FALSE)
                          predMCNFcVcn <- as.character(yTestMCN)
                          names(predMCNFcVcn) <- rownames(newdata)

                      } else {

                          predMCNFcVcn <- as.numeric(yTesMN)
                          names(predMCNFcVcn) <- rownames(newdata)

                      }

                  } else if(is.matrix(fitted(object))) {

                      if(mode(fitted(object)) == "character") {
                          predMCNFcVcn  <- object@suppLs[[".char2numF"]](yTesMN,
                                                                         c2nL = FALSE)
                      } else
                          predMCNFcVcn <- yTesMN

                      dimnames(predMCNFcVcn) <- list(rownames(newdata), "pred")

                  }

                  return(predMCNFcVcn)

              }

          }) ## predict

setGeneric("getSummaryDF",
           function(object, ...) {standardGeneric("getSummaryDF")})

setMethod("getSummaryDF", "opls",
          function(object) {
                  return(object@summaryDF)
          })

setGeneric("getPcaVarVn",
           function(object, ...) {standardGeneric("getPcaVarVn")})

setMethod("getPcaVarVn", "opls",
          function(object) {
                  return(object@pcaVarVn)
          })

setGeneric("getScoreMN",
           function(object, ...) {standardGeneric("getScoreMN")})

setMethod("getScoreMN", "opls",
          function(object, orthoL = FALSE) {
              if(orthoL)
                  return(object@orthoScoreMN)
              else
                  return(object@scoreMN)
          })

setGeneric("getLoadingMN",
           function(object, ...) {standardGeneric("getLoadingMN")})

setMethod("getLoadingMN", "opls",
          function(object, orthoL = FALSE) {
              if(orthoL)
                  return(object@orthoLoadingMN)
              else
                  return(object@loadingMN)
          })

setGeneric("getWeightMN",
           function(object, ...) {standardGeneric("getWeightMN")})

setMethod("getWeightMN", "opls",
          function(object, orthoL = FALSE) {
              if(orthoL)
                  return(object@orthoWeightMN)
              else
                  return(object@weightMN)
          })

setGeneric("getVipVn",
           function(object, ...) {standardGeneric("getVipVn")})

setMethod("getVipVn", "opls",
          function(object, orthoL = FALSE) {
              if(orthoL)
                  return(object@orthoVipVn)
              else
                  return(object@vipVn)
          })

setGeneric("getSubsetVi",
           function(object, ...) {standardGeneric("getSubsetVi")})

setMethod("getSubsetVi", "opls",
          function(object) {
                  return(object@subsetVi)
          })






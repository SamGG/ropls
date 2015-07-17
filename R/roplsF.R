roplsF <- function(xMN,
                   yMCN = NULL,
                   predI = NA,
                   orthoI = 0,
                   plotVc = c("correlation",
                       "none",
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
                       "xy-weight")[8],

                   algoC = c("default", "nipals", "svd")[1],
                   crossvalI = 7,
                   fileFig.pdfC = NULL,
                   fileInfo.txtC = NULL,
                   log10L = FALSE,
                   permI = 0,
                   scaleC = c("center",
                       "pareto",
                       "standard")[3],
                   testVi = NULL,
                   verboseC = c("all",
                       "none",
                       "overview",
                       "summary")[1],

                   parAsColVcn = NA,
                   parCexN = 1,
                   parCompVi = c(1, 2),
                   parDevNewL = TRUE,
                   parEllipsesL = NA,
                   parLabVc = NA,
                   parTitleL = TRUE,
                   parTopLoadI = 3) {


    ## Option
    ##---------------

    strAsFacL <- options()$stringsAsFactors
    options(stingsAsFactors = FALSE)


    ## Constants
    ##---------------

    epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16

    if(!is.null(yMCN)) {
        c2nVc <- sort(unique(drop(yMCN))) ## NAs are removed
        c2nVn <- as.numeric(factor(c2nVc))
        names(c2nVn) <- c2nVc

        n2cVc <- names(c2nVn)
        names(n2cVc) <- c2nVn

        c2nLs <- list(c2nVn = c2nVn,
                      n2cVc = n2cVc)
    }


    ## Log file (in case of integration into Galaxy)
    ##----------------------------------------------

    if(!is.null(fileInfo.txtC))
        sink(fileInfo.txtC, append = TRUE)


    ## Checking
    ##---------------

    if(!is.matrix(xMN))
        stop("'xMN' must be a matrix", call. = FALSE)

    if(mode(xMN) != 'numeric')
        stop("'xMN' must be of 'numeric' type", call. = FALSE)

    if(any(apply(xMN, 2, function(colVn) all(is.na(colVn)))))
        stop("'xMN' contains columns with 'NA' only", call. = FALSE)

    if(!is.logical(log10L))
        stop("'log10L' must be a logical", call. = FALSE)

    if(permI < 0 || (permI - floor(permI)) > 1e-10)
        stop("'permI' must be an integer", call. = FALSE)

    if(permI > 0 && (is.null(yMCN) || ncol(yMCN) > 1))
        stop("Permutation testing available only when 'yMCN' has one column only", call. = FALSE)

    if(permI > 0 && !is.null(testVi))
        stop("Permutation testing and external validation cannot be performed simultaneously", call. = FALSE)
    ## compatibility between permutation testing and external cross-validation to be verified

    if(!is.null(yMCN) && !is.matrix(yMCN))
        stop("'yMCN' must be a matrix", call. = FALSE)

    if(!is.null(yMCN) && !(mode(yMCN) %in% c("character", "numeric")))
        stop("'yMCN' must be of 'character' or 'numeric' type", call. = FALSE)

    if(!is.null(yMCN) && nrow(xMN) != nrow(yMCN))
        stop("'xMN' and 'yMCN' must have the same number of rows", call. = FALSE)

    if(!(algoC %in% c('default', 'nipals', 'svd')))
        stop("'algoC' must be either 'default', 'nipals', or 'svd'", call. = FALSE)

    if(algoC == "default")
        algoC <- ifelse(is.null(yMCN) && !any(is.na(c(xMN))), "svd", "nipals")

    if(!is.null(yMCN) && algoC != "nipals")
        stop("'nipals' algorithm must be used for (O)PLS(-DA)", call. = FALSE)

    if((is.na(orthoI) || orthoI > 0) && ncol(yMCN) > 1)
        stop("OPLS(-DA) only available for 'yMCN' with single column only", call. = FALSE)

    if(!(length(scaleC) == 1 && scaleC %in% c('none', 'center', 'pareto', 'standard')))
        stop("'scaleC' must be either 'none', 'center', 'pareto', or 'standard'", call. = FALSE)

    if(!is.null(testVi) && (is.null(yMCN) || ncol(yMCN) > 1))
        stop("external validation only available for 'yMCN' with single column only'", call. = FALSE)

    if(!is.null(testVi) &&
       !(mode(testVi) == 'character' && testVi == 'odd') &&
       !all(testVi %in% 1:nrow(xMN)))
        stop("'testVi' must be either set to 'odd' or an integer vector of 'xMN' row numbers", call. = FALSE)

    if(!is.null(yMCN) &&
       (is.na(orthoI) || orthoI > 0) &&
       any(is.na(yMCN)))
        stop("'yMCN' must not contain missing values for OPLS(-DA)", call. = FALSE)

    if(length(verboseC) != 1 || !(verboseC %in% c('all', 'none', 'overview', 'summary')))
        stop("'verboseC' must be either 'all', 'none', 'overview', or 'summary'", call. = FALSE)

    if(crossvalI > nrow(xMN))
        stop("'crossvalI' must be less than the row number of 'xMN'", call. = FALSE)

    if(is.na(orthoI) || orthoI > 0)
        if(is.na(predI) || predI > 1) {
            predI <- 1
            warning("OPLS: number of predictive components ('predI' argument) set to 1", call. = FALSE)
        }

    if(!is.na(predI))
        if(predI > min(nrow(xMN), ncol(xMN))) {
            predI <- min(nrow(xMN), ncol(xMN))
            warning("'predI' set to the minimum of the 'xMN' matrix dimensions: ", predI, call. = FALSE)
        }



    ##------------------------------------
    ##   Computations
    ##------------------------------------

    cvaVc <- paste("cv", 1:crossvalI, sep = "")
    cvaLs <- vector(length = length(cvaVc), mode = "list")
    names(cvaLs) <- cvaVc
    cvaOutLs <- split(1:nrow(xMN), rep(1:crossvalI, length = nrow(xMN)))

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

    if(!is.null(testVi) &&
       mode(testVi) == "character" &&
       testVi == "odd") {

        if(mode(yMCN) == "numeric")
            testVi <- seq(1, nrow(xMN), by = 2)
        else {
            testVi <- integer()
            for(claC in unique(drop(yMCN)))
                testVi <- c(testVi,
                           which(drop(yMCN) == claC)[seq(2, sum(drop(yMCN) == claC), by = 2)])
            testVi <- sort(testVi)
        }
    }


    ## Filtering out zero variance variables
    ##--------------------------------------

    xVarIndLs <- list()
    xVarIndLs[[1]] <- 1:nrow(xMN)

    if(!is.null(testVi)) {
        xVarIndLs[[1]] <- setdiff(1:nrow(xMN), testVi)
    } else if(!is.null(yMCN) && ncol(yMCN) == 1 && nrow(xMN) >= 2 * crossvalI)
        for(cvkN in 1:crossvalI)
            xVarIndLs <- c(xVarIndLs, list(setdiff(1:nrow(xMN), cvaOutLs[[cvkN]])))

    xVarVarLs <- lapply(xVarIndLs,
                        function(xVarVi) {
                            apply(xMN[xVarVi, ],
                                  2,
                                  function(colVn) var(colVn, na.rm = TRUE))
                        })

    xZeroVarVi <- integer()
    for(k in 1:length(xVarVarLs))
        xZeroVarVi <- union(xZeroVarVi, which(xVarVarLs[[k]] < epsN))

    if(length(xZeroVarVi) > 0) {
        names(xZeroVarVi) <- colnames(xMN)[xZeroVarVi]
        xMN <- xMN[, -xZeroVarVi]
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

    ropLs <- .coreF(xMN = xMN,
                    yMCN = yMCN,
                    orthoI = orthoI,
                    predI = predI,
                    scaleC = scaleC,
                    algoC = algoC,
                    crossvalI = crossvalI,
                    testVi = testVi,
                    c2nLs = c2nLs,
                    xZeroVarVi = xZeroVarVi)

    tCompMN <- NULL
    pCompMN <- NULL
    topLoadMN <- NULL

    permMN <- NULL



    ## RMSER

    ## if(!is.null(yMCN) && ncol(yMCN) == 1 && !is.null(ropLs[["testVi"]])) {

    ##     refVi <- setdiff(1:nrow(xMN), ropLs[["testVi"]])
    ##     ynrMCN <- yMCN

    ##     ynrMCN[refVi, ] <- sample(ynrMCN[refVi, ])

    ##     rndLs <- .coreF(xMN = xMN,
    ##                     yMCN = ynrMCN,
    ##                     orthoI = ropLs[["orthoI"]],
    ##                     predI = ropLs[["predI"]],
    ##                     scaleC = scaleC,
    ##                     algoC = algoC,
    ##                     crossvalI = crossvalI,
    ##                     testVi = ropLs[["testVi"]],
    ##                     c2nLs = c2nLs,
    ##                     xZeroVarVi = ropLs[["xZeroVarVi"]])

    ##     ropLs[["summaryDF"]][, "RMSER"] <- rndLs[["summaryDF"]][, "RMSEP"]

    ## }

    ## Permutation testing (Szymanska et al, 2012)

    if(permI > 0) {

        modSumVc <- colnames(ropLs[["summaryDF"]])

        permMN <- matrix(0,
                        nrow = 1 + permI,
                        ncol = length(modSumVc),
                        dimnames = list(NULL, modSumVc))

        perSimVn <- numeric(1 + permI)
        perSimVn[1] <- 1


        permMN[1, ] <- as.matrix(ropLs[["summaryDF"]])

        if(is.null(fileInfo.txtC))
            bar <- txtProgressBar(1, permI, style = 3)


        for(k in 1:permI) {

            if(is.null(fileInfo.txtC))
                setTxtProgressBar(bar, k)

            yVcn <- drop(yMCN)
            if(is.null(ropLs[["testVi"]])) {
                yPerVcn <- sample(yVcn)
            } else {
                yPerVcn <- numeric(nrow(xMN))
                refVi <- setdiff(1:nrow(xMN), ropLs[["testVi"]])
                yPerVcn[refVi] <- sample(yVcn[refVi])
                yPerVcn[ropLs[["testVi"]]] <- yVcn[ropLs[["testVi"]]]
            }
            yPerMCN <- matrix(yPerVcn, ncol = 1)

            perLs <- .coreF(xMN = xMN,
                            yMCN = yPerMCN,
                            orthoI = ropLs[["orthoI"]],
                            predI = ropLs[["predI"]],
                            scaleC = scaleC,
                            algoC = algoC,
                            crossvalI = crossvalI,
                            testVi = ropLs[["testVi"]],
                            c2nLs = c2nLs,
                            xZeroVarVi = ropLs[["xZeroVarVi"]])

            permMN[1 + k, ] <- as.matrix(perLs[["summaryDF"]])

            perSimVn[1 + k] <- .similarityF(yMCN, yPerMCN, c2nLs = c2nLs, charL = mode(yMCN) == "character")

        }

        if(is.null(fileInfo.txtC))
            close(bar)

        permMN <- cbind(permMN, sim = perSimVn)

        perPvaVn <- c(pR2Y = (1 + length(which(permMN[-1, "R2Y(cum)"] >= permMN[1, "R2Y(cum)"]))) / (nrow(permMN) - 1),
                      pQ2 = (1 + length(which(permMN[-1, "Q2(cum)"] >= permMN[1, "Q2(cum)"]))) / (nrow(permMN) - 1))
        ropLs[["summaryDF"]][, "pR2Y"] <- perPvaVn["pR2Y"]
        ropLs[["summaryDF"]][, "pQ2"] <- perPvaVn["pQ2"]

    }

    if(is.null(yMCN)) {
        ropTypC <- "PCA"
    } else {
        if(ncol(yMCN) > 1 || mode(yMCN) == "numeric")
            ropTypC <- "PLS"
        else
            ropTypC <- "PLS-DA"
    }
    if(ropLs[["orthoI"]] > 0)
        ropTypC <- paste("O", ropTypC, sep = "")


    ropLs <- c(list(typC = ropTypC),
               ropLs,
               list(permMN = permMN,
                    cvaLs = cvaLs))


    ##------------------------------------
    ##   Numerical results
    ##------------------------------------

    if(verboseC == "all") {

        ## General information

        message("'roplsF' call: ", as.character(Sys.time()))

        message("Number of observations: ", nrow(xMN))
        message("Number of X variables: ", ncol(xMN))
        message("Number of 'near zero' excluded X variables: ", length(ropLs[["xZeroVarVi"]]))

    }

    totN <- length(c(xMN))
    nasN <- sum(is.na(c(xMN)))

    if(!is.null(yMCN)) {
        if(verboseC == "all")
            message("Number of Y variables: ", ncol(yMCN))
        totN <- totN + length(c(yMCN))
        nasN <- nasN + sum(is.na(c(yMCN)))
    }

    if(verboseC == "all") {

        message("Number of missing values: ", nasN, " (", round(nasN / totN * 100), "%)")

        ## Raw summary
        ##------------

        message("Summary of the ", parTopLoadI, " increasing variance spaced raw variables:")

    }

    if(ncol(xMN) > parTopLoadI) {
        xVarVn <- apply(xMN, 2, var)
        names(xVarVn) <- 1:length(xVarVn)
        xVarVn <- sort(xVarVn)
        xVarSorVin <- as.numeric(names(xVarVn[seq(1, length(xVarVn), length = parTopLoadI)]))
        xRawMN <- xMN[, xVarSorVin]
    } else
        xRawMN <- xMN

    if(verboseC == "all")
        print(summary(xRawMN))


    if(ncol(xMN) < 100) {

        if(verboseC == "all")
            message("Correlations between the X-variables:")

        xCorMN <- cor(xMN, use = "pairwise.complete.obs")
        xCorMN[lower.tri(xCorMN, diag = TRUE)] <- 0

        if(ncol(xMN) > parTopLoadI) {

            xCorNexDF <- which(abs(xCorMN) > sort(abs(xCorMN), decreasing = TRUE)[parTopLoadI + 1],
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

        rm(xCorMN)

        if(verboseC == "all")
            print(signif(xCorDisMN, 2))

    }

    if(verboseC == "all") {

        if(scaleC != "none")
            message("X: mean-centering",
                    switch(scaleC,
                           pareto = " and pareto scaling",
                           standard = " and unit-variance scaling"))

        if(is.null(yMCN)) {

            message(ropLs[["typC"]], " ('", algoC, "' algorithm)")
            message("Number of components: ", ropLs[["predI"]])

        } else {

            message("Number of Y variables: ", ncol(ropLs[["yPreMN"]]))
            message("Y: mean-centering",
                    switch(scaleC,
                           pareto = " and pareto scaling",
                           standard = " and unit-variance scaling"))

            if(ropLs[["orthoI"]] > 0) {
                message(ropLs[["typC"]], " ('nipals' algorithm)")
                message("Number of orthogonal components: ", ropLs[["orthoI"]])
            } else
                message("PLS ('nipals' algorithm)")

            message("Number of predictive components: ", ropLs[["predI"]])

        }


        message("Number of reference observations: ",
                nrow(ropLs[["xModelMN"]]),
                " (",
                round(nrow(ropLs[["xModelMN"]]) / nrow(xMN) * 100),
                "%)")

        message("Correlations between variables and components:")

    }

    if(ropLs[["predI"]] + ropLs[["orthoI"]] < 2) {

        if(length(plotVc) > 1 || plotVc != "none") {
            warning("A single component model has been selected by cross-validation: Only the 'overview' plot is available", call. = FALSE)
            plotVc <- "overview"
        }

        ropLs[["tCompMN"]] <- ropLs[["tMN"]]
        ropLs[["pCompMN"]] <- ropLs[["pMN"]]

    } else {

        if(ropLs[["orthoI"]] > 0) {
            if(parCompVi[2] > ropLs[["orthoI"]] + 1)
                stop("Selected orthogonal component for plotting (ordinate) exceeds the total number of orthogonal components of the model", call. = FALSE)
            ropLs[["tCompMN"]] <- cbind(ropLs[["tMN"]][, 1], ropLs[["tOrthoMN"]][, parCompVi[2] - 1])
            ropLs[["pCompMN"]] <- cbind(ropLs[["pMN"]][, 1], ropLs[["pOrthoMN"]][, parCompVi[2] - 1])
            colnames(ropLs[["pCompMN"]]) <- colnames(ropLs[["tCompMN"]]) <- c("h1", paste("o", parCompVi[2] - 1, sep = ""))
        } else {
            if(max(parCompVi) > ropLs[["predI"]])
                stop("Selected component for plotting as ordinate exceeds the total number of predictive components of the model", call. = FALSE)
            ropLs[["tCompMN"]] <- ropLs[["tMN"]][, parCompVi, drop = FALSE]
            ropLs[["pCompMN"]] <- ropLs[["pMN"]][, parCompVi, drop = FALSE]
        }

    }

    cxtCompMN <- cor(ropLs[["xModelMN"]], ropLs[["tCompMN"]],
                     use = "pairwise.complete.obs")

    if(!is.null(ropLs[["yModelMN"]]))
        cytCompMN <- cor(ropLs[["yModelMN"]], ropLs[["tCompMN"]], use = "pairwise.complete.obs")


    if(parTopLoadI * 4 < ncol(ropLs[["xModelMN"]])) {

        pexVin <- integer(parTopLoadI * ncol(ropLs[["pCompMN"]]) * 2) ## 'ex'treme values

        for(k in 1:ncol(ropLs[["pCompMN"]])) {

            pkVn <-  ropLs[["pCompMN"]][, k]

            pexVin[1:(2 * parTopLoadI) + 2 * parTopLoadI * (k - 1)] <- c(order(pkVn)[1:parTopLoadI],
                                                                         rev(order(pkVn, decreasing = TRUE)[1:parTopLoadI]))

        }

    } else
        pexVin <- 1:ncol(ropLs[["xModelMN"]])


    pxtCompMN <- cbind(ropLs[["pCompMN"]],
                       cxtCompMN)

    if(ncol(ropLs[["pCompMN"]]) == 1) {
       colnames(pxtCompMN)[2] <- paste0("cor_", colnames(pxtCompMN)[2])
    } else
        colnames(pxtCompMN)[3:4] <- paste0("cor_", colnames(pxtCompMN)[3:4])

    ropLs[["topLoadMN"]] <- pxtCompMN

    ropLs[["topLoadMN"]] <- ropLs[["topLoadMN"]][pexVin, , drop = FALSE]

    if(parTopLoadI * 4 < ncol(ropLs[["xModelMN"]]) &&
       ncol(ropLs[["pCompMN"]]) > 1) {

        ropLs[["topLoadMN"]][(2 * parTopLoadI + 1):(4 * parTopLoadI), c(1, 3)] <- NA
        ropLs[["topLoadMN"]][1:(2 * parTopLoadI), c(2, 4)] <- NA

    }

    if(verboseC == "all")
        print(signif(ropLs[["topLoadMN"]], 2))

    if(verboseC %in% c("all", "overview")) {

        message("Model overview:")

        optDigN <- options()[["digits"]]
        options(digits = 3)
        print(ropLs[["modelDF"]])
        options(digits = optDigN)

    }

    if(verboseC %in% c("all", "summary")) {

        message("Model summary:")

        optDigN <- options()[["digits"]]
        options(digits = 3)
        print(ropLs[["summaryDF"]])
        options(digits = optDigN)

    }

    if(length(plotVc) > 1 || plotVc != "none") {


        ##------------------------------------
        ##   Graphics
        ##------------------------------------


        if("summary" %in% plotVc)
            plotVc <- c(ifelse(permI > 0, "permutation", "overview"),
                       "outlier",
                       "x-score",
                       "x-loading")


        ## Checking arguments
        ##-------------------

        if(!all(plotVc %in% c('correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight')))
            stop("'plotVc' elements must be either 'correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight'", call. = FALSE)

        if('predict-test' %in% plotVc && is.null(ropLs[["testVi"]]))
            stop("For the 'predict-test' graphic to be generated, 'testVi' must not be kept to NULL", call. = FALSE)

        if(!any(is.na(parLabVc))) {
            if(length(parLabVc) != nrow(xMN))
                stop("'parLabVc' vector length must be equal to the number of 'xMN' rows")
            if(mode(parLabVc) != "character")
                stop("'parLabVc' must be of 'character' type")
        }

        if(!any(is.na(parAsColVcn))) {
            if(length(parAsColVcn) != nrow(xMN))
                stop("'parAsColVcn' vector length must be equal to the number of 'xMN' rows")
            if(!(mode(parAsColVcn) %in% c("character", "numeric")))
                stop("'parAsColVcn' must be of 'character' or 'numeric' type")
        }

        if(permI == 0 && 'permutation' %in% plotVc)
            stop("'permI' must be > 0 for 'permutation' graphic to be plotted", call. = FALSE)

        if(ropLs[["orthoI"]] > 0)
            if(parCompVi[1] != 1) {
                parCompVi[1] <- 1
                warning("OPLS: first component to display ('parCompVi' first value) set to 1", call. = FALSE)
            }

        if("xy-weight" %in% plotVc &&
           (is.null(yMCN) || is.na(ropLs[["orthoI"]]) || ropLs[["orthoI"]] > 0))
            stop("'xy-weight graphic can be displayed only for PLS(-DA) models", call. = FALSE)

        if(any(grepl('predict', plotVc)) && (is.null(yMCN) || ncol(yMCN) != 1))
            stop("'yMCN' must have a single column for 'predict' graphics", call. = FALSE)

        if(is.na(parEllipsesL)) {
            if((is.null(yMCN) && !any(is.na(parAsColVcn))) || ## PCA case
               (!is.null(yMCN) && ncol(yMCN) == 1 && mode(yMCN) == "character")) { ## (O)PLS-DA cases
                parEllipsesL <- TRUE
            } else
                parEllipsesL <- FALSE
        } else
            if(parEllipsesL && !is.null(yMCN) && mode(yMCN) != "character")
                stop("'yMCN' must have a single column of character type for mahalanobis ellipses to be plotted", call. = FALSE)


        ## Observation and variable names and colors
        ##------------------------------------------

        ## obsLabVc

        if(!any(is.na(parLabVc))) {
            obsLabVc <- parLabVc
        } else if(!is.null(rownames(xMN))) {
            obsLabVc <- rownames(xMN)
        } else
            obsLabVc <- as.character(1:nrow(xMN))

        if(!is.null(ropLs[["testVi"]])) {
            tesLabVc <- obsLabVc[ropLs[["testVi"]]]
            obsLabVc <- obsLabVc[-ropLs[["testVi"]]]
        } else
            tesLabVc <- ""

        ## obsColVc

        if(!any(is.na(parAsColVcn))) {
            obsColVc <- .colorF(parAsColVcn)[["colVc"]]
            obsLegVc <- parAsColVcn
        } else if(!is.null(yMCN) && ncol(yMCN) == 1) {
            obsColVc <- .colorF(c(yMCN))[["colVc"]]
            obsLegVc <- c(yMCN)
        } else {
            obsColVc <- rep("black", nrow(xMN))
            obsLegVc <- NULL
        }

        if(!is.null(ropLs[["testVi"]])) {
            tesColVc <- obsColVc[ropLs[["testVi"]]]
            obsColVc <- obsColVc[-ropLs[["testVi"]]]
            if(!is.null(obsLegVc)) {
                tesLegVc <- obsLegVc[ropLs[["testVi"]]]
                obsLegVc <- obsLegVc[-ropLs[["testVi"]]]
            }
        }


        ## Layout
        ##-------

        if(!parDevNewL && length(plotVc) != 1)
            stop("'plotVc' must be of length 1 when 'parDevNewL' is set to FALSE", call. = FALSE)

        if(parDevNewL) {
            layRowN <- ceiling(sqrt(length(plotVc)))
            if(is.null(fileFig.pdfC))
                dev.new()
            else
                pdf(fileFig.pdfC)
            layout(matrix(1:layRowN^2, byrow = TRUE, nrow = layRowN))
        }

        layL <- !parDevNewL || length(plotVc) > 1


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

        for(ploC in plotVc)
            .plotF(ploC,
                   ropLs = ropLs,
                   c2nLs = c2nLs,
                   cxtCompMN = cxtCompMN,
                   cytCompMN = cytCompMN,
                   obsColVc = obsColVc,
                   obsLabVc = obsLabVc,
                   obsLegVc = obsLegVc,
                   layL = layL,
                   parCexN = parCexN,
                   parEllipsesL = parEllipsesL,
                   parTitleL = parTitleL,
                   parTopLoadI = parTopLoadI,
                   parCompVi = parCompVi,
                   pexVin = pexVin,
                   plotVc = plotVc,
                   tesColVc = tesColVc,
                   tesLabVc = tesLabVc,
                   tesLegVc = tesLegVc,
                   xMN = xMN,
                   xRawMN = xRawMN,
                   yMCN = yMCN)

        if(layL)
            par(font=1, font.axis=1, font.lab=1, lwd=1,
                mar=c(5.1, 4.1, 4.1, 2.1),
                pch=1)

        if(!is.null(fileFig.pdfC))
            dev.off()

    } ## if(length(plotVc) > 1 || plotVc != "none")


    ropLs[["xModelMN"]] <- NULL
    ropLs[["yModelMN"]] <- NULL

    if("yPreMN" %in% names(ropLs)) {

        if(mode(yMCN) == "character") {
            ropLs[["yPredMCN"]] <- .char2numF(ropLs[["yPreMN"]], c2nLs = c2nLs,
                                   c2nL = FALSE)
        } else
            ropLs[["yPredMCN"]] <- ropLs[["yPreMN"]]

        ropLs[["yPreMN"]] <- NULL

    }

    if("yTesMN" %in% names(ropLs)) {

        if(mode(yMCN) == "character") {
            ropLs[["yTestMCN"]] <- .char2numF(ropLs[["yTesMN"]], c2nLs = c2nLs,
                                              c2nL = FALSE)
        } else
            ropLs[["yTestMCN"]] <- ropLs[["yTesMN"]]

        ropLs[["yTesMN"]] <- NULL

    }

    ## Warning messages
    ##-----------------

    if(verboseC == "all")
        warnings()

    if(!is.null(fileInfo.txtC))
        sink(NULL)


    ## Returning
    ##----------

    options(stringsAsFactors = strAsFacL)

    return(invisible(ropLs))


} ## end of roplsF


.char2numF <- function(inpMCN,
                       c2nLs,
                       c2nL = TRUE) {

    if(c2nL)
        outMCN <- matrix(as.vector(c2nLs[["c2nVn"]][drop(inpMCN)]), ncol = 1)
    else
        outMCN <- matrix(as.vector(c2nLs[["n2cVc"]][as.character(round(drop(.DQ2F(inpMCN, c2nLs = c2nLs))))]), ncol = 1)

    dimnames(outMCN) <- dimnames(inpMCN)
    return(outMCN)

}


## Transforms a character or numeric vector into colors
.colorF <- function(namVcn) {

    ## 16 color palette without 'gray'
    palVc <- c("black", "red", "green3", "blue", "cyan", "magenta", "#FF7F00", "#6A3D9A", "#B15928", "aquamarine4", "yellow4", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#FFFF99")

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


## Core algorithms for PCA, PLS(-DA), and OPLS(-DA)
.coreF <- function(xMN,
                   yMCN,
                   orthoI,
                   predI,
                   scaleC,
                   algoC,
                   crossvalI,
                   testVi,
                   c2nLs,
                   xZeroVarVi) {

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
    yPreMN <- NULL   ## (O)PLS only
    yTesMN <- NULL   ## (O)PLS only
    tOrthoMN <- NULL   ## OPLS only
    pOrthoMN <- NULL   ## OPLS only
    wOrthoMN <- NULL   ## OPLS only


    ## Missing values
    ##---------------


    naxVi <- which(is.na(c(xMN)))
    naxL <- length(naxVi) > 0

    if(algoC == "svd" && length(which(is.na(c(xMN)))) > 0) {
        minN <- min(c(xMN[!is.na(xMN)])) / 2
        xMN[is.na(xMN)] <- minN
        warning("Missing values set to ", round(minN, 1), " (half minimum value) for 'svd' algorithm to be used", call. = FALSE)
    }

    if(!is.null(yMCN)) {
        nayVi <- which(is.na(c(yMCN)))
        nayL <- length(nayVi) > 0
    }

    ## yMCN 'character' to 'numeric' conversion + .errorF function

    yMN <- yMCN

    if(!is.null(yMCN)) {

        if(mode(yMCN) == "character")
            yMN <- .char2numF(yMCN, c2nLs = c2nLs)

        ## training and a test partition

        if(!is.null(testVi)) {

            xTesMN <- xMN[testVi, , drop = FALSE]
            xMN <- xMN[-testVi, , drop = FALSE]
            yMN <- yMN[-testVi, , drop = FALSE]

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

    if(is.na(orthoI)) {
        orthoI <- 14
        predI <- 1
        autNcoL <- TRUE
    } else if(is.na(predI)) {
        if(orthoI > 0)
            predI <- 1
        else {
            if(is.null(yMCN))
                predI <- min(nrow(xMN), ncol(xMN))
            else
                predI <- 15
        }
        autNcpL <- TRUE
    }


    ##------------------------------------
    ##   Preprocessing
    ##------------------------------------


    ## X variable variances
    ##---------------------

    xVarVn <- apply(xMN, 2, function(colVn) var(colVn, na.rm = TRUE))


    ## X-Scaling
    ##---------------

    if(scaleC != "none") {

        xMeanVn <- apply(xMN, 2, function(colVn) mean(colVn, na.rm = TRUE))

        if(scaleC == "center")
            xSdVn <- rep(1, times = ncol(xMN))
        else if(scaleC == "pareto")
            xSdVn <- apply(xMN, 2, function(colVn) sqrt(sd(colVn, na.rm = TRUE)))
        else if(scaleC == "standard")
            xSdVn <- apply(xMN, 2, function(colVn) sd(colVn, na.rm = TRUE))

        xMN <- scale(xMN, center = xMeanVn, scale = xSdVn)


    } else {

        xMeanVn <- rep(0, times = ncol(xMN))
        xSdVn <- rep(1, times = ncol(xMN))

    }

    if(!is.null(colnames(xMN))) {
        xvaNamVc <- colnames(xMN)
    } else
        xvaNamVc <- paste("x", 1:ncol(xMN), sep = "")

    predIamVc <- paste("h", 1:predI, sep = "")

    pMN <- matrix(0,
                  nrow = ncol(xMN),
                  ncol = predI,
                  dimnames = list(xvaNamVc, predIamVc))

    tMN <- uMN <- matrix(0,
                         nrow = nrow(xMN),
                         ncol = predI,
                         dimnames = list(obsNamVc, predIamVc))

    ssxTotN <- sum(xMN^2, na.rm = TRUE)

    if(is.null(yMCN)) {


        ##------------------------------------
        ##   PCA
        ##------------------------------------


        varVn <- numeric(predI)
        names(varVn) <- predIamVc
        vSumVn <- sum(apply(xMN, 2, function(y) var(y, na.rm = TRUE))) ## xMN is centered

        modelDF <- as.data.frame(matrix(0,
                                      nrow = predI,
                                      ncol = 3,
                                      dimnames = list(predIamVc, c("R2X", "R2X(cum)", "Iter."))))

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

                   tMN <- tMN[, 1:predI]
                   pMN <- pMN[, 1:predI]
                   varVn <- varVn[1:predI]
                   rownames(tMN) <- obsNamVc
                   rownames(pMN) <- xvaNamVc
                   names(varVn) <- colnames(pMN) <- colnames(tMN) <- predIamVc

                   modelDF[, "R2X"] <- round(varVn / vSumVn, 3)

               }) ## svd

        modelDF[, "R2X(cum)"] <- cumsum(modelDF[, "R2X"])

        if(autNcpL) {
            vSelVl <- cumsum(varVn) / vSumVn > 0.5
            if(sum(vSelVl) == 0)
                predI <- max(length(varVn), 2)
            else
                predI <- max(which(vSelVl)[1], 2)
            tMN <- tMN[, 1:predI]
            pMN <- pMN[, 1:predI]
            varVn <- varVn[1:predI]
            modelDF <- modelDF[1:predI, ]
        }

        summaryDF <- modelDF[predI, c("R2X(cum)"), drop = FALSE]


    } else { ## if(is.null(yMCN))


        ## Y-Scaling
        ##---------------

        if(scaleC != "none") {

            yMeanVn <- apply(yMN, 2, function(colVn) mean(colVn, na.rm = TRUE))

            if(scaleC == "center")
                ySdVn <- rep(1, times = ncol(yMN))
            else if(scaleC == "pareto")
                ySdVn <- apply(yMN, 2, function(colVn) sqrt(sd(colVn, na.rm = TRUE)))
            else if(scaleC == "standard")
                ySdVn <- apply(yMN, 2, function(colVn) sd(colVn, na.rm = TRUE))

            yMN <- scale(yMN, center = yMeanVn, scale = ySdVn)

        } else {

            yMeanVn <- rep(0, times = ncol(yMN))
            ySdVn <- rep(1, times = ncol(yMN))

        }

        if(!is.null(colnames(yMN))) {
            yvaNamVc <- colnames(yMN)
        } else
            yvaNamVc <- paste("y", 1:ncol(yMN), sep = "")


        wMN <- pMN
        uMN <- tMN

        cMN <- matrix(0,
                      nrow = ncol(yMN),
                      ncol = predI,
                      dimnames = list(yvaNamVc, predIamVc))




        ## Cross-validation variables

        cvfNamVc <- paste("cv", 1:crossvalI, sep = "")
        cvfOutLs <- split(1:nrow(xMN), rep(1:crossvalI, length = nrow(xMN)))
        prkVn <- numeric(crossvalI)

        ## rules

        ru1ThrN <- ifelse(nrow(xMN) > 100, yes = 0, no = 0.05)

        ## ssyTotN <- rs0N <- sum((yMN - mean(yMN, na.rm = TRUE))^2, na.rm = TRUE)
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
                                          dimnames = list(predIamVc, c("R2X", "R2X(cum)", "R2Y", "R2Y(cum)", "Q2", "Q2(cum)", "Signif.", "Iter."))))
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

                ## cross-validation

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
                } else {
                    modelDF[hN, "Signif."] <- "R1"
                }

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
                    stop("No model was built because the first predictive component was already not significant;\nSelect a number of predictive components of 2 if you want the algorithm to compute a model despite this.", call. = FALSE)

                wMN <- wMN[, 1:hN, drop = FALSE]
                tMN <- tMN[, 1:hN, drop = FALSE]
                pMN <- pMN[, 1:hN, drop = FALSE]
                cMN <- cMN[, 1:hN, drop = FALSE]
                uMN <- uMN[, 1:hN, drop = FALSE]

                predIamVc <- predIamVc[1:hN]

                predI <- hN

                modelDF <- modelDF[1:hN, ]

            }

            summaryDF <- modelDF[predI, c("R2X(cum)", "R2Y(cum)", "Q2(cum)")]


            ## Rotation matrix (W*)

            if(predI == 1) {

                rMN <- wMN

            } else {

                pwMN <- crossprod(pMN, wMN)
                rMN <- wMN %*% solve(pwMN)
                colnames(rMN) <- predIamVc

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


            if(!is.null(testVi)) {
                yActMCN <- yMCN[-testVi, , drop = FALSE]
            } else
                yActMCN <- yMCN

            if(mode(yMCN) == "character") {
                yActMN <- .char2numF(yActMCN, c2nLs = c2nLs)
            } else
                yActMN <- yActMCN


            summaryDF[, "RMSEE"] <- sqrt(.errorF(yActMN, yPreMN)^2 * nrow(yActMN) / (nrow(yActMN) - (1 + predI))) ## for SIMCA compatibility


            if(!is.null(testVi)) { ## tRaining/tEst partition

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
                    yTesScaMN <- yTesScaMN[!is.na(yMCN[testVi, ]), , drop = FALSE]

                yTesMN <- scale(scale(yTesScaMN,
                                       FALSE,
                                       1 / ySdVn),
                                 -yMeanVn,
                                 FALSE) ## predicted values
                attr(yTesMN, "scaled:center") <- NULL
                attr(yTesMN, "scaled:scale") <- NULL

                if(mode(yMCN) == "character") {
                    yTestMCN <- .char2numF(yTesMN, c2nLs = c2nLs, c2nL = FALSE)
                } else
                    yTestMCN <- yTesMN

                yTesActMCN <- yMCN[testVi, , drop = FALSE] ## actual values
                if(mode(yMCN) == "character") {
                    yTesActMN <- .char2numF(yTesActMCN, c2nLs = c2nLs)
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

            orthoIamVc <- paste("o", 1:orthoI, sep = "")

            tOrthoMN <- matrix(0,
                            nrow = nrow(xMN),
                            ncol = orthoI,
                            dimnames = list(obsNamVc, orthoIamVc))
            wOrthoMN <- pOrthoMN <- matrix(0,
                                     nrow = ncol(xMN),
                                     ncol = orthoI,
                                     dimnames = list(xvaNamVc, orthoIamVc))
            cOrthoMN <- matrix(0,
                            nrow = ncol(yMN),
                            ncol = orthoI,
                            dimnames = list(yvaNamVc, orthoIamVc))

            modelDF <- as.data.frame(matrix(NA,
                                          nrow = 3 + orthoI,
                                          ncol = 7,
                                          dimnames = list(c("h1", "rot", orthoIamVc, "sum"), c("R2X", "R2X(cum)", "R2Y", "R2Y(cum)", "Q2", "Q2(cum)", "Signif."))))
            for(j in 1:ncol(modelDF))
                mode(modelDF[, j]) <- ifelse(colnames(modelDF)[j] == "Signif.", "character", "numeric")

            xcvLs <- c(lapply(cvfOutLs, function(Vi) xMN[-Vi, , drop = FALSE]),
                       list(xMN))
            xcvTesLs <- c(lapply(cvfOutLs, function(Vi) xMN[Vi, , drop = FALSE]))
            ycvLs <- c(lapply(cvfOutLs, function(Vi) yMN[-Vi, , drop = FALSE]),
                       list(yMN))

            breL <- FALSE

            for(noN in 1:(orthoI + 1)) {

                if(breL)
                    break

                for(cvN in 1:length(xcvLs)) {

                    if(ncol(ycvLs[[cvN]]) > 1) {

                        ## step -|1 [case vector y (p121) | matrix Y (p127)]

                        if(naxL || nayL)
                            wwMN <- apply(ycvLs[[cvN]],
                                          2,
                                          function(colVn) {
                                              wwjVn <- numeric(ncol(xcvLs[[cvN]]))
                                              for(j in 1:ncol(xcvLs[[cvN]])) {
                                                  comVl <- complete.cases(xcvLs[[cvN]][, j]) & complete.cases(colVn)
                                                  wwjVn[j] <- crossprod(xcvLs[[cvN]][comVl,j], colVn[comVl]) / drop(crossprod(colVn[comVl]))
                                              }
                                              wwjVn
                                          })
                        else
                            wwMN <- apply(ycvLs[[cvN]],
                                          2,
                                          function(colVn)
                                          crossprod(xcvLs[[cvN]], colVn) / drop(crossprod(colVn)))

                        ## step -|2

                        wwSvdLs <- svd(wwMN)
                        wwNcpVin <- which(wwSvdLs[["d"]]^2 > epsN * sum(wwSvdLs[["d"]]^2))

                        twMN <- wwSvdLs[["u"]][, wwNcpVin, drop = FALSE] %*% diag(wwSvdLs[["d"]][wwNcpVin], nrow = length(wwNcpVin))

                    }


                    ## step -|4

                    uOldVn <- ycvLs[[cvN]][, 1, drop = FALSE]

                    repeat {

                        ## step 1|5

                        if(naxL || nayL) {
                            wVn <- numeric(ncol(xcvLs[[cvN]]))
                            for(j in 1:ncol(xcvLs[[cvN]])) {
                                comVl <- complete.cases(xcvLs[[cvN]][, j]) &
                                    complete.cases(uOldVn)
                                wVn[j] <- crossprod(xcvLs[[cvN]][comVl, j], uOldVn[comVl]) / drop(crossprod(uOldVn[comVl]))
                            }
                        } else
                            wVn <- crossprod(xcvLs[[cvN]], uOldVn) / drop(crossprod(uOldVn))

                        ## step 2|6

                        wVn <- wVn / sqrt(drop(crossprod(wVn)))

                        ## step 3|7

                        if(naxL) {
                            tVn <- numeric(nrow(xcvLs[[cvN]]))
                            for(i in 1:nrow(xcvLs[[cvN]])) {
                                comVl <- complete.cases(xcvLs[[cvN]][i, ])
                                tVn[i] <- crossprod(xcvLs[[cvN]][i, comVl], wVn[comVl]) / drop(crossprod(wVn[comVl]))
                            }
                        } else
                            tVn <- xcvLs[[cvN]] %*% wVn

                        ## step 4|8

                        if(nayL) {
                            cVn <- numeric(ncol(ycvLs[[cvN]]))
                            for(j in 1:ncol(ycvLs[[cvN]])) {
                                comVl <- complete.cases(ycvLs[[cvN]][, j])
                                cVn[j] <- crossprod(ycvLs[[cvN]][comVl, j], tVn[comVl]) / drop(crossprod(tVn[comVl]))
                            }
                        } else
                            cVn <- crossprod(ycvLs[[cvN]], tVn) / drop(crossprod(tVn))

                        ## step 5|9

                        if(nayL) {
                            uVn <- numeric(nrow(xcvLs[[cvN]]))
                            for(i in 1:nrow(xcvLs[[cvN]])) {
                                comVl <- complete.cases(ycvLs[[cvN]][i, ])
                                uVn[i] <- crossprod(ycvLs[[cvN]][i, comVl], cVn[comVl]) / drop(crossprod(cVn[comVl]))
                            }
                        } else
                            uVn <- ycvLs[[cvN]] %*% cVn / drop(crossprod(cVn))

                        if(nayL) {
                            comVl <- complete.cases(uOldVn)
                            dscN <- drop(sqrt(crossprod((uVn[comVl] - uOldVn[comVl] / uVn[comVl]))))
                        } else
                            dscN <- drop(sqrt(crossprod((uVn - uOldVn) / uVn)))

                        if(ncol(ycvLs[[cvN]]) == 1 ||
                           dscN < 1e-10) {

                            break

                        } else {

                            uOldVn <- uVn

                        }

                    } ## end of repeat

                    ## step 6|

                    if(naxL) {
                        pVn <- numeric(ncol(xcvLs[[cvN]]))
                        for(j in 1:ncol(xcvLs[[cvN]])) {
                            comVl <- complete.cases(xcvLs[[cvN]][, j])
                            pVn[j] <- crossprod(xcvLs[[cvN]][comVl, j], tVn[comVl]) / drop(crossprod(tVn[comVl]))
                        }
                    } else
                        pVn <- crossprod(xcvLs[[cvN]], tVn) / drop(crossprod(tVn))

                    ## step 7|

                    if(ncol(ycvLs[[cvN]]) > 1)
                        for(j in 1:ncol(twMN))
                            wOrthoVn <- pVn - drop(crossprod(twMN[, j, drop = FALSE], pVn)) / drop(crossprod(twMN[, j, drop = FALSE])) * twMN[, j, drop = FALSE]
                    else
                        wOrthoVn <- pVn - drop(crossprod(wVn, pVn)) / drop(crossprod(wVn)) * wVn

                    ## step 8|

                    wOrthoVn <- wOrthoVn / sqrt(drop(crossprod(wOrthoVn)))

                    ## step 9|

                    if(naxL) {
                        tOrthoVn <- numeric(nrow(xcvLs[[cvN]]))
                        for(i in 1:nrow(xcvLs[[cvN]])) {
                            comVl <- complete.cases(xcvLs[[cvN]][i, ])
                            tOrthoVn[i] <- crossprod(xcvLs[[cvN]][i, comVl], wOrthoVn[comVl]) / drop(crossprod(wOrthoVn[comVl]))
                        }
                    } else
                        tOrthoVn <- xcvLs[[cvN]] %*% wOrthoVn

                    ## if(nayL) {
                    ##     cOrthoVn <- numeric(ncol(ynMN))
                    ##     for(j in 1:ncol(ynMN)) {
                    ##         comVl <- complete.cases(ynMN[, j])
                    ##         cOrthoVn[j] <- crossprod(ynMN[comVl, j], tOrthoVn[comVl]) / drop(crossprod(tOrthoVn[comVl]))
                    ##     }
                    ## } else
                    ##     cOrthoVn <- crossprod(ynMN, tOrthoVn) / drop(crossprod(tOrthoVn))

                    ## step 10|

                    if(naxL) {
                        pOrthoVn <- numeric(ncol(xcvLs[[cvN]]))
                        for(j in 1:ncol(xcvLs[[cvN]])) {
                            comVl <- complete.cases(xcvLs[[cvN]][, j])
                            pOrthoVn[j] <- crossprod(xcvLs[[cvN]][comVl, j], tOrthoVn[comVl]) / drop(crossprod(tOrthoVn[comVl]))
                        }
                    } else
                        pOrthoVn <- crossprod(xcvLs[[cvN]], tOrthoVn) / drop(crossprod(tOrthoVn))

                    ## step 12|

                    if(noN == 1)
                        rowI <- 1
                    else
                        rowI <- 1 + noN

                    if(cvN <= crossvalI) { ## cross-validation

                        if(any(is.na(xcvTesLs[[cvN]]))) {
                            prxVn <- numeric(nrow(xcvTesLs[[cvN]]))
                            for(r in 1:length(prxVn)) {
                                comVl <- complete.cases(xcvTesLs[[cvN]][r, ])
                                prxVn[r] <- crossprod(xcvTesLs[[cvN]][r, comVl], wVn[comVl]) / drop(crossprod(wVn[comVl]))
                            }
                            prkVn[cvN] <- sum((yMN[cvfOutLs[[cvN]], , drop = FALSE] - prxVn %*% t(cVn))^2, na.rm = TRUE)
                        } else
                            prkVn[cvN] <- sum((yMN[cvfOutLs[[cvN]], , drop = FALSE] - xcvTesLs[[cvN]] %*% wVn %*% t(cVn))^2, na.rm = TRUE)

                        if(naxL) {
                            tOrthoTesVn <- numeric(nrow(xcvTesLs[[cvN]]))
                            for(i in 1:nrow(xcvTesLs[[cvN]])) {
                                comVl <- complete.cases(xcvTesLs[[cvN]][i, ])
                                tOrthoTesVn[i] <- crossprod(xcvTesLs[[cvN]][i, comVl], wOrthoVn[comVl]) / drop(crossprod(wOrthoVn[comVl]))
                            }
                        } else
                            tOrthoTesVn <- xcvTesLs[[cvN]] %*% wOrthoVn

                        xcvTesLs[[cvN]] <- xcvTesLs[[cvN]] - tcrossprod(tOrthoTesVn, pOrthoVn)

                        if(cvN == crossvalI) {
                            q2N <- 1 - sum(prkVn) / rs0N
                            if(noN == 1)
                                modelDF["h1", "Q2(cum)"] <- modelDF["h1", "Q2"] <- q2N
                            else {
                                modelDF[rowI, "Q2(cum)"] <- q2N - modelDF["h1", "Q2"]
                                modelDF[rowI, "Q2"] <- q2N - sum(modelDF[1:(rowI - 1), "Q2"], na.rm = TRUE)
                            }
                        }

                    } else { ## cvN == crossvalI + 1 (full matrix)

                        if(noN == 1) {
                            if(naxL)
                                modelDF[rowI, "R2X(cum)"] <- modelDF[rowI, "R2X"] <- sum((tcrossprod(tVn, pVn)[!is.na(xMN)])^2) / ssxTotN
                            else
                                modelDF[rowI, "R2X(cum)"] <- modelDF[rowI, "R2X"] <- sum(tcrossprod(tVn, pVn)^2) / ssxTotN
                            if(nayL)
                                modelDF[rowI, "R2Y(cum)"] <- modelDF[rowI, "R2Y"] <- sum((tcrossprod(tVn, cVn)[!is.na(yMN)])^2) / ssyTotN
                            else
                                modelDF[rowI, "R2Y(cum)"] <- modelDF[rowI, "R2Y"] <- sum(tcrossprod(tVn, cVn)^2) / ssyTotN
                        } else {
                            if(nayL)
                                modelDF[rowI, "R2Y(cum)"] <- modelDF[rowI, "R2Y"] <- sum((tcrossprod(tVn, cVn)[!is.na(yMN)])^2) / ssyTotN
                            else
                                modelDF[rowI, "R2Y(cum)"] <- sum(tcrossprod(tVn, cVn)^2) / ssyTotN - modelDF["h1", "R2Y(cum)"]
                            if(noN == 2)
                                modelDF[rowI, "R2Y"] <- modelDF[rowI, "R2Y(cum)"]
                            else
                                modelDF[rowI, "R2Y"] <- modelDF[rowI, "R2Y(cum)"] - modelDF[3, "R2Y(cum)"]
                        }

                        ## cOrthoMN[, noN] <- cOrthoVn
                        if(noN <= orthoI) {
                            pOrthoMN[, noN] <- pOrthoVn
                            tOrthoMN[, noN] <- tOrthoVn
                            wOrthoMN[, noN] <- wOrthoVn
                            if(naxL)
                                modelDF[2 + noN, "R2X"] <- sum((tcrossprod(tOrthoVn, pOrthoVn)[!is.na(xMN)])^2) / ssxTotN
                            else
                                modelDF[2 + noN, "R2X"] <- sum(tcrossprod(tOrthoVn, pOrthoVn)^2) / ssxTotN
                        }

                        if(modelDF[rowI, "R2Y"] < 0.01)
                            modelDF[rowI, "Signif."] <- "N4"
                        else if(modelDF[rowI, "Q2"] < ru1ThrN)
                            modelDF[rowI, "Signif."] <- "NS"
                        else
                            modelDF[rowI, "Signif."] <- "R1"


                        if(autNcoL && modelDF[rowI, "Signif."] != "R1" && rowI > 2) {
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
                        xcvLs[[cvN]] <- xcvLs[[cvN]] - tcrossprod(tOrthoVn, pOrthoVn)

                } ## for(cvN in 1:length(xcvLs)) {

            } ## for(noN in 1:(orthoI + 1)) {

            rm(xcvLs)
            rm(xcvTesLs)
            rm(ycvLs)

            if(naxL) {
                modelDF["rot", "R2X(cum)"] <- sum((tcrossprod(tMN, pMN)[!is.na(xMN)])^2) / ssxTotN
            } else
                modelDF["rot", "R2X(cum)"] <- sum(tcrossprod(tMN, pMN)^2) / ssxTotN
            if(nayL) {
                modelDF["sum", "R2Y(cum)"] <- modelDF["rot", "R2Y(cum)"] <- sum((tcrossprod(tMN, cMN)[!is.na(yMN)])^2) / ssyTotN
            } else
                modelDF["sum", "R2Y(cum)"] <- modelDF["rot", "R2Y(cum)"] <- sum(tcrossprod(tMN, cMN)^2) / ssyTotN

            modelDF["rot", "R2X"] <- modelDF["rot", "R2X(cum)"] - modelDF["h1", "R2X(cum)"]
            modelDF[2 + 1:orthoI, "R2X(cum)"] <- cumsum(modelDF[2 + 1:orthoI, "R2X"])
            modelDF["sum", "R2X(cum)"] <- modelDF["rot", "R2X(cum)"] + modelDF[2 + orthoI, "R2X(cum)"]
            modelDF["sum", "Q2(cum)"] <- modelDF["rot", "Q2(cum)"] <- sum(modelDF[, "Q2"], na.rm = TRUE)

            if(autNcoL) {

                orthoI <- noN - 3

                if(orthoI == 0)
                    stop("No model was built because the first orthogonal component was already not significant;\nSelect a number of orthogonal components of 1 if you want the algorithm to compute a model despite this.", call. = FALSE)

                pOrthoMN <- pOrthoMN[, 1:orthoI, drop = FALSE]
                tOrthoMN <- tOrthoMN[, 1:orthoI, drop = FALSE]
                wOrthoMN <- wOrthoMN[, 1:orthoI, drop = FALSE]

                orthoIamVc <- orthoIamVc[1:orthoI]
                modelDF <- modelDF[c(1:(orthoI + 2), nrow(modelDF)), ]

                if(naxL)
                    modelDF["rot", "R2X(cum)"] <- sum((tcrossprod(tMN, pMN)[!is.na(xMN)])^2) / ssxTotN
                else
                    modelDF["rot", "R2X(cum)"] <- sum(tcrossprod(tMN, pMN)^2) / ssxTotN
                if(nayL)
                    modelDF["sum", "R2Y(cum)"] <- modelDF["rot", "R2Y(cum)"] <- sum((tcrossprod(tMN, cMN)[!is.na(yMN)])^2) / ssyTotN
                else
                    modelDF["sum", "R2Y(cum)"] <- modelDF["rot", "R2Y(cum)"] <- sum(tcrossprod(tMN, cMN)^2) / ssyTotN
                modelDF["rot", "R2X"] <- modelDF["rot", "R2X(cum)"] - modelDF["h1", "R2X(cum)"]
                modelDF["sum", "R2X(cum)"] <- modelDF["rot", "R2X(cum)"] + modelDF[2 + orthoI, "R2X(cum)"]
                modelDF["sum", "Q2(cum)"] <- modelDF["rot", "Q2(cum)"] <- sum(modelDF[, "Q2"], na.rm = TRUE)

            }

            summaryDF <- modelDF["sum", c("R2X(cum)", "R2Y(cum)", "Q2(cum)")]

            ## Rotation matrix (W*)

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

            if(!is.null(testVi)) {
                yActMCN <- yMCN[-testVi, , drop = FALSE]
            } else
                yActMCN <- yMCN

            if(mode(yMCN) == "character") {
                yActMN <- .char2numF(yActMCN, c2nLs = c2nLs)
            } else
                yActMN <- yActMCN

            summaryDF[, "RMSEE"] <- sqrt(.errorF(yActMN, yPreMN)^2 * nrow(yActMN) / (nrow(yActMN) - (1 + predI + orthoI)))


            if(!is.null(testVi)) { ## tRaining/tEst partition

                xteMN <- scale(xTesMN, xMeanVn, xSdVn)

                for(noN in 1:orthoI) {
                    if(naxL) {
                        xtoMN <- matrix(0, nrow = nrow(xteMN), ncol = 1)
                        for(i in 1:nrow(xtoMN)) {
                            comVl <- complete.cases(xteMN[i, ])
                            xtoMN[i, ] <- crossprod(xteMN[i, comVl], wOrthoMN[comVl, noN]) / drop(crossprod(wOrthoMN[comVl, noN]))
                        }
                    } else
                        xtoMN <- xteMN %*% wOrthoMN[, noN]

                    xteMN <- xteMN - tcrossprod(xtoMN, pOrthoMN[, noN])
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
                    yTesScaMN <- yTesScaMN[!is.na(yMCN[testVi, ]), , drop = FALSE]

                yTesMN <- scale(scale(yTesScaMN,
                                      FALSE,
                                      1 / ySdVn),
                                -yMeanVn,
                                FALSE)
                attr(yTesMN, "scaled:center") <- NULL
                attr(yTesMN, "scaled:scale") <- NULL

                if(mode(yMCN) == "character") {
                    yTestMCN <- .char2numF(yTesMN,
                                           c2nLs = c2nLs, c2nL = FALSE)
                } else
                    yTestMCN <- yTesMN

                yTesActMCN <- yMCN[testVi, , drop = FALSE] ## actual values
                if(mode(yMCN) == "character") {
                    yTesActMN <- .char2numF(yTesActMCN,
                                            c2nLs = c2nLs)
                } else
                    yTesActMN <- yTesActMCN

                summaryDF[, "RMSEP"] <- .errorF(c(yTesMN), c(yTesActMN))

            } else
                yTestMCN <- NULL


        } ## end of OPLS


        ## VIP (identical for both PLS and O-PLS)
        ## from the VIP.R script by Wehrens and Mevik
        ## http://mevik.net/work/software/VIP.R
        ## Chong I.-G. and Jun C.-H. (2005)
        ## Wold S., Sjostrom M. and Eriksson L. (2001)

        if(orthoI == 0 && (nrow(cMN) == 1)) { ## PLS with single quantitative yMCN
            ## VIP computation only implemented for single-response models

            ssqVn <- c(cMN)^2 * colSums(tMN^2)
            wn2Vn <- colSums(wMN^2)
            vipVn <- sqrt(rowSums(sweep(wMN^2,
                                        2,
                                        ssqVn / wn2Vn,
                                        "*")) * nrow(wMN) / sum(ssqVn))

        } else if(orthoI > 0) ## OPLS and OPLS-DA (single quant. or qual. yMCN)
            vipVn <- sqrt(ncol(xMN)) * abs(wMN[, 1])

    }

    summaryDF[, "pre"] <- predI
    summaryDF[, "ort"] <- orthoI


    ##------------------------------------
    ##   Returning
    ##------------------------------------


    retLs <- list(predI = predI,
                  xMeanVn = xMeanVn,
                  xSdVn = xSdVn,
                  xZeroVarVi = xZeroVarVi,
                  testVi = testVi,
                  tMN = tMN,
                  pMN = pMN,
                  varVn = varVn,
                  modelDF = modelDF,
                  summaryDF = summaryDF,
                  yMeanVn = yMeanVn,
                  ySdVn = ySdVn,
                  wMN = wMN,
                  cMN = cMN,
                  uMN = uMN,
                  rMN = rMN,
                  bMN = bMN,
                  vipVn = vipVn,
                  yPreMN = yPreMN,
                  yTesMN = yTesMN,
                  orthoI = orthoI,
                  tOrthoMN = tOrthoMN,
                  pOrthoMN = pOrthoMN,
                  wOrthoMN = wOrthoMN,
                  xModelMN = xMN,
                  yModelMN = yMN)


} ## .coreF


## Discriminant Q2 (Westerhuis et al, 2008)
.DQ2F <- function(inpMN,
                  c2nLs) {

    epsN <- .Machine[["double.eps"]] ## [1] 2.22e-16

    outMN <- matrix(as.numeric(sapply(drop(inpMN),
                                      function(inpN)
                                      ifelse(is.na(inpN),
                                             return(NA),
                                             ifelse(inpN < head(c2nLs[["c2nVn"]], 1) - 0.5,
                                                    return(head(c2nLs[["c2nVn"]], 1) - 0.5 + epsN),
                                                    ifelse(inpN > tail(c2nLs[["c2nVn"]], 1) + 0.5,
                                                           return(tail(c2nLs[["c2nVn"]], 1) + 0.5 - epsN),
                                                           inpN))))),
                    ncol = 1)
    dimnames(outMN) <- dimnames(inpMN)
    return(outMN)

}


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


.errorF <- function(x, y)
    sqrt(mean(drop((x - y)^2), na.rm = TRUE))


## Plots the figure legend
.legendF <- function(namOrLegVcn,
                     locCMN = "topright",
                     txtCexN = 0.7) {
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

}


.log10F <- function(inpMN) {

    if(length(which(inpMN < 0)) > 0)
        stop("Negative values in the table to be log10 transformed", call. = FALSE)

    zerMN <- inpMN == 0

    inpMN[zerMN] <- 1

    return(log10(inpMN))

} ## .log10F

.plotF <- function(ploC,
                   ropLs,
                   c2nLs,
                   cxtCompMN,
                   cytCompMN,
                   obsColVc,
                   obsLabVc,
                   obsLegVc,
                   layL,
                   parCexN,
                   parEllipsesL,
                   parTitleL,
                   parTopLoadI,
                   parCompVi,
                   pexVin,
                   plotVc,
                   tesColVc,
                   tesLabVc,
                   tesLegVc,
                   xMN,
                   xRawMN,
                   yMCN) {

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

            if(ropLs[["orthoI"]] > 0)
                yLabC <- paste("with tOrtho",
                               parCompVi[2] - 1,
                               sep = "")

            yLimVn <- xLimVn <- c(-1, 1)

            ploMN <- cxtCompMN

            if(!is.null(yMCN))
                ploMN <- rbind(ploMN,
                               cytCompMN)

        } else if(substr(ploC, 1, 7) == "predict") {

            maiC <- paste("Predicted vs Actual",
                          paste(" (",
                                unlist(strsplit(ploC, "-"))[2],
                                ")",
                                sep = ""),
                          sep = "")
            xLabC <- "predicted"
            yLabC <- "actual"

            if(grepl("train", ploC)) {
                ploNamVc <- obsLabVc
                ploColVc <- obsColVc
            } else {
                ploNamVc <- tesLabVc
                ploColVc <- tesColVc
            }

            ypMN <- eval(parse(text = paste("ropLs[['y", switch(unlist(strsplit(ploC, "-"))[2], train = "Pre", test = "Tes"), "MN']]", sep = ""))) ## predicted

            if(is.null(ropLs[["testVi"]]))
                yaMCN <- yMCN ## actual
            else {
                if(grepl("train", ploC))
                    yaMCN <- yMCN[-ropLs[["testVi"]], , drop = FALSE]
                else
                    yaMCN <- yMCN[ropLs[["testVi"]], , drop = FALSE]
            }

            if(mode(yMCN) == "character") {
                yaMN <- .char2numF(yaMCN, c2nLs = c2nLs)
            } else
                yaMN <- yaMCN

            ploMN <- cbind(ypMN,
                           yaMN) ## to be modified (when ncol(yPreMCN) > 1)

        } else if(ploC == "x-loading") {

            maiC <- "Loadings"

            xLabC <- paste("p",
                           parCompVi[1],
                           " (",
                           round(ropLs[["modelDF"]][parCompVi[1], "R2X"] * 100),
                           "%)",
                           sep = "")

            yLabC <- paste("p",
                           parCompVi[2],
                           " (",
                           round(ropLs[["modelDF"]][parCompVi[2], "R2X"] * 100),
                           "%)",
                           sep = "")

            ploMN <- ropLs[["pCompMN"]]

            if(!is.null(yMCN) && ropLs[["orthoI"]] > 0) {
                yLabC <- paste("pOrtho",
                               parCompVi[2] - 1,
                               " (",
                               round(ropLs[["modelDF"]][parCompVi[2] - 1, "R2X"] * 100),
                               "%)",
                               sep = "")
            }

        } else if(ploC == "x-score") {

            maiC <- paste0("Scores (", ropLs[["typC"]], ")")

            xLabC <- paste("t",
                           parCompVi[1],
                           " (",
                           round(ropLs[["modelDF"]][parCompVi[1], "R2X"] * 100),
                           "%)",
                           sep = "")

            yLabC <- paste("t",
                           parCompVi[2],
                           " (",
                           round(ropLs[["modelDF"]][parCompVi[2], "R2X"] * 100),
                           "%)",
                           sep = "")

            ploMN <- ropLs[["tCompMN"]]

            if(!is.null(yMCN) && ropLs[["orthoI"]] > 0) {
                yLabC <- paste("to", parCompVi[2] - 1, sep = "")
            }

            xLimVn <- c(-1, 1) * max(sqrt(var(ploMN[, 1]) * hotFisN), max(abs(ploMN[, 1])))
            yLimVn <- c(-1, 1) *max(sqrt(var(ploMN[, 2]) * hotFisN), max(abs(ploMN[, 2])))

            ploColVc <- obsColVc

        } else if(ploC == "xy-score") {

            maiC <- "XY-Scores"
            xLabC <- paste("t", parCompVi[1], sep = "")
            yLabC <- paste("u/c", parCompVi[1], sep = "")

            ploMN <- cbind(ropLs[["tMN"]][, parCompVi[1]], ropLs[["uMN"]][, parCompVi[1]] / ropLs[["cMN"]][parCompVi[1]])

            ploColVc <- obsColVc

        } else if(ploC == "xy-weight") {

            maiC <- "Weights"
            xLabC <- paste("w*c", parCompVi[1], sep = "")
            yLabC <- paste("w*c", parCompVi[2], sep = "")

            ploMN <- rbind(ropLs[["rMN"]][, parCompVi],
                           ropLs[["cMN"]][, parCompVi])

            pchVn <- rep(17, times = nrow(ploMN))
            ploColVc <- rep("black", times = nrow(ploMN))

            pchVn[(nrow(ropLs[["rMN"]]) + 1):nrow(ploMN)] <- 15
            ploColVc[(nrow(ropLs[["rMN"]]) + 1):nrow(ploMN)] <- "red"

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

            corPexVin <- pexVin
            corPchVn <- rep(18, ncol(ropLs[["xModelMN"]]))
            corNamVc <- colnames(ropLs[["xModelMN"]])


            if(!is.null(yMCN)) {
                corPexVin <- c(corPexVin, (ncol(ropLs[["xModelMN"]]) + 1):nrow(ploMN))
                corPchVn <- c(corPchVn, rep(15, ncol(ropLs[["yModelMN"]])))
                corNamVc <- c(corNamVc, colnames(ropLs[["yModelMN"]]))
            }

            points(ploMN,
                   pch = corPchVn)

            points(ploMN[corPexVin, ],
                   pch = corPchVn[corPexVin],
                   col = "red")

            text(ploMN[corPexVin, ],
                 cex = parCexN,
                 col = "red",
                 labels = corNamVc[corPexVin],
                 pos = 3)

            if(!is.null(yMCN) && length(plotVc) == 1)
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
                   pch = 18)

            points(ploMN[pexVin, ],
                   pch = 18,
                   col = "red")

            pexLabVc <- colnames(ropLs[["xModelMN"]])[pexVin]
            pexLabVc[duplicated(pexLabVc)] <- ""

            text(ploMN[pexVin, ],
                 cex = parCexN,
                 col = "red",
                 labels = pexLabVc,
                 pos = rep(c(4, 2, 3, 1), each = parTopLoadI))

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

            mtext(paste("R2X", round(ropLs[["summaryDF"]][, "R2X(cum)"], 3), sep = "\n"),
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

            if(!is.null(yMCN)) {

                mtext(paste("R2Y", round(ropLs[["summaryDF"]][, "R2Y(cum)"], 3), sep = "\n"),
                      at = pu1N * ifelse(layL, 1, 0.8),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("Q2Y", round(ropLs[["summaryDF"]][, "Q2(cum)"], 3), sep = "\n"),
                      at = pu1N * ifelse(layL, 0.6, 0.4),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("RMSEE", round(ropLs[["summaryDF"]][, "RMSEE"], 3), sep = "\n"),
                      at =  -pu1N * ifelse(layL, 0.6, 0.4),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("pre", ropLs[["predI"]], sep = "\n"),
                      at = -pu1N * ifelse(layL, 0.92, 0.7),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                if(ropLs[["orthoI"]] > 0)
                    mtext(paste("ort", ropLs[["orthoI"]], sep = "\n"),
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
                 labels = c(rownames(ropLs[["rMN"]]), rownames(ropLs[["cMN"]])))

            if(!layL)
                legend("topleft",
                       col = c("black", "red"),
                       legend = c("X", "Y"),
                       text.col = c("black", "red"))

        }

    } ## end of ploPclF()

    if(is.null(ropLs[["tCompMN"]]) && ploC %in% c("correlation",
                                      "outlier",
                                      "x-loading",
                                      "x-score",
                                      "xy-weight"))
        warning("No ", ploC, " plotting", call. = FALSE)

    ## Hotteling's T2 (Tenenhaus98, p86)
    ##----------------------------------

    if(!is.null(ropLs[["tCompMN"]]))
        hotFisN <- (nrow(ropLs[["tCompMN"]]) - 1) * 2 * (nrow(ropLs[["tCompMN"]])^2 - 1) / (nrow(ropLs[["tCompMN"]]) * nrow(ropLs[["tCompMN"]]) * (nrow(ropLs[["tCompMN"]]) - 2)) * qf(0.95, 2, nrow(ropLs[["tCompMN"]]) - 2)


    radVn <- seq(0, 2 * pi, length.out = 100)


    if(ploC == "outlier") {

        ## Observation diagnostics
        ## see Hubert2005 p66
        ##------------------------

        mahInvCovMN <- solve(cov(ropLs[["tCompMN"]]))

        pcaResMN <- cbind(sdsVn = apply(ropLs[["tCompMN"]],
                              1,
                              function(x) sqrt(t(as.matrix(x)) %*% mahInvCovMN %*% as.matrix(x))),
                          odsVn = apply(ropLs[["xModelMN"]] - tcrossprod(ropLs[["tCompMN"]], ropLs[["pCompMN"]]),
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

        if(is.null(yMCN)) {

            barplot(ropLs[["modelDF"]][, "R2X"] * 100,
                    main = "Variance explained",
                    names.arg = rownames(ropLs[["modelDF"]]),
                    xlab = "PC",
                    ylab = "% of total variance")


        } else {

            if(ropLs[["orthoI"]] == 0)
                modBarDF <- ropLs[["modelDF"]]
            else
                modBarDF <- ropLs[["modelDF"]][!(rownames(ropLs[["modelDF"]]) %in% c("rot", "sum")), ]

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
                    col = .colorF(c("R2Y(cum)", "Q2(cum)"))[["colVc"]])

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

        plot(c(min(ropLs[["permMN"]][, "sim"]), 1),
             c(min(ropLs[["permMN"]][, c("R2Y(cum)", "Q2(cum)")]), 1),
             main = paste0("pR2Y = ",
                 ropLs[["summaryDF"]][, "pR2Y"],
                 ", pQ2 = ",
                 ropLs[["summaryDF"]][, "pQ2"]),
             type = "n",
             xlab=expression(Similarity(bold(y), bold(y[perm]))),
             ylab = "")

        points(ropLs[["permMN"]][, "sim"], ropLs[["permMN"]][, "Q2(cum)"],
               col = "black",
               pch = 18)
        abline(h = ropLs[["permMN"]][1, "Q2(cum)"],
               col = "black")
        points(ropLs[["permMN"]][, "sim"], ropLs[["permMN"]][, "R2Y(cum)"],
               col = "red",
               pch = 18)
        abline(h = ropLs[["permMN"]][1, "R2Y(cum)"],
               col = "red")
        .legendF(c("R2Y", "Q2Y"),
                 "bottomright")

    } ## permutation


    if(ploC == "x-variance") {

        par(las=2)

        boxplot(xRawMN,
                main = "X variances (min, med, max)",
                names=rep("", 3),
                xaxt="n",
                yaxt="n")

        axis(at=axTicks(2),
             side=2)

        mtext(substr(colnames(xRawMN), 1, 9),
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


.similarityF <- function(x, y, c2nLs, charL = FALSE) {

    if(charL) {

        x <- .char2numF(x, c2nLs = c2nLs)
        y <- .char2numF(y, c2nLs = c2nLs)

    }

    return(cor(x, y, use = "pairwise.complete.obs"))

} ## .similarityF

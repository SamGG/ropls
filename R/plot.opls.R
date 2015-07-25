plot.opls <- function(x,
                      y,
                      plotVc = c("correlation",
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
                      parCexN = 1,
                      parCompVi = c(1, 2),
                      parDevNewL = TRUE,
                      parEllipsesL = NA,
                      parLabVc = NA,
                      parTitleL = TRUE,
                      file.pdfC = NULL,
                      ...) {

    if("summary" %in% plotVc)
        plotVc <- c(ifelse(!is.null(x[["suppLs"]][["permMN"]]), "permutation", "overview"),
                    "outlier",
                    "x-score",
                    "x-loading")


    ## Checking arguments
    ##-------------------

    if(!all(plotVc %in% c('correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight')))
        stop("'plotVc' elements must be either 'correlation', 'outlier', 'overview', 'permutation', 'predict-train', 'predict-test', 'x-loading', 'x-score', 'x-variance', 'xy-score', 'xy-weight'", call. = FALSE)

    if('predict-test' %in% plotVc && is.null(x[["subset"]]))
        stop("For the 'predict-test' graphic to be generated, 'subset' must not be kept to NULL", call. = FALSE)

    if(!any(is.na(parLabVc))) {
        if(length(parLabVc) != nrow(x[["scoreMN"]]))
            stop("'parLabVc' vector length must be equal to the number of 'x' rows")
        if(mode(parLabVc) != "character")
            stop("'parLabVc' must be of 'character' type")
    }

    if(!any(is.na(parAsColFcVn))) {
        if(length(parAsColFcVn) != nrow(x[["scoreMN"]]))
            stop("'parAsColFcVn' vector length must be equal to the number of 'x' rows")
        if(!(mode(parAsColFcVn) %in% c("character", "numeric")))
            stop("'parAsColFcVn' must be of 'character' or 'numeric' type")
    }

    if(is.null(x[["suppLs"]][["permMN"]]) && 'permutation' %in% plotVc)
        stop("'permI' must be > 0 for 'permutation' graphic to be plotted", call. = FALSE)

    if(x[["summaryDF"]][, "ort"] > 0)
        if(parCompVi[1] != 1) {
            parCompVi[1] <- 1
            warning("OPLS: first component to display ('parCompVi' first value) set to 1", call. = FALSE)
        }

    if("xy-weight" %in% plotVc &&
       substr(x[["typeC"]], 1, 3) != "PLS")
       ## (is.null(yMCN) || is.na(x[["summaryDF"]][, "ort"]) || x[["summaryDF"]][, "ort"] > 0))
        stop("'xy-weight graphic can be displayed only for PLS(-DA) models", call. = FALSE)

    if(any(grepl('predict', plotVc)) && is.matrix(x[["fitted"]]) && ncol(x[["fitted"]]) > 1)
        ## if(any(grepl('predict', plotVc)) && (is.null(yMCN) || ncol(yMCN) != 1))
        stop("'predict' graphics available for single response models only", call. = FALSE)

    if(is.na(parEllipsesL)) {
        if((x[["typeC"]] == "PCA" && !all(is.na(parAsColFcVn))) || ## PCA case
           grepl("-DA$", x[["typeC"]])) { ## (O)PLS-DA cases
            parEllipsesL <- TRUE
        } else
            parEllipsesL <- FALSE
    } else if(parEllipsesL && !grepl("-DA$", x[["typeC"]]) && all(is.na(parAsColFcVn)))
        stop("Ellipses can be plotted for PCA (or PLS regression) only if the 'parAsColFcVn' is not 'NA'",
             call. = FALSE)

    if(x[["summaryDF"]][, "pre"] + x[["summaryDF"]][, "ort"] < 2) {

        if(length(plotVc) > 1 || plotVc != "overview")
            stop("Only the 'overview' plot is available for single component models", call. = FALSE)

        tCompMN <- x[["scoreMN"]]
        pCompMN <- x[["loadingMN"]]

    } else {

        if(x[["summaryDF"]][, "ort"] > 0) {
            if(parCompVi[2] > x[["summaryDF"]][, "ort"] + 1)
                stop("Selected orthogonal component for plotting (ordinate) exceeds the total number of orthogonal components of the model", call. = FALSE)
            tCompMN <- cbind(x[["scoreMN"]][, 1], x[["orthoScoreMN"]][, parCompVi[2] - 1])
            pCompMN <- cbind(x[["loadingMN"]][, 1], x[["orthoLoadingMN"]][, parCompVi[2] - 1])
            colnames(pCompMN) <- colnames(tCompMN) <- c("h1", paste("o", parCompVi[2] - 1, sep = ""))
        } else {
            if(max(parCompVi) > x[["summaryDF"]][, "pre"])
                stop("Selected component for plotting as ordinate exceeds the total number of predictive components of the model", call. = FALSE)
            tCompMN <- x[["scoreMN"]][, parCompVi, drop = FALSE]
            pCompMN <- x[["loadingMN"]][, parCompVi, drop = FALSE]
        }

    }

    ## if(ncol(tCompMN) > 1) {

    ##     mahInvCovMN <- solve(cov(tCompMN))

    ##     pcaResMN <- cbind(sdsVn = apply(tCompMN,
    ##                           1,
    ##                           function(x) sqrt(t(as.matrix(x)) %*% mahInvCovMN %*% as.matrix(x))),
    ##                       odsVn = apply(x[["suppLs"]][["xModelMN"]] - tcrossprod(tCompMN, pCompMN),
    ##                           1,
    ##                           function(x) sqrt(drop(crossprod(x[complete.cases(x)])))))

    ## } else
    ##     pcaResMN <- NULL

    cxtCompMN <- cor(x[["suppLs"]][["xModelMN"]], tCompMN,
                     use = "pairwise.complete.obs")

    if(!is.null(x[["suppLs"]][["yModelMN"]]))
        cytCompMN <- cor(x[["suppLs"]][["yModelMN"]], tCompMN, use = "pairwise.complete.obs")


    if(x[["topLoadI"]] * 4 < ncol(x[["suppLs"]][["xModelMN"]])) {

        pexVi <- integer(x[["topLoadI"]] * ncol(pCompMN) * 2) ## 'ex'treme values

        for(k in 1:ncol(pCompMN)) {

            pkVn <-  pCompMN[, k]

            pexVi[1:(2 * x[["topLoadI"]]) + 2 * x[["topLoadI"]] * (k - 1)] <- c(order(pkVn)[1:x[["topLoadI"]]],
                                                                         rev(order(pkVn, decreasing = TRUE)[1:x[["topLoadI"]]]))

        }

    } else
        pexVi <- 1:ncol(x[["suppLs"]][["xModelMN"]])

    pxtCompMN <- cbind(pCompMN,
                       cxtCompMN)

    if(ncol(pCompMN) == 1) {
       colnames(pxtCompMN)[2] <- paste0("cor_", colnames(pxtCompMN)[2])
    } else
        colnames(pxtCompMN)[3:4] <- paste0("cor_", colnames(pxtCompMN)[3:4])

    topLoadMN <- pxtCompMN

    topLoadMN <- topLoadMN[pexVi, , drop = FALSE]

    if(x[["topLoadI"]] * 4 < ncol(x[["suppLs"]][["xModelMN"]]) &&
       ncol(pCompMN) > 1) {

        topLoadMN[(2 * x[["topLoadI"]] + 1):(4 * x[["topLoadI"]]), c(1, 3)] <- NA
        topLoadMN[1:(2 * x[["topLoadI"]]), c(2, 4)] <- NA

    }


    ## Observation and variable names and colors
    ##------------------------------------------

    ## obsLabVc

    if(!any(is.na(parLabVc))) {
        obsLabVc <- parLabVc
    } else if(!is.null(rownames(tCompMN))) {
        obsLabVc <- rownames(tCompMN)
    } else
        obsLabVc <- as.character(1:nrow(tCompMN))

    if(!is.null(x[["subset"]])) {
        obsLabVc <- obsLabVc[x[["subset"]]]
        tesLabVc <- obsLabVc[-x[["subset"]]]
    } else
        tesLabVc <- ""

    ## obsColVc

    if(!any(is.na(parAsColFcVn))) {
        obsColVc <- .colorF(as.vector(parAsColFcVn))[["colVc"]]
        obsLegVc <- as.vector(parAsColFcVn)
    } else if(!is.null(x[["suppLs"]][["yMCN"]]) && ncol(x[["suppLs"]][["yMCN"]]) == 1) {
        obsColVc <- .colorF(c(x[["suppLs"]][["yMCN"]]))[["colVc"]]
        obsLegVc <- c(x[["suppLs"]][["yMCN"]])
    } else {
        obsColVc <- rep("black", nrow(tCompMN))
        obsLegVc <- NULL
    }

    if(!is.null(x[["subset"]])) {
        obsColVc <- obsColVc[x[["subset"]]]
        tesColVc <- obsColVc[-x[["subset"]]]
        if(!is.null(obsLegVc)) {
            obsLegVc <- obsLegVc[x[["subset"]]]
            tesLegVc <- obsLegVc[-x[["subset"]]]
        }
    }


    ## Layout
    ##-------

    if(!parDevNewL && length(plotVc) != 1)
        stop("'plotVc' must be of length 1 when 'parDevNewL' is set to FALSE", call. = FALSE)

    if(parDevNewL) {
        layRowN <- ceiling(sqrt(length(plotVc)))
        if(is.null(file.pdfC))
            dev.new()
        else
            pdf(file.pdfC)
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
               opLs = x,
               obsColVc = obsColVc,
               obsLabVc = obsLabVc,
               obsLegVc = obsLegVc,
               layL = layL,
               parCexN = parCexN,
               parEllipsesL = parEllipsesL,
               parTitleL = parTitleL,
               parCompVi = parCompVi,
               plotVc = plotVc,
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


} ## plot.opls


.plotF <- function(ploC,
                   opLs,
                   obsColVc,
                   obsLabVc,
                   obsLegVc,
                   layL,
                   parCexN,
                   parEllipsesL,
                   parTitleL,
                   parCompVi,
                   plotVc,
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

            if(opLs[["summaryDF"]][, "ort"] > 0)
                yLabC <- paste("with tOrtho",
                               parCompVi[2] - 1,
                               sep = "")

            yLimVn <- xLimVn <- c(-1, 1)

            ploMN <- cxtCompMN

            if(opLs[["typeC"]] != "PCA")
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

            ypMN <- eval(parse(text = paste("opLs[['suppLs']][['y", switch(unlist(strsplit(ploC, "-"))[2], train = "Pre", test = "Tes"), "MN']]", sep = ""))) ## predicted

            if(is.null(opLs[["subset"]]))
                yaMCN <- opLs[["suppLs"]][["yMCN"]] ## actual
            else {
                if(grepl("train", ploC))
                    yaMCN <- opLs[["suppLs"]][["yMCN"]][opLs[["subset"]], , drop = FALSE]
                else
                    yaMCN <- opLs[["suppLs"]][["yMCN"]][-opLs[["subset"]], , drop = FALSE]
            }

            if(mode(opLs[["suppLs"]][["yMCN"]]) == "character") {
                yaMN <- opLs[["suppLs"]][[".char2numF"]](yaMCN)
            } else
                yaMN <- yaMCN

            ploMN <- cbind(ypMN,
                           yaMN) ## to be modified (when ncol(yPreMCN) > 1)

        } else if(ploC == "x-loading") {

            maiC <- "Loadings"

            xLabC <- paste0("p",
                            parCompVi[1],
                            " (",
                            round(opLs[["modelDF"]][parCompVi[1], "R2X"] * 100),
                            "%)")

            yLabC <- paste0("p",
                            parCompVi[2],
                            " (",
                            round(opLs[["modelDF"]][parCompVi[2], "R2X"] * 100),
                            "%)")

            ploMN <- pCompMN

            if(!is.null(opLs[["suppLs"]][["yMCN"]]) && opLs[["summaryDF"]][, "ort"] > 0)
                yLabC <- paste0("pOrtho",
                                parCompVi[2] - 1,
                                " (",
                                round(opLs[["modelDF"]][parCompVi[2] - 1, "R2X"] * 100),
                                "%)")


        } else if(ploC == "x-score") {

            maiC <- paste0("Scores (", opLs[["typeC"]], ")")

            xLabC <- paste0("t",
                            parCompVi[1],
                            " (",
                            round(opLs[["modelDF"]][parCompVi[1], "R2X"] * 100),
                            "%)")

            yLabC <- paste0("t",
                            parCompVi[2],
                            " (",
                            round(opLs[["modelDF"]][parCompVi[2], "R2X"] * 100),
                           "%)")

            ploMN <- tCompMN

            if(grepl("^OPLS", opLs[["typeC"]]))
                yLabC <- paste0("to", parCompVi[2] - 1)

            xLimVn <- c(-1, 1) * max(sqrt(var(ploMN[, 1]) * hotFisN), max(abs(ploMN[, 1])))
            yLimVn <- c(-1, 1) *max(sqrt(var(ploMN[, 2]) * hotFisN), max(abs(ploMN[, 2])))

            ploColVc <- obsColVc

        } else if(ploC == "xy-score") {

            maiC <- "XY-Scores"
            xLabC <- paste("t", parCompVi[1], sep = "")
            yLabC <- paste("u/c", parCompVi[1], sep = "")

            ploMN <- cbind(opLs[["scoreMN"]][, parCompVi[1]], opLs[["uMN"]][, parCompVi[1]] / opLs[["cMN"]][parCompVi[1]])

            ploColVc <- obsColVc

        } else if(ploC == "xy-weight") {

            maiC <- "Weights"
            xLabC <- paste0("w*c", parCompVi[1])
            yLabC <- paste0("w*c", parCompVi[2])

            ploMN <- rbind(opLs[["rotationMN"]][, parCompVi],
                           opLs[["cMN"]][, parCompVi])

            pchVn <- rep(17, times = nrow(ploMN))
            ploColVc <- rep("black", times = nrow(ploMN))

            pchVn[(nrow(opLs[["rotationMN"]]) + 1):nrow(ploMN)] <- 15
            ploColVc[(nrow(opLs[["rotationMN"]]) + 1):nrow(ploMN)] <- "red"

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
            ## corPchVn <- rep(18, ncol(opLs[["suppLs"]][["xModelMN"]]))
            ## corNamVc <- colnames(opLs[["suppLs"]][["xModelMN"]])


            if(opLs[["typeC"]] != "PCA") {
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

            if(opLs[["typeC"]] != "PCA" && length(plotVc) == 1)
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

            points(ploMN[pexVi, ],
                   pch = 18,
                   col = "red")

            ## pexLabVc <- colnames(opLs[["suppLs"]][["xModelMN"]])[pexVi]
            pexLabVc <- rownames(opLs[["loadingMN"]])[pexVi]
            pexLabVc[duplicated(pexLabVc)] <- ""

            text(ploMN[pexVi, ],
                 cex = parCexN,
                 col = "red",
                 labels = pexLabVc,
                 pos = rep(c(4, 2, 3, 1), each = opLs[["topLoadI"]]))

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

            mtext(paste("R2X", round(opLs[["summaryDF"]][, "R2X(cum)"], 3), sep = "\n"),
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

            if(!is.null(opLs[["suppLs"]][["yMCN"]])) {

                mtext(paste("R2Y", round(opLs[["summaryDF"]][, "R2Y(cum)"], 3), sep = "\n"),
                      at = pu1N * ifelse(layL, 1, 0.8),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("Q2Y", round(opLs[["summaryDF"]][, "Q2(cum)"], 3), sep = "\n"),
                      at = pu1N * ifelse(layL, 0.6, 0.4),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("RMSEE", round(opLs[["summaryDF"]][, "RMSEE"], 3), sep = "\n"),
                      at =  -pu1N * ifelse(layL, 0.6, 0.4),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                mtext(paste("pre", opLs[["summaryDF"]][, "pre"], sep = "\n"),
                      at = -pu1N * ifelse(layL, 0.92, 0.7),
                      cex = cexRqcN,
                      font = 1,
                      line = 3,
                      side = 1)

                if(opLs[["summaryDF"]][, "ort"] > 0)
                    mtext(paste("ort", opLs[["summaryDF"]][, "ort"], sep = "\n"),
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
                 labels = c(rownames(opLs[["rotationMN"]]), rownames(opLs[["cMN"]])))

            if(!layL)
                legend("topleft",
                       col = c("black", "red"),
                       legend = c("X", "Y"),
                       text.col = c("black", "red"))

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
                          odsVn = apply(opLs[["suppLs"]][["xModelMN"]] - tcrossprod(tCompMN, pCompMN),
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

        if(opLs[["typeC"]] == "PCA") {

            barplot(opLs[["modelDF"]][, "R2X"] * 100,
                    main = "Variance explained",
                    names.arg = rownames(opLs[["modelDF"]]),
                    xlab = "PC",
                    ylab = "% of total variance")


        } else {

            if(opLs[["summaryDF"]][, "ort"] == 0) {
                modBarDF <- opLs[["modelDF"]]
            } else
                modBarDF <- opLs[["modelDF"]][!(rownames(opLs[["modelDF"]]) %in% c("rot", "sum")), ]

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

        plot(c(min(opLs[["suppLs"]][["permMN"]][, "sim"]), 1),
             c(min(opLs[["suppLs"]][["permMN"]][, c("R2Y(cum)", "Q2(cum)")]), 1),
             main = paste0("pR2Y = ",
                 opLs[["summaryDF"]][, "pR2Y"],
                 ", pQ2 = ",
                 opLs[["summaryDF"]][, "pQ2"]),
             type = "n",
             xlab=expression(Similarity(bold(y), bold(y[perm]))),
             ylab = "")

        points(opLs[["suppLs"]][["permMN"]][, "sim"], opLs[["suppLs"]][["permMN"]][, "Q2(cum)"],
               col = "black",
               pch = 18)
        abline(h = opLs[["suppLs"]][["permMN"]][1, "Q2(cum)"],
               col = "black")
        points(opLs[["suppLs"]][["permMN"]][, "sim"], opLs[["suppLs"]][["permMN"]][, "R2Y(cum)"],
               col = "red",
               pch = 18)
        abline(h = opLs[["suppLs"]][["permMN"]][1, "R2Y(cum)"],
               col = "red")
        .legendF(c("R2Y", "Q2Y"),
                 "bottomright")

    } ## permutation


    if(ploC == "x-variance") {

        par(las=2)

        boxplot(opLs[["suppLs"]][["xSubIncVarMN"]],
                main = "X variances (min, med, max)",
                names=rep("", 3),
                xaxt="n",
                yaxt="n")

        axis(at=axTicks(2),
             side=2)

        mtext(substr(colnames(opLs[["suppLs"]][["xSubIncVarMN"]]), 1, 9),
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

} ## .legendF


.similarityF <- function(x, y,
                         .char2numF,
                         charL = FALSE) {

    if(charL) {
        x <- .char2numF(x)
        y <- .char2numF(y)
    }

    return(cor(x, y, use = "pairwise.complete.obs"))

} ## .similarityF


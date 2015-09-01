print.opls <- function(x, .sinkC = NULL, ...) {

    message(x[["typeC"]])

    message(x[["descriptionMC"]]["samples", ],
            " samples x ",
            x[["descriptionMC"]]["X_variables", ],
            " variables",
            ifelse(!is.null(x[["yPredMCN"]]),
                   paste0(" and ", ncol(x[["yPredMCN"]]), " response", ifelse(ncol(x[["yPredMCN"]]) > 1, "s", "")),
                   ""))

    message(x[["descriptionMC"]]["missing_values", ], " NAs")

    message(x[["descriptionMC"]]["near_zero_excluded_X_variables", ],
            " excluded variables (near zero variance)")

    message(x[["suppLs"]][["scaleC"]], " x", ifelse(x[["typeC"]] != "PCA", " and y", ""), " scaling")

    optDigN <- options()[["digits"]]
    options(digits = 3)
    print(x[["summaryDF"]])
    options(digits = optDigN)

} ## print.opls


summary.opls <- function(object, ...)
    structure(object, class="summary.opls")


print.summary.opls <- function (x, ...) {

    message("\n1) Data set:\n")

    message(x[["descriptionMC"]]["samples", ],
            " samples x ",
            x[["descriptionMC"]]["X_variables", ],
            " variables",
            ifelse(!is.null(x[["yPredMCN"]]),
                   paste0(" and ", ncol(x[["yPredMCN"]]), " response", ifelse(ncol(x[["yPredMCN"]]) > 1, "s", "")),
                   ""))

    message(x[["descriptionMC"]]["missing_values", ], " NAs")

    message(x[["descriptionMC"]]["near_zero_excluded_X_variables", ], " excluded variables (near zero variance)")

    message(x[["suppLs"]][["scaleC"]], " x", ifelse(x[["typeC"]] != "PCA", " and y", ""), " scaling")

    message("Summary of the ", x[["topLoadI"]], " increasing variance spaced raw variables:")
    print(summary(x[["suppLs"]][["xSubIncVarMN"]]))

    if(!is.null(x[["suppLs"]][["xCorMN"]])) {
        message("Correlations between the X-variables:")
        print(signif(x[["suppLs"]][["xCorMN"]], 2))
    }

    message("\n2) Model: ", x[["typeC"]], "\n")

    message("Correlations between variables and first 2 components:")

    if(x[["summaryDF"]][, "pre"] + x[["summaryDF"]][, "ort"] < 2) {

        warning("A single component model has been selected by cross-validation", call. = FALSE)

        tCompMN <- x[["scoreMN"]]
        pCompMN <- x[["loadingMN"]]

    } else {

        if(x[["summaryDF"]][, "ort"] > 0) {
            tCompMN <- cbind(x[["scoreMN"]][, 1], x[["orthoScoreMN"]][, 1])
            pCompMN <- cbind(x[["loadingMN"]][, 1], x[["orthoLoadingMN"]][, 1])
            colnames(pCompMN) <- colnames(tCompMN) <- c("h1", "o1")
        } else {
            tCompMN <- x[["scoreMN"]][, 1:2]
            pCompMN <- x[["loadingMN"]][, 1:2]
        }

    }

    cxtCompMN <- cor(x[["suppLs"]][["xModelMN"]], tCompMN,
                     use = "pairwise.complete.obs")

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
    print(signif(topLoadMN, 2))

    message("")

    optDigN <- options()[["digits"]]
    options(digits = 3)
    print(x[["modelDF"]])
    options(digits = optDigN)

}



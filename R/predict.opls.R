predict.opls <- function(object, newdata, ...) {

    if(object[["typeC"]] == "PCA")
        stop("Predictions currently available for (O)PLS(-DA) models only (not PCA)",
             call. = FALSE)

    if(missing(newdata)) {

        return(object[["fitted"]])

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

        if(ncol(newdata) != as.numeric(object[["descriptionMC"]]["X_variables", ])) {
            if(length(object[["xZeroVarVi"]]) == 0) {
                stop("'newdata' number of variables is ",
                     ncol(newdata),
                     " whereas the number of variables used for model training was ",
                     as.numeric(object[["descriptionMC"]]["X_variables", ]),
                     ".",
                     call. = FALSE)
            } else if(ncol(newdata) - as.numeric(object[["descriptionMC"]]["X_variables", ]) ==
               as.numeric(object[["descriptionMC"]]["near_zero_excluded_X_variables", ])) {
                warning(as.numeric(object[["descriptionMC"]]["near_zero_excluded_X_variables", ]),
                        " near zero variance variables excluded during the model training will be removed from 'newdata'.",
                        call. = FALSE)
                newdata <- newdata[, -object[["xZeroVarVi"]], drop = FALSE]
            } else {
                stop("'newdata' number of variables (",
                     ncol(newdata),
                     ") does not correspond to the number of initial variables (",
                     as.numeric(object[["descriptionMC"]]["X_variables", ]),
                     ") minus the number of near zero variance variables excluded during the training (",
                     as.numeric(object[["descriptionMC"]]["near_zero_excluded_X_variables", ]),
                     ").",
                     call. = FALSE)
            }
        }

        xteMN <- scale(newdata, object[["xMeanVn"]], object[["xSdVn"]])

        if(object[["summaryDF"]][, "ort"] > 0) {

            for(noN in 1:object[["summaryDF"]][, "ort"]) {
                if(object[["suppLs"]][["naxL"]]) {
                    xtoMN <- matrix(0, nrow = nrow(xteMN), ncol = 1)
                    for(i in 1:nrow(xtoMN)) {
                        comVl <- complete.cases(xteMN[i, ])
                        xtoMN[i, ] <- crossprod(xteMN[i, comVl], object[["orthoWeightMN"]][comVl, noN]) / drop(crossprod(object[["orthoWeightMN"]][comVl, noN]))
                    }
                } else
                    xtoMN <- xteMN %*% object[["orthoWeightMN"]][, noN]

                xteMN <- xteMN - tcrossprod(xtoMN, object[["orthoLoadingMN"]][, noN])
            }

        }

        if(object[["suppLs"]][["naxL"]]) {
            yTesScaMN <- matrix(0, nrow = nrow(xteMN), ncol = ncol(object[["coefficients"]]),
                                dimnames = list(rownames(xteMN), colnames(object[["coefficients"]])))
            for(j in 1:ncol(yTesScaMN))
                for(i in 1:nrow(yTesScaMN)) {
                    comVl <- complete.cases(xteMN[i, ])
                    yTesScaMN[i, j] <- crossprod(xteMN[i, comVl], object[["coefficients"]][comVl, j])
                }
        } else
            yTesScaMN <- xteMN %*% object[["coefficients"]]

        ## if(object[["suppLs"]][["nayL"]])
        ##     yTesScaMN <- yTesScaMN[!is.na(yMCN[testVi, ]), , drop = FALSE]

        yTesMN <- scale(scale(yTesScaMN,
                              FALSE,
                              1 / object[["ySdVn"]]),
                        -object[["yMeanVn"]],
                        FALSE)
        attr(yTesMN, "scaled:center") <- NULL
        attr(yTesMN, "scaled:scale") <- NULL

        if(is.factor(object[["fitted"]])) {

            yTestMCN <- object[["suppLs"]][[".char2numF"]](yTesMN,
                                                           c2nL = FALSE)
            predMCNFcVcn <- as.character(yTestMCN)
            names(predMCNFcVcn) <- rownames(newdata)
            predMCNFcVcn <- factor(predMCNFcVcn, levels = object[["suppLs"]][["yLevelVc"]])

        } else if(is.vector(object[["fitted"]])) {

            if(is.character(object[["fitted"]])) {

                yTestMCN <- object[["suppLs"]][[".char2numF"]](yTesMN,
                                                               c2nL = FALSE)
                predMCNFcVcn <- as.character(yTestMCN)
                names(predMCNFcVcn) <- rownames(newdata)

            } else {

                predMCNFcVcn <- as.numeric(yTesMN)
                names(predMCNFcVcn) <- rownames(newdata)

            }

        } else if(is.matrix(object[["fitted"]])) {

            if(mode(object[["fitted"]]) == "character") {
               predMCNFcVcn  <- object[["suppLs"]][[".char2numF"]](yTesMN,
                                                               c2nL = FALSE)
            } else
                predMCNFcVcn <- yTesMN

            dimnames(predMCNFcVcn) <- list(rownames(newdata), "pred")

        }

        return(predMCNFcVcn)

    }

} ## predict.opls

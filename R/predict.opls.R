predict.opls <- function(object, newdata, ...) {

    ## TO DO: PCA prediction

    if(missing(newdata)) {

        return(object[["fitted"]])

    } else {

        if(class(newdata) != "matrix" || mode(newdata) != "numeric")
            stop("'newdata' must be a numerical matrix")

        if(length(object[["xZeroVarVi"]]))
            newdata <- newdata[, -object[["xZeroVarVi"]]]

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
            predMNFcVcn <- as.character(yTestMCN)
            names(predMNFcVcn) <- rownames(newdata)
            predMNFcVcn <- factor(predMNFcVcn, levels = object[["suppLs"]][["yLevelVc"]])

        } else if(is.vector(object[["fitted"]])) {

            predMNFcVcn <- as.vector(yTestMCN)
            names(predMNFcVcn) <- rownames(newdata)

        } else if(is.matrix(object[["fitted"]])) {

            predMNFcVcn <- yTestMCN
            rownames(predMNFcVcn) <- rownames(newdata)

        } else
            stop() ## should not occur


        return(predMNFcVcn)

    }

} ## predict.opls

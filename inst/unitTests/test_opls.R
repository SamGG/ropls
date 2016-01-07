test_str <- function() {

    data(foods)

    strF(foods)
    strF(foods, border = 3)

    fatML <- matrix(TRUE, nrow = 1, ncol = 1000)
    strF(fatML, bigMarkC = "'")

    testMC <- matrix("a", nrow = 10, ncol = 10)
    strF(testMC)

    data(sacurine)

    strF(sacurine[["dataMatrix"]])

    strF(sacurine[["sampleMetadata"]])

}

test_plot <- function() {

    data(sacurine)

    for(typeC in c("correlation", "outlier", "overview",
               "permutation", "predict-train","predict-test",
               "summary", "x-loading", "x-score", "x-variance",
               "xy-score", "xy-weight")) {

        if(grepl("predict", typeC))
            subset <- "odd"
        else
            subset <- NULL

        opLs <- opls(sacurine[["dataMatrix"]],
                     sacurine[["sampleMetadata"]][, "gender"],
                     predI = ifelse(typeC != "xy-weight", 1, 2),
                     orthoI = ifelse(typeC != "xy-weight", 1, 0),
                     permI = ifelse(typeC == "permutation", 10, 0),
                     subset = subset,
                     printL = FALSE, plotL = FALSE)

        plot(opLs, typeVc = typeC)

    }

}

test_print <- function() {

    data(sacurine)
    pcaLs <- opls(sacurine[["dataMatrix"]], predI = 2, printL = FALSE, plotL = FALSE)
    print(pcaLs)
    plsLs <- opls(sacurine[["dataMatrix"]],
                  sacurine[["sampleMetadata"]][, "gender"],
                  predI = 2, printL = FALSE, plotL = FALSE)
    print(plsLs)


}

test_PCA <- function() {

    data(foods)

    fooMN <- as.matrix(foods[, colnames(foods) != "Country"])
    rownames(fooMN) <- foods[, "Country"]
    foo.pca <- opls(fooMN, permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(foo.pca[["pcaVarVn"]][1],
                       6.19,
                       tolerance = 1e-3)
    checkEqualsNumeric(foo.pca[["scoreMN"]][1, 1],
                       1.309494,
                       tolerance = 1e-5)

}

test_PLS_single <- function() {

    data(cornell) ## see Tenenhaus, 1998

    corPlsLs <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                     matrix(cornell[, "y"], ncol = 1),
                     permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(corPlsLs[["scoreMN"]][1, 1],
                       2.0513,
                       tolerance = 1e-3)
    cornell.pls <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                     cornell[, "y"],
                     permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(cornell.pls[["vipVn"]][1],
                       1.125,
                       tolerance = 1e-3)

}

test_PLS_multiple <- function() {

    data(lowarp) ## see Eriksson et al. (2001); presence of NAs

    lowarp.pls <- opls(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
                       as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
                                        grepl("^st", colnames(lowarp))]),
                       permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(lowarp.pls[["weightStarMN"]][2, 1],
                       -0.0861,
                       tolerance = 1e-3)
    checkEqualsNumeric(lowarp.pls[["uMN"]]["s5", "h2"],
                       0.596598855,
                       tolerance = 1e-8)

}

test_sacurine_PCA <- function() {

    data(sacurine)

    sac.pca <- opls(sacurine[["dataMatrix"]],
                    printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.pca[["summaryDF"]]["h8", "R2X(cum)"],
                       0.501,
                       tolerance = 1e-3)
    sac.pcaSvd <- opls(sacurine[["dataMatrix"]],
                       algoC = "svd",
                       printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.pcaSvd[["scoreMN"]][1, 1],
                       -8.744009,
                       tolerance = 1e-5)

}

test_sacurine_PLSDA <- function() {

    data(sacurine)

    sac.plsda <- opls(sacurine[["dataMatrix"]],
                      matrix(sacurine[["sampleMetadata"]][, "gender"], ncol = 1),
                      permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.plsda[["summaryDF"]]["h3", "Q2(cum)"],
                       0.5835585,
                       tolerance = 1e-7)
    sac.plsda <- opls(sacurine[["dataMatrix"]],
                      sacurine[["sampleMetadata"]][, "gender"],
                      permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.plsda[["vipVn"]][5],
                       0.7034828,
                       tolerance = 1e-7)

    ## splitting the dataset between a reference and a test subsets

    sac.plsdaCrv <- opls(sacurine[["dataMatrix"]],
                         sacurine[["sampleMetadata"]][, "gender"],
                         subset = setdiff(1:nrow(sacurine[["dataMatrix"]]), 1:10),
                         permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.plsdaCrv[["summaryDF"]]["h3", "RMSEP"],
                       0.2745173,
                       tolerance = 1e-6)
    ## for subset = "odd"
    ##     R2X(cum)  R2Y(cum)   Q2(cum)     RMSEE     RMSEP pre ort
    ## h2 0.1802609 0.7666581 0.5618918 0.2446342 0.3422395   2   0

}

test_sacurine_OPLSDA <- function() {

    data(sacurine)

    sac.oplsda <- opls(sacurine[["dataMatrix"]],
                       matrix(sacurine[["sampleMetadata"]][, "gender"],
                              ncol = 1),
                       predI = 1,
                       orthoI = 1,
                       permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.oplsda[["orthoScoreMN"]][1, 1],
                       4.980604,
                       tolerance = 1e-7)
    ##     R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort
    ## sum    0.185    0.668   0.554 0.289   1   1

    ## permutation testing

    sac.oplsdaPer <- opls(sacurine[["dataMatrix"]],
                          sacurine[["sampleMetadata"]][, "gender"],
                          predI = 1,
                          orthoI = 1,
                          permI = 10,
                          printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.oplsdaPer[["suppLs"]][["permMN"]][1, 1],
                       0.1845708,
                       tolerance = 1e-6)

}


test_sacurine_PLSDA_pareto <- function() {

    data(sacurine)

    sac.plsdaPar <- opls(sacurine[["dataMatrix"]],
                         matrix(sacurine[["sampleMetadata"]][, "gender"], ncol = 1),
                         scaleC = "pareto",
                         permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.plsdaPar[["summaryDF"]]["h2", "Q2(cum)"],
                       0.5210747,
                       tolerance = 1e-7)
    sac.plsdaPar <- opls(sacurine[["dataMatrix"]],
                         sacurine[["sampleMetadata"]][, "gender"],
                         scaleC = "pareto",
                         permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.plsdaPar[["vipVn"]][5],
                       0.7710519,
                       tolerance = 1e-7)

}

test_sacurine_PLS_predict <- function() {

    data(sacurine)

    sac.opls <- opls(sacurine[["dataMatrix"]],
                     sacurine[["sampleMetadata"]][, "age"],
                     predI = 1,
                     orthoI = 1,
                     permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(predict(sac.opls)[107],
                       42.97439,
                       tolerance = 1e-5)

    sac.opls <- opls(sacurine[["dataMatrix"]],
                     sacurine[["sampleMetadata"]][, "gender"],
                     predI = 1,
                     orthoI = 1,
                     permI = 0, printL = FALSE, plotL = FALSE)
    pred107Fc <- "M"
    names(pred107Fc) <- "HU_125"
    pred107Fc <- factor(pred107Fc, levels = c("M", "F"))
    checkEquals(predict(sac.opls)[107],
                pred107Fc)

    trainVi <- 1:floor(nrow(sacurine[["dataMatrix"]]) / 2)

    sac.pls.tr <- opls(sacurine[["dataMatrix"]],
                       sacurine[["sampleMetadata"]][, "age"],
                       subset = trainVi,
                       permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(as.numeric(sac.pls.tr[["descriptionMC"]]["samples", 1]),
                       length(trainVi),
                       tolerance = 0)
    checkEqualsNumeric(predict(sac.pls.tr)[90],
                       42.94277,
                       tolerance = 1e-5)
    testVi <- setdiff(1:nrow(sacurine[["dataMatrix"]]), trainVi)
    predVn <- predict(sac.pls.tr, sacurine[["dataMatrix"]][testVi, ])
    checkEqualsNumeric(predVn["HU_207"],
                       36.38991,
                       tolerance = 1e-5)

}

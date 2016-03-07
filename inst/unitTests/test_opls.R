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

    ## (O)PLS-DA: Turning display of ellipses off in case parAsColFcVn is numeric
    plot(opLs,
         parAsColFcVn = sacurine[["sampleMetadata"]][, "age"])

    ## Converting 'parAsColFcVn' character into a factor (with warning)
    plot(opLs,
         parAsColFcVn = as.character(sacurine[["sampleMetadata"]][, "gender"]))

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

    data(foods) ## contains 3 NA

    fooMN <- as.matrix(foods[, colnames(foods) != "Country"])
    rownames(fooMN) <- foods[, "Country"]

    foo.pca <- opls(fooMN, permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getPcaVarVn(foo.pca)[1],
                       6.19,
                       tolerance = 1e-3)
    checkEqualsNumeric(getScoreMN(foo.pca)[1, 1],
                       1.309494,
                       tolerance = 1e-5)

}

test_PLS_single <- function() {

    data(cornell) ## see Tenenhaus, 1998

    corPlsLs <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                     matrix(cornell[, "y"], ncol = 1),
                     permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getScoreMN(corPlsLs)[1, 1],
                       2.0513,
                       tolerance = 1e-3)
    cornell.pls <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                        cornell[, "y"],
                        permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getVipVn(cornell.pls)[1],
                       1.125,
                       tolerance = 1e-3)

}

test_PLS_multiple <- function() {

    data(lowarp) ## see Eriksson et al. (2001); presence of NAs

    lowarp.pls <- opls(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
                       as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
                                        grepl("^st", colnames(lowarp))]),
                       permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(lowarp.pls@weightStarMN[2, 1],
                       -0.0861,
                       tolerance = 1e-3)
    checkEqualsNumeric(lowarp.pls@uMN["s5", "p2"],
                       0.596598855,
                       tolerance = 1e-8)

}

test_sacurine_PCA <- function() {

    data(sacurine)

    sac.pca <- opls(sacurine[["dataMatrix"]],
                    printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getSummaryDF(sac.pca)["Total", "R2X(cum)"],
                       0.501,
                       tolerance = 1e-3)
    sac.pcaSvd <- opls(sacurine[["dataMatrix"]],
                       algoC = "svd",
                       printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getScoreMN(sac.pcaSvd)[1, 1],
                       -8.744009,
                       tolerance = 1e-5)

}

test_sacurine_PLSDA <- function() {

    data(sacurine)

    sac.plsda <- opls(sacurine[["dataMatrix"]],
                      matrix(sacurine[["sampleMetadata"]][, "gender"], ncol = 1),
                      permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getSummaryDF(sac.plsda)["Total", "Q2(cum)"],
                       0.584,
                       tolerance = 1e-3)
    sac.plsda <- opls(sacurine[["dataMatrix"]],
                      sacurine[["sampleMetadata"]][, "gender"],
                      permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getVipVn(sac.plsda)[5],
                       0.7034828,
                       tolerance = 1e-7)

    ## splitting the dataset between a reference and a test subsets

    sac.plsdaCrv <- opls(sacurine[["dataMatrix"]],
                         sacurine[["sampleMetadata"]][, "gender"],
                         subset = setdiff(1:nrow(sacurine[["dataMatrix"]]), 1:10),
                         permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getSummaryDF(sac.plsdaCrv)["Total", "RMSEP"],
                       0.275,
                       tolerance = 1e-3)
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
    checkEqualsNumeric(getScoreMN(sac.oplsda, orthoL = TRUE)[1, 1],
                       4.980604,
                       tolerance = 1e-7)
    ##     R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort
    ## sum    0.185    0.668   0.554 0.289   1   1

    checkEqualsNumeric(getVipVn(sac.oplsda)[2],
                       1.551154,
                       tolerance = 1e-6)

    checkEqualsNumeric(getVipVn(sac.oplsda, orthoL = TRUE)[3],
                       1.642173,
                       tolerance = 1e-6)

    ## permutation testing

    sac.oplsdaPer <- opls(sacurine[["dataMatrix"]],
                          sacurine[["sampleMetadata"]][, "gender"],
                          predI = 1,
                          orthoI = 1,
                          permI = 10,
                          printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.oplsdaPer@suppLs[["permMN"]][1, 1],
                       0.185,
                       tolerance = 1e-3)

}


test_sacurine_PLSDA_pareto <- function() {

    data(sacurine)

    sac.plsdaPar <- opls(sacurine[["dataMatrix"]],
                         matrix(sacurine[["sampleMetadata"]][, "gender"], ncol = 1),
                         scaleC = "pareto",
                         permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getSummaryDF(sac.plsdaPar)["Total", "Q2(cum)"],
                       0.521,
                       tolerance = 1e-3)
    sac.plsdaPar <- opls(sacurine[["dataMatrix"]],
                         sacurine[["sampleMetadata"]][, "gender"],
                         scaleC = "pareto",
                         permI = 0, printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(getVipVn(sac.plsdaPar)[5],
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
    checkEqualsNumeric(as.numeric(sac.pls.tr@descriptionMC["samples", 1]),
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

test_subset <- function() {
    ## bug fixed following T. Souza remark

    tiagoMN <- c(10, 11, 12, 11, 10, 11, 10, 11, 10, 11, 10, 11, 10, 11, 10, 11, 10, 11, 10, 8, 10, 11, 10, 12, 10, 11, 10, 11, 10, 11, 11, 10, 11, 10, 11, 10, 11, 10, 12, 10, 5, 6, 5, 6, 7, 6, 5, 6, 5, 6, 15, 16, 15, 16, 15, 16, 15, 14, 15, 15, 16, 16, 16, 14, 16, 15, 16, 14, 16, 14, 2, 3, 2, 3, 22, 3, 4, 3, 2, 3, 3, 7, 3, 7, 3, 7, 4, 7, 3, 7, 3, 7, 3, 7, 3, 4, 3, 7, 3, 7, 2, 8, 2, 8, 20, 8, 2, 4, 2, 8, 7, 3, 7, 3, 4, 3, 7, 3, 7, 7, 3, 7, 3, 7, 3, 7, 3, 7, 6, 8, 25, 29, 25, 24, 25, 29, 25, 29, 25, 25, 23, 27, 23, 27, 23, 27, 23, 27, 23, 25, 20, 21, 22, 21, 20, 21, 20, 21, 20, 21)
    dim(tiagoMN) <- c(20, 8)
    dimnames(tiagoMN) <- list(1:20, letters[1:8])
    tiagoFc <- rep(c("P", "F"), each = 10)

    opls(tiagoMN, tiagoFc, predI = 1, orthoI = 1, subset = "odd")

    oplExc <- opls(tiagoMN, tiagoFc, subset = seq(1, 20, by = 2))
    checkException(plot(oplExc, parLabVc = paste0("s", seq(1, 20, by = 2))),
                   silent = TRUE)

}

test_plsda_multiclass <- function() {

    data(sacurine)
    ageVn <- sacurine[["sampleMetadata"]][, "age"]
    ageVc <- ifelse(ageVn < 25, "young",
                    ifelse(ageVn < 40, NA,
                           ifelse(ageVn < 45, "normal",
                                  ifelse(ageVn < 57, NA,
                                         "old"))))

    sacMN <- sacurine[["dataMatrix"]][!is.na(ageVc), ]
    ageVc <- ageVc[!is.na(ageVc)]

    oplsda.mul <- opls(sacMN, ageVc, predI = 2)

    checkEqualsNumeric(getSummaryDF(oplsda.mul)[, "Q2(cum)"],
                       0.0394,
                       tolerance = 1e-3)

    checkEqualsNumeric(length(getVipVn(oplsda.mul, orthoL = TRUE)),
                       0,
                       tolerance = 1e-10)

    checkException(opls(sacMN, ageVc, orthoI = NA),
                   silent = TRUE)

    oplsda.mul.par <- opls(sacMN, ageVc,
                           predI = 2,
                           scaleC = "pareto")

    checkEqualsNumeric(oplsda.mul.par@modelDF["p2", "R2Y(cum)"],
                       0.442,
                       tolerance = 1e-3)

    checkEqualsNumeric(getSummaryDF(oplsda.mul.par)[, "Q2(cum)"],
                       0.0531,
                       tolerance = 1e-3)

}

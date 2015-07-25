test_PCA <- function() {

    data(foods)

    fooMN <- as.matrix(foods[, colnames(foods) != "Country"])
    rownames(fooMN) <- foods[, "Country"]
    foo.pca <- opls(fooMN, printL = FALSE, plotL = FALSE)
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
                     printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(corPlsLs[["scoreMN"]][1, 1],
                       2.0513,
                       tolerance = 1e-3)
    cornell.pls <- opls(as.matrix(cornell[, grep("x", colnames(cornell))]),
                     cornell[, "y"],
                     printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(cornell.pls[["vipVn"]][1],
                       1.125,
                       tolerance = 1e-3)

}

test_PLS_multiple <- function() {

    data(lowarp) ## see Eriksson et al. (2001); presence of NAs

    lowarp.pls <- opls(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
                       as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
                                        grepl("^st", colnames(lowarp))]),
                       printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(lowarp.pls[["rotationMN"]][2, 1],
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
                      matrix(sacurine[["sampleMetadata"]][, "genderFc"], ncol = 1),
                      printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.plsda[["summaryDF"]]["h3", "Q2(cum)"],
                       0.5835585,
                       tolerance = 1e-7)
    sac.plsda <- opls(sacurine[["dataMatrix"]],
                      sacurine[["sampleMetadata"]][, "genderFc"],
                      printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.plsda[["vipVn"]][5],
                       0.7034828,
                       tolerance = 1e-7)

    ## splitting the dataset between a reference and a test subsets

    sac.plsdaCrv <- opls(sacurine[["dataMatrix"]],
                         sacurine[["sampleMetadata"]][, "genderFc"],
                         subset = setdiff(1:nrow(sacurine[["dataMatrix"]]), 1:10),
                         printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.plsdaCrv[["summaryDF"]]["h3", "RMSEP"],
                       0.2745173,
                       tolerance = 1e-6)

}

test_sacurine_OPLSDA <- function() {

    data(sacurine)

    sac.oplsda <- opls(sacurine[["dataMatrix"]],
                       matrix(sacurine[["sampleMetadata"]][, "genderFc"],
                              ncol = 1),
                       orthoI = 1,
                       printL = FALSE, plotL = FALSE)
    checkEqualsNumeric(sac.oplsda[["orthoScoreMN"]][1, 1],
                       4.980604,
                       tolerance = 1e-7)

    ## permutation testing

    sac.oplsdaPer <- opls(sacurine[["dataMatrix"]],
                        sacurine[["sampleMetadata"]][, "genderFc"],
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
                         matrix(sacurine[["sampleMetadata"]][, "genderFc"], ncol = 1),
                         scaleC = "pareto",
                         printL = FALSE, plotL = FALSE)
    ## plotVc = "none", verboseC = "none")
    checkEqualsNumeric(sac.plsdaPar[["summaryDF"]]["h2", "Q2(cum)"],
                       0.5210747,
                       tolerance = 1e-7)
    sac.plsdaPar <- opls(sacurine[["dataMatrix"]],
                         sacurine[["sampleMetadata"]][, "genderFc"],
                         scaleC = "pareto",
                         printL = FALSE, plotL = FALSE)
                         ## plotVc = "none", scaleC = "pareto", verboseC = "none")
    checkEqualsNumeric(sac.plsdaPar[["vipVn"]][5],
                       0.7710519,
                       tolerance = 1e-7)

}

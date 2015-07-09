test_PCA <- function() {

    data(foods)

    fooMN <- as.matrix(foods[, colnames(foods) != "Country"])
    rownames(fooMN) <- foods[, "Country"]
    fooPcaLs <- roplsF(fooMN, plotVc = "none", verboseC = "none")
    checkEqualsNumeric(fooPcaLs[["varVn"]][1],
                       6.19,
                       tolerance = 1e-3)
    checkEqualsNumeric(fooPcaLs[["tMN"]][1, 1],
                       1.309494,
                       tolerance = 1e-5)

}

test_PLS_single <- function() {

    data(cornell) ## see Tenenhaus, 1998

    corPlsLs <- roplsF(as.matrix(cornell[, grep("x", colnames(cornell))]),
                       matrix(cornell[, "y"], ncol = 1),
                       plotVc = "none",
                       verboseC = "none")
    checkEqualsNumeric(corPlsLs[["tCompMN"]][1, 1],
                       2.0513,
                       tolerance = 1e-3)
    checkEqualsNumeric(corPlsLs[["vipVn"]][1],
                       1.125,
                       tolerance = 1e-3)

}

test_PLS_multiple <- function() {

    data(lowarp) ## see Eriksson et al. (2001); presence of NAs

    lowPlsLs <- roplsF(as.matrix(lowarp[, c("glas", "crtp", "mica", "amtp")]),
                       as.matrix(lowarp[, grepl("^wrp", colnames(lowarp)) |
                                        grepl("^st", colnames(lowarp))]),
                       plotVc = "none", verboseC = "none")
    checkEqualsNumeric(lowPlsLs[["rMN"]][2, 1],
                       -0.0861,
                       tolerance = 1e-3)
    checkEqualsNumeric(lowPlsLs[["uMN"]]["s5", "h2"],
                       0.596598855,
                       tolerance = 1e-8)

}

test_sacurine_PCA <- function() {

    data(sacurine)

    sacPcaNipLs <- roplsF(sacurine[["profileMN"]],
                          plotVc = "none", verboseC = "none")
    checkEqualsNumeric(sacPcaNipLs[["summaryDF"]]["h8", "R2X(cum)"],
                       0.501,
                       tolerance = 1e-3)
    sacPcaSvdLs <- roplsF(sacurine[["profileMN"]],
                          algoC = "svd",
                          plotVc = "none", verboseC = "none")
    checkEqualsNumeric(sacPcaSvdLs[["tMN"]][1, 1],
                       -8.744009,
                       tolerance = 1e-5)

}

test_sacurine_PLSDA <- function() {

    data(sacurine)

    sacPlsLs <- roplsF(sacurine[["profileMN"]],
                       matrix(sacurine[["sampleDF"]][, "genderVc"], ncol = 1),
                       plotVc = "none", verboseC = "none")
    checkEqualsNumeric(sacPlsLs[["summaryDF"]]["h3", "Q2(cum)"],
                       0.5835585,
                       tolerance = 1e-7)
    checkEqualsNumeric(sacPlsLs[["vipVn"]][5],
                       0.7034828,
                       tolerance = 1e-7)

    ## splitting the dataset between a reference and a test subsets

    sacPlsCrvLs <- roplsF(sacurine[["profileMN"]],
                          matrix(sacurine[["sampleDF"]][, "genderVc"],
                                 ncol = 1),
                          plotVc = "none", testVi = 1:10, verboseC = "none")
    checkEqualsNumeric(sacPlsCrvLs[["summaryDF"]]["h3", "RMSEP"],
                       0.2745173,
                       tolerance = 1e-6)

}

test_sacurine_OPLSDA <- function() {

    data(sacurine)

    sacOplLs <- roplsF(sacurine[["profileMN"]],
                       matrix(sacurine[["sampleDF"]][, "genderVc"],
                              ncol = 1),
                       orthoN = 1,
                       plotVc = "none", verboseC = "none")
    checkEqualsNumeric(sacOplLs[["tOrthoMN"]][1, 1],
                       4.980604,
                       tolerance = 1e-7)

    ## permutation testing

    sacOplPerLs <- roplsF(sacurine[["profileMN"]],
                          matrix(sacurine[["sampleDF"]][, "genderVc"],
                                 ncol = 1),
                          predN = 1,
                          orthoN = 1,
                          permN = 10,
                          plotVc = "none", verboseC = "none")
    checkEqualsNumeric(sacOplPerLs[["permMN"]][1, 1],
                       0.1845708,
                       tolerance = 1e-6)

}


test_sacurine_PLSDA_pareto <- function() {

    data(sacurine)

    sacPlsScaLs <- roplsF(sacurine[["profileMN"]],
                       matrix(sacurine[["sampleDF"]][, "genderVc"], ncol = 1),
                       plotVc = "none", scaleC = "pareto", verboseC = "none")
    checkEqualsNumeric(sacPlsScaLs[["summaryDF"]]["h2", "Q2(cum)"],
                       0.5210747,
                       tolerance = 1e-7)
    checkEqualsNumeric(sacPlsScaLs[["vipVn"]][5],
                       0.7710519,
                       tolerance = 1e-7)

}

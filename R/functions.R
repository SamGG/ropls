#' Printed summary of an R object
#'
#' Displays the class, mode, size and first...last values of the object
#'
#'
#' @param inpMF Input matrix, dataframe or vector
#' @param borderN Number of border (first and last) rows and columns to display
#' @param bigMarkC Big mark separator for summary results
#' @return This function has no output.
#' @author Etienne Thevenot (CEA)
#' @seealso \code{\link{str}}
#' @examples
#'
#' data(sacurine)
#' strF(sacurine[['dataMatrix']])
#' strF(sacurine[['sampleMetadata']])
#'
#' @export strF
strF <- function(inpMF,
                 borderN = 2,
                 bigMarkC = ",") {


    topF <- function() {

        switch(typC,

               vector = {

                   topDF <- data.frame(length = format(length(inpMF), big.mark = bigMarkC),
                                       class = class(inpMF),
                                       mode = mode(inpMF),
                                       typeof = typeof(inpMF),
                                       size = format(object.size(inpMF), units = "Mb"))

                   if(numL)
                       topDF <- cbind.data.frame(topDF,
                                                 data.frame(min = formatC(min(inpMF, na.rm = TRUE),
                                                                digits = 2, format = "g"),
                                                            mean = formatC(mean(inpMF, na.rm = TRUE),
                                                                digits = 2, format = "g"),
                                                            median = formatC(median(inpMF, na.rm = TRUE),
                                                                digits = 2, format = "g"),
                                                            max = formatC(max(inpMF, na.rm = TRUE),
                                                                digits = 2, format = "g")))



               }, ## vector

               matrix = {

                   topDF <- data.frame(dim = paste(format(nrow(inpMF), big.mark = bigMarkC), format(ncol(inpMF), big.mark = bigMarkC), sep = " x "),
                                       class = class(inpMF),
                                       mode = mode(inpMF),
                                       typeof = typeof(inpMF),
                                       size = format(object.size(inpMF), units = "Mb"),
                                       NAs = length(which(is.na(inpMF))))

                   if(numL)
                       topDF <- cbind.data.frame(topDF,
                                                 data.frame(min = formatC(min(inpMF, na.rm = TRUE),
                                                                digits = 2, format = "g"),
                                                            mean = formatC(mean(inpMF, na.rm = TRUE),
                                                                digits = 2, format = "g"),
                                                            median = formatC(median(inpMF, na.rm = TRUE),
                                                                digits = 2, format = "g"),
                                                            max = formatC(max(inpMF, na.rm = TRUE),
                                                                digits = 2, format = "g")))


               }, ## matrix

               data.frame = {

                   claVc <- sapply(inpMF, data.class)

                   if(length(claVc) > 2 * borderN)
                       claVc <- c(head(claVc, borderN),
                                  "...",
                                  tail(claVc, borderN))

                   if(!is.null(names(claVc))) {

                       if(length(claVc) > 2 * borderN)
                           names(claVc)[borderN + 1] <- "..."

                       claDF <- as.data.frame(t(claVc))

                       rownames(claDF) <- ""

                       print(claDF)

                   } else {

                       claVc <- paste(claVc, collapse = " ")
                       class(claVc) <- "table"

                       message(claVc)

                   }

                   topDF <- data.frame(nRow = format(dim(inpMF)[1], big.mark = bigMarkC),
                                       nCol = format(dim(inpMF)[2], big.mark = bigMarkC),
                                       size = format(object.size(inpMF), units = "Mb"),
                                       NAs = length(which(is.na(inpMF))))
               }) ## data.frame

        rownames(topDF) <- ""

        print(topDF)

    } ## topF


    tabF <- function() {

        if(typC %in% c("data.frame", "matrix")) {

            tabDF <- inpMF

            dimAbbVl <- dim(tabDF) > 2 * borderN

            if(is.data.frame(tabDF)) {
                if(dimAbbVl[2]) {
                    bordColVi <- c(1:borderN,
                                   (ncol(tabDF) - borderN + 1):ncol(tabDF))
                } else
                    bordColVi <- 1:ncol(tabDF)

                for(borderI in bordColVi)
                    if(is.factor(tabDF[, borderI]))
                        tabDF[, borderI] <- as.character(tabDF[, borderI])
            }

            if(all(dimAbbVl)) {

                tabDF <- rbind(cbind(tabDF[1:borderN, 1:borderN, drop = FALSE],
                                     ... = rep("...", times = borderN),
                                     tabDF[1:borderN, (ncol(tabDF) - borderN + 1):ncol(tabDF), drop = FALSE]),
                               rep("...", 2 * borderN + 1),
                               cbind(tabDF[(nrow(tabDF) - borderN + 1):nrow(tabDF), 1:borderN, drop = FALSE],
                                     ... = rep("...", times = borderN),
                                     tabDF[(nrow(tabDF) - borderN + 1):nrow(tabDF), (ncol(tabDF) - borderN + 1):ncol(tabDF), drop = FALSE]))

                if(is.matrix(inpMF)) {

                    if(!is.null(rownames(inpMF)) && any(duplicated(rownames(tabDF)))) {
                        rownames(tabDF) <- make.names(rownames(tabDF),
                                                      unique = TRUE)
                    } else if(is.null(rownames(inpMF)))
                        rownames(tabDF) <- c(1:borderN,
                                             "...",
                                             (nrow(inpMF) - borderN + 1):nrow(inpMF))

                    if(is.null(colnames(inpMF)))
                        colnames(tabDF) <- c(1:borderN,
                                             "...",
                                             (ncol(inpMF) - borderN + 1):ncol(inpMF))

                    tabDF <- as.data.frame(tabDF)

                }

                rownames(tabDF)[borderN + 1] <- "..."

            } else if(dimAbbVl[1]) {

                if(is.data.frame(tabDF) && is.null(colnames(tabDF)))
                    colnames(tabDF) <- 1:ncol(tabDF)

                tabDF <- rbind(tabDF[1:borderN, , drop = FALSE],
                               rep("...", ncol(tabDF)),
                               tabDF[(nrow(tabDF) - borderN + 1):nrow(tabDF), , drop = FALSE])

                if(is.matrix(inpMF)) {

                    if(!is.null(rownames(inpMF)) && any(duplicated(rownames(tabDF)))) {
                        rownames(tabDF) <- make.names(rownames(tabDF),
                                                      unique = TRUE)
                    } else if(is.null(rownames(inpMF)))
                        rownames(tabDF) <- c(1:borderN,
                                             "...",
                                             (nrow(inpMF) - borderN + 1):nrow(inpMF))

                    if(is.null(colnames(inpMF)))
                        colnames(tabDF) <- 1:ncol(inpMF)

                    tabDF <- as.data.frame(tabDF)

                }

                rownames(tabDF)[borderN + 1] <- "..."

            } else if(dimAbbVl[2]) {

                tabDF <- cbind(tabDF[, 1:borderN, drop = FALSE],
                               ... = rep("...", times = nrow(tabDF)),
                               tabDF[, (ncol(tabDF) - borderN + 1):ncol(tabDF), drop = FALSE])

                if(is.matrix(inpMF)) {

                    if(!is.null(rownames(inpMF)) && any(duplicated(rownames(tabDF))))
                        rownames(tabDF) <- make.names(rownames(tabDF),
                                                      unique = TRUE)

                    if(is.null(colnames(inpMF)))
                        colnames(tabDF) <- c(1:borderN,
                                             "...",
                                             (ncol(inpMF) - borderN + 1):ncol(inpMF))

                    tabDF <- as.data.frame(tabDF)

                }

            }

        } ## 'data.frame' and 'matrix'

        if(typC == "vector") {

            tabDF <- inpMF

            if(numL)
                tabDF <- round(tabDF, 3)

            if(length(tabDF) > 2 * borderN) {

                tabDF <- c(head(tabDF, borderN),
                           "...",
                           tail(tabDF, borderN))

                if(!is.null(names(inpMF)))
                    names(tabDF)[borderN + 1] <- "..."

            }

            if(!is.null(names(inpMF))) {

                tabDF <- as.data.frame(t(tabDF))

                rownames(tabDF) <- ""

            } else {

                tabDF <- paste(tabDF, collapse = " ")

                message(tabDF)

                return(invisible(NULL))

            }

        } ## vector

        print(tabDF)

    } ## tabF


    if(any(class(inpMF) %in% c("character", "integer", "logical", "numeric", "double"))) {
        typC <- "vector"
     } else
        typC <- class(inpMF)

    numL <- mode(inpMF) %in% c("numeric", "integer", "double")

    if(!(typC %in% c("vector", "matrix", "data.frame"))) {
        str(inpMF)
        return(invisible(NULL))
    }

    topF()

    tabF()

} ## strF


#' fromW4M
#'
#' Creating a ExpressionSet object from the 3 'dataMatrix.tsv',
#' 'sampleMetadata.tsv' and 'variableMetadata.tsv' tabulated files
#'
#' @param dirC Character: directory containing the 3 .tsv files
#' @param namePatternC Character: optional file name pattern common to all three
#' file names (e.g., when you want to distinguish between two sets of files
#' within the same directory)
#' @param fileTableNamesVc Vector of characters: if your file names do not
#' contain the standard 'dataMatrix', 'sampleMetadata', and 'variableMetadata'
#' patterns (e.g. if you use 'profile', 'observation', and 'feature' instead),
#' please indicate them here
#' @param verboseL Logical: should comments be printed?
#' @return ExpressionSet instance
#' @author Etienne Thevenot, \email{etienne.thevenot@@cea.fr}
#' @examples sacSet <- fromW4M(file.path(path.package("ropls"), "extdata"))
#' @rdname fromW4M
#' @export
fromW4M <- function(dirC,
                    namePatternC = "",
                    fileTableNamesVc = c("dataMatrix",
                                         "sampleMetadata",
                                         "variableMetadata"),
                    verboseL = TRUE) {

    tabVc <- c("dataMatrix",
               "sampleMetadata",
               "variableMetadata")

    if(!file.exists(dirC))
        stop("Directory '", dirC, "' was not found.",
             call.=FALSE)

    filVc <- character(length(tabVc))
    names(filVc) <- tabVc

    filAllVc <- list.files(dirC,
                           pattern = "^.*\\.tsv$")

    ## restricting to files with pattern
    if(namePatternC != "")
        filAllVc <- grep(namePatternC, filAllVc, value = TRUE)

    ## restricting to one file for each table
    for(tabC in tabVc) {
        namC <- fileTableNamesVc[tabVc == tabC]
        filC <- grep(namC, filAllVc, value = TRUE)
        if(length(filC) == 0) {
            stop("No file found for the ", tabC, " with ",
                 ifelse(namePatternC != "",
                        paste0("'", namC, "' pattern and "), ""),
                 "a name including '", namC, "' in the '", dirC,
                 "' directory.", call. = FALSE)
        } else if(length(filC) > 1) {
            stop("Several files found for the ", tabC, " with ",
                 ifelse(namePatternC != "", paste0("'", namC, "' pattern and "),
                        ""), "a name including '", namC, "' in the '",
                 dirC, "' directory.", call. = FALSE)
        } else {
            filVc[tabC] <- filC
            ## R standards for row and column names in matrices and data frames
            .checkRformatF(dirC, filC, verboseL)
        }
    }

    ## Loading data

    for(tabC in tabVc) {

        tabDF <- read.table(file.path(dirC, filVc[tabC]),
                            check.names = FALSE,
                            header = TRUE,
                            row.names = 1,
                            sep = "\t",
                            stringsAsFactors = FALSE)
        switch(tabC,
               dataMatrix = {
                   datMN <- as.matrix(tabDF)
               },
               sampleMetadata = {
                   samDF <- tabDF
               },
               variableMetadata = {
                   varDF <- tabDF
               })

    }

    chkL <- .checkW4mFormatF(t(datMN), samDF, varDF)

    if(chkL) {
        TRUE
    } else
        "Problem with the sample or variable names in the tables (see above)"

    eset <- ExpressionSet(assayData = datMN,
                          phenoData = new("AnnotatedDataFrame",
                              data = samDF),
                          featureData = new("AnnotatedDataFrame",
                              data = varDF),
                          experimentData = new("MIAME",
                              title = namePatternC))

    validObject(eset)

    return(eset)

}

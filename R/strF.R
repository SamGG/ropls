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
                                       size = format(object.size(inpMF), big.mark = bigMarkC))

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
                                       size = format(object.size(inpMF), big.mark = bigMarkC),
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

                   if(length(claVc) <= 2) {
                       claVc <- claVc
                   } else
                       claVc <- c(head(claVc, borderN),
                                  "...",
                                  tail(claVc, borderN))

                   if(!is.null(names(claVc))) {

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
                                       size = format(object.size(inpMF), big.mark = bigMarkC),
                                       NAs = length(which(is.na(inpMF))))
               }) ## data.frame

        rownames(topDF) <- ""

        print(topDF)

    } ## topF


    tabF <- function() {

        if(typC %in% c("data.frame", "matrix")) {

            tabDF <- inpMF

            if(class(inpMF) == "matrix")
                tabDF <- as.data.frame(tabDF)

            dimAbbVl <- rep(TRUE, 2)

            ## checking arguments dimensions

            if(all(dim(tabDF) < 4)) {
                stop("Table row or col dimension must be at least 4 to be abbreviated")
            } else
                dimAbbVl <- dimAbbVl & 6 <= dim(tabDF)

            if(all(dimAbbVl)) {

                tabClaVc <- sapply(tabDF, data.class)

                tabDF <- rbind(cbind(tabDF[1:borderN, 1:borderN],
                                     ... = rep("...", times = borderN),
                                     tabDF[1:borderN, (ncol(tabDF) - borderN + 1):ncol(tabDF)]),
                               ifelse(c(tabClaVc[1:borderN], "...", tabClaVc[(ncol(tabDF) - borderN + 1):ncol(tabDF)]) == "factor", NA, "..."),
                               cbind(tabDF[(nrow(tabDF) - borderN + 1):nrow(tabDF), 1:borderN],
                                     ... = rep("...", times = borderN),
                                     tabDF[(nrow(tabDF) - borderN + 1):nrow(tabDF), (ncol(tabDF) - borderN + 1):ncol(tabDF)]))

                rownames(tabDF)[borderN + 1] <- "..."

            } else if(dimAbbVl[1]) {

                tabClaVc <- sapply(tabDF, data.class)

                tabDF <- rbind(tabDF[1:borderN, ],
                               ifelse(tabClaVc == "factor", NA, "..."),
                               tabDF[(nrow(tabDF) - borderN + 1):nrow(tabDF), ])

                rownames(tabDF)[borderN + 1] <- "..."

            } else
                tabDF <- cbind(tabDF[, 1:borderN],
                               ... = rep("...", times = nrow(tabDF)),
                               tabDF[, (ncol(tabDF) - borderN + 1):ncol(tabDF)])

        } ## 'data.frame' and 'matrix'

        if(typC == "vector") {

            tabDF <- inpMF

            if(numL)
                tabDF <- round(tabDF, 3)

            if(length(tabDF) <= 2) {
                tabDF <- tabDF
            } else
                tabDF <- c(head(tabDF, borderN),
                           "...",
                           tail(tabDF, borderN))

            if(!is.null(names(inpMF))) {

                names(tabDF)[borderN + 1] <- "..."

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
        borderN <- min(floor(length(inpMF) / 2) - 1, borderN)
    } else {
        typC <- class(inpMF)
        borderN <- min(floor(nrow(inpMF) / 2) - 1, borderN)
    }

    numL <- mode(inpMF) %in% c("numeric", "integer", "double")


    if(!(typC %in% c("vector", "matrix", "data.frame"))) {
        str(inpMF)
        return(invisible(NULL))
    }

    borderN <- min(floor(nrow(inpMF) / 2) - 1, borderN)


    topF()

    tabF()

} ## strF

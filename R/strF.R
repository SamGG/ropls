strF <- function(inpMF,
                 borN = 2,
                 mrkC = ",") {


    topF <- function() {

        switch(typC,

               vector = {

                   topDF <- data.frame(length = format(length(inpMF), big.mark = mrkC),
                                       class = class(inpMF),
                                       mode = mode(inpMF),
                                       typeof = typeof(inpMF),
                                       size = format(object.size(inpMF), big.mark = mrkC))

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

                   topDF <- data.frame(dim = paste(format(nrow(inpMF), big.mark = mrkC), format(ncol(inpMF), big.mark = mrkC), sep = " x "),
                                       class = class(inpMF),
                                       mode = mode(inpMF),
                                       typeof = typeof(inpMF),
                                       size = format(object.size(inpMF), big.mark = mrkC),
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
                       claVc <- c(head(claVc, borN),
                                  "...",
                                  tail(claVc, borN))

                   if(!is.null(names(claVc))) {

                       names(claVc)[borN + 1] <- "..."

                       claDF <- as.data.frame(t(claVc))

                       rownames(claDF) <- ""

                       print(claDF)

                   } else {

                       claVc <- paste(claVc, collapse = " ")
                       class(claVc) <- "table"

                       message(claVc)

                   }

                   topDF <- data.frame(nRow = format(dim(inpMF)[1], big.mark = mrkC),
                                       nCol = format(dim(inpMF)[2], big.mark = mrkC),
                                       size = format(object.size(inpMF), big.mark = mrkC),
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

            if(all(dim(tabDF) < 4))
                stop("Table row or col dimension must be at least 4 to be abbreviated")
            else
                dimAbbVl <- dimAbbVl & 6 <= dim(tabDF)

            if(all(dimAbbVl)) {

                tabClaVc <- sapply(tabDF, data.class)

                tabDF <- rbind(cbind(tabDF[1:borN, 1:borN],
                                     ... = rep("...", times = borN),
                                     tabDF[1:borN, (ncol(tabDF) - borN + 1):ncol(tabDF)]),
                               ifelse(c(tabClaVc[1:borN], "...", tabClaVc[(ncol(tabDF) - borN + 1):ncol(tabDF)]) == "factor", NA, "..."),
                               cbind(tabDF[(nrow(tabDF) - borN + 1):nrow(tabDF), 1:borN],
                                     ... = rep("...", times = borN),
                                     tabDF[(nrow(tabDF) - borN + 1):nrow(tabDF), (ncol(tabDF) - borN + 1):ncol(tabDF)]))

                rownames(tabDF)[borN + 1] <- "..."

            } else if(dimAbbVl[1]) {

                tabDF <- rbind(tabDF[1:borN, ],
                               rep("...", times = ncol(tabDF)),
                               tabDF[(nrow(tabDF) - borN + 1):nrow(tabDF), ])

                rownames(tabDF)[borN + 1] <- "..."

            } else
                tabDF <- cbind(tabDF[, 1:borN],
                               ... = rep("...", times = nrow(tabDF)),
                               tabDF[, (ncol(tabDF) - borN + 1):ncol(tabDF)])

        } ## 'data.frame' and 'matrix'

        if(typC == "vector") {

            tabDF <- inpMF

            if(numL)
                tabDF <- round(tabDF, 3)

            if(length(tabDF) <= 2) {
                tabDF <- tabDF
            } else
                tabDF <- c(head(tabDF, borN),
                           "...",
                           tail(tabDF, borN))

            if(!is.null(names(inpMF))) {

                names(tabDF)[borN + 1] <- "..."

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
        borN <- min(floor(length(inpMF) / 2) - 1, borN)
    } else {
        typC <- class(inpMF)
        borN <- min(floor(nrow(inpMF) / 2) - 1, borN)
    }

    numL <- mode(inpMF) %in% c("numeric", "integer", "double")


    if(!(typC %in% c("vector", "matrix", "data.frame"))) {
        str(inpMF)
        return(invisible(NULL))
    }

    borN <- min(floor(nrow(inpMF) / 2) - 1, borN)


    topF()

    tabF()


} ## strF

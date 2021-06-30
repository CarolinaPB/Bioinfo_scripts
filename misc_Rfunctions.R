
# Get details out of attribute field GFF
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}


# Read gff file into data.table
gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = data.table(read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                              header=FALSE, comment.char="#", nrows = nrows,
                              colClasses=c("character", "character", "character", "integer",
                                           "integer",
                                           "character", "character", "character", "character")))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  return(gff)
}

# get list with values from boxplot like whisker, quantiles, etc
ggplot2_boxplot_stats <- function(x){

    quartiles <- as.numeric(quantile(x,
                                     probs = c(0.25, 0.5, 0.75)))

    names(quartiles) <- c("25th percentile",
                          "50th percentile\n(median)",
                          "75th percentile")

    IQR <- diff(quartiles[c(1,3)])

    upper_whisker <- max(x[x < (quartiles[3] + 1.5 * IQR)])
    lower_whisker <- min(x[x > (quartiles[1] - 1.5 * IQR)])

    upper_dots <- x[x > (quartiles[3] + 1.5*IQR)]
    lower_dots <- x[x < (quartiles[1] - 1.5*IQR)]

    return(list("quartiles" = quartiles,
                "25th percentile" = as.numeric(quartiles[1]),
                "50th percentile\n(median)" = as.numeric(quartiles[2]),
                "75th percentile" = as.numeric(quartiles[3]),
                "IQR" = IQR,
                "upper_whisker" = upper_whisker,
                "lower_whisker" = lower_whisker,
                "upper_dots" = upper_dots,
                "lower_dots" = lower_dots))
  }



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
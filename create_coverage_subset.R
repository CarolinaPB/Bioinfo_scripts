library(data.table)
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="coverage file", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="coverage", 
              help="prefix to name the output file", metavar="character")
); 


opt_parser = OptionParser(description="\nTo use with the output coverage table from mafCoverage\nWrites a using the chromosomes of both genomes as reference and as query. You need to pay attention to which one you want as which",option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please add a file -f <file.txt>", call.=FALSE)
}

tab <- fread(opt$file, fill=T)

# tab <- fread("coverage.txt")
tab <- tab[,1:4]
tab_subset <- tab[!startsWith(`# referenceSpecies/Chr`, "ChrUn_random_")][!startsWith(`querySpecies/Chr`, "ChrUn_random_")][order(`# referenceSpecies/Chr`, -coverage)][coverage>=0.01]

top2_tab <- tab_subset[, head(.SD, 2), by=`# referenceSpecies/Chr`]

filename <- paste0(opt$prefix, "_coverage_subset.txt")
filename_sorted <- paste0(opt$prefix, "_coverage_subset_sorted.txt")

fwrite(top2_tab, filename, sep = "\t")


command <- paste0("sort -k1,1 -k4,4r -V ", filename, " > ", filename_sorted)
system(command)
system(paste0("rm ", filename))

#!/usr/bin/env Rscript

#packagesList <- c("data.table", "optparse")
#newPackages <- packagesList[!(packagesList %in% installed.packages()[,"Package"])]
#if(length(newPackages) > 0) {
#  install.packages(newPackages)
#}
library(data.table)
library("optparse")

option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="structural var vcf file", metavar="character"),
  make_option(c("-i", "--intermediate"), type="character", default="vcf_subset.txt", 
              help="intermediate output file name [default= %default]", metavar="character"),
  make_option(c("-c", "--svcaller"), type="character", default=NULL, 
              help="SV caller can be svim or svim-asm", metavar="character")
); 


opt_parser = OptionParser(description="\nThis script is optimized to use with SVIM vcf output\nThis script will output a summary of the structural variants found in a given vcf to stdout.\nYou can redirect the output to a new file.\nThis script creates an intermediate file with three columns: chromosome name, name of structural variant and size of SV",option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$vcf)){
  print_help(opt_parser)
  stop("Please add a file -f <file.txt>", call.=FALSE)
}

command <- paste0("sh /lustre/nobackup/WUR/ABGC/moiti001/TOOLS/scripts/get_struct_var_vcf_subset.sh -vcf=", opt$vcf, " -o=", opt$intermediate," -c=", opt$svcaller)
system(command)


struct <- fread(opt$intermediate)
struct <- fread("vcf_subset.txt")
colnames(struct) <- c("chromosome", "var", "var_size")

# Total number of each variant

print("TOTAL NUMBER OF STRUCTURAL VARIANTS")
struct[,.N, by=var]

# Min and max size of variant

print("MIN AND MAX SIZE OF EACH VARIANT")
struct_maxsize <- struct[,max(abs(as.numeric(var_size))),by=var]
struct_min <- struct[,min(abs(as.numeric(var_size))), by=var]
struct_size_limits <- merge(struct_min, struct_maxsize, by="var")
colnames(struct_size_limits) <- c("var", "min_size", "max_size")
struct_size_limits

# Number of each variant per chr

print("NUMBER OF EACH VARIANT PER CHR")
struct[, .N, by=c("chromosome","var")]

# Number of variants found by chr

print("NUMBER OF VARIANTS FOUND IN EACH CHR")
struct[, .N, by=chromosome]

# Order by number of variants
print("NUMBER OF VARIANTS FOUND IN EACH CHR ORDERED BY NUM OF VAR")
struct[, .N, by=chromosome][order(-N)]

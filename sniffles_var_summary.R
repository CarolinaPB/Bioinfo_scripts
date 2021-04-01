library(data.table)
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="sniffles VCF", metavar="character")
  # make_option(c("-p", "--prefix"), type="character", default="coverage", 
  #             help="prefix to name the output file", metavar="character")
); 


opt_parser = OptionParser(description="Get SV summary table out of sniffles VCF",option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please add a file -f <file.txt>", call.=FALSE)
}

command <- paste0("sh /lustre/nobackup/WUR/ABGC/moiti001/TOOLS/scripts/get_sniffles_summary.sh -vcf=", opt$file)
system(command)

# system(paste0("vcf-query -f '%CHROM\t%INFO/SVTYPE\t%POS\t%INFO/END\t%INFO/CHR2\t%INFO/SVLEN\t%INFO/RE\n'" , opt$file, " >> sniffles_summary.txt"))

sniffles_summary <- fread("sniffles_summary.txt")




# Total number of each variant

print("TOTAL NUMBER OF STRUCTURAL VARIANTS")
sniffles_summary[,.N, by=SV]

# Min and max size of variant

print("MIN AND MAX SIZE OF EACH VARIANT")
struct_maxsize <- sniffles_summary[,max(abs(as.numeric(Len))),by=SV]
struct_min <- sniffles_summary[,min(abs(as.numeric(Len))), by=SV]
struct_size_limits <- merge(struct_min, struct_maxsize, by="SV")
colnames(struct_size_limits) <- c("var", "min_size", "max_size")
struct_size_limits

# Number of each variant per chr

print("NUMBER OF EACH VARIANT PER CHR")
sniffles_summary[, .N, by=c("Chr","SV")]

# Number of variants found by chr

print("NUMBER OF VARIANTS FOUND IN EACH CHR")
sniffles_summary[, .N, by=Chr]

# Order by number of variants
print("NUMBER OF VARIANTS FOUND IN EACH CHR ORDERED BY NUM OF VAR")
sniffles_summary[, .N, by=Chr][order(-N)]

# Number of each GT per variant
print("GENOTYPE PER SV")
sniffles_summary[,.N, by=c("SV", "GT")][order(SV)]

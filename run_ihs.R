library(vcfR)
library(rehh)
library(data.table)
library(R.utils)
library(optparse)

option_list = list(
  make_option(c("-r", "--rds"), type="character", default=NULL, 
              help="first rds to be used", metavar="character"),
  make_option(c("-p", "--pop"), type="character", default=NULL, 
              help="name of population 1", metavar="character"),
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


print(paste0("Running iHS on ", opt$pop))


wgscan <- readRDS(file = opt$rds)


#iHS ANALYISIS
###########################################################################
# Chromosome list for your species
# Excluding chr29-not available for galGal6 ref genome
chrs <- c(seq(1,28), seq(30,33), "W", "Z")
hap_file = opt$vcf

# calculate genome-wide iHS values
print("calculating ihs")
palette("default")
wgscan.ihs <- ihh2ihs(wgscan)

print("saving ihs")
saveRDS(wgscan.ihs, file= paste0("iHS_", opt$pop, ".rds"))

pdf(paste0("iHS_", opt$pop,".pdf"))
manhattanplot(wgscan.ihs,
              pval = TRUE,
              threshold = 4,
              main = paste0("iHS population (", opt$pop, ")"))
dev.off()

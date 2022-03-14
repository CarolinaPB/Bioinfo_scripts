library(vcfR)
library(rehh)
library(data.table)
library(R.utils)
library(optparse)


option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="vcf to be used", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output prefix", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(paste0("Running iHS on ", opt$vcf))
print(paste0("output prefix: ", opt$output))

#iHS ANALYISIS
###########################################################################
# Chromosome list for your species
# Excluding chr29-not available for galGal6 ref genome
chrs <- c(seq(1,28), seq(30,33), "W", "Z")
hap_file = opt$vcf

for(i in chrs) {

    hh <- data2haplohh(hap_file = hap_file,
                       chr.name = i,
                       polarize_vcf = FALSE,
                       vcf_reader = "data.table")
    scan <- scan_hh(hh)
    if (i == "1") {
      wgscan <- scan
    } else {
      wgscan <- rbind(wgscan, scan)
  }

}

print("saving wgscan")
saveRDS(wgscan, file=paste0("wgscan_", opt$output, ".rds"))
write.csv(wgscan,paste0("wgscan_", opt$output, ".csv"))

# calculate genome-wide iHS values
print("calculating ihs")
palette("default")
wgscan.ihs <- ihh2ihs(wgscan)

print("saving ihs")
saveRDS(wgscan.ihs, file= paste0("iHS_", opt$output, ".rds"))

pdf(paste0("iHS_", opt$output,".pdf"))
manhattanplot(wgscan.ihs,
              pval = TRUE,
              threshold = 4,
              main = paste0("iHS Ethiopian chicken populations (", opt$output, ")"))
dev.off()

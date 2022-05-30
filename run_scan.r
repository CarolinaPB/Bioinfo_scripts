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

# TODO: add the chromosomes for your species 
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

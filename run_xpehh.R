library(vcfR)
library(rehh)
library(data.table)
library(R.utils)
library(optparse)


option_list = list(
  make_option(c("-f", "--rds1"), type="character", default=NULL, 
              help="first rds to be used", metavar="character"),
  make_option(c("-s", "--rds2"), type="character", default=NULL, 
           help="second rds to be used", metavar="character"),
  make_option(c("-p", "--pop1"), type="character", default=NULL, 
              help="name of population 1", metavar="character"),
  make_option(c("-g", "--pop2"), type="character", default=NULL, 
              help="name of population 2", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$pop1)

print(paste0("Running xpehh ", opt$pop1, " vs ", opt$pop2))
print(paste0("Files: ", opt$rds1, ", ", opt$rds2))

wgscan_1 <- readRDS(file = opt$rds1)
wgscan_2 <- readRDS(file = opt$rds2)


# Run xpehh
agroeco.xpehh<-ies2xpehh(wgscan_1,wgscan_2, opt$p,opt$g)
saveRDS(agroeco.xpehh, file= paste0("xpehh_", opt$output, ".rds"))

jpeg(paste0("xpehh_", opt$p, "_", opt$g,".jpeg"))
manhattanplot(agroeco.xpehh,
              main = paste0("XPEHH populations (", opt$p, " - ", opt$g, ")"))
dev.off()

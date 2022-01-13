# Plot association plots (GWAS)
options(warn=-1)
suppressMessages(library(qqman))
args = commandArgs(trailingOnly=TRUE)

# arg[1] - name of file to save
# arg[2] - name of table to use where the following columns must be present and have these column names
# CHR
# BP 
# P
# SNP
# arg[3] - title of plot

pdf(file = args[1])



gwas=read.table(args[2], header=TRUE)
manhattan(gwas, main = args[3], col = c("black", "deeppink2"), cex = 0.5)



dev.off()
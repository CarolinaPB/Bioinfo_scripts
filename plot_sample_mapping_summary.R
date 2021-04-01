#!/usr/bin/env Rscript

library(data.table)

# Plot histogram of mapping stats: mean coverage, mapping rate and mean mapping quality. 
# Meant to work with the output of create_qualimap_summary.sh 
# input file should have 4 columns in this order:
#   sample
#   mean coverage
#   mapping rate
#   mean mapping quality
# 
# Usage: Rscript plot_sample_mapping_summary.R <qualimap_summary.txt>
# Results will be saved to mapping_stats directoryy

args = commandArgs(trailingOnly=TRUE)

ifelse(!dir.exists(file.path("mapping_stats")), dir.create(file.path("mapping_stats")), "Saving plots to mapping_stats")

mapping_stats <- fread(args[1])
# mapping_stats <- fread("sample_quality_summary.tsv")
colnames(mapping_stats) <- c("sample", "mean_coverage", "mapping_rate", "mean_mapping_qual")


mapping_stats[,("mean_coverage"):=lapply(.SD, gsub, pattern="X", replacement=""), .SDcols="mean_coverage"]
mapping_stats[,("mapping_rate"):=lapply(.SD, gsub, pattern="%", replacement=""), .SDcols="mapping_rate"]

png("mapping_stats/mean_coverage.png",width=600, height=350)
hist(as.numeric(mapping_stats$mean_coverage), breaks=30, main="Mean coverage", xlab="mean coverage (X)", ylab = "# samples")
dev.off()

png("mapping_stats/mapping_rate.png",width=600, height=350)
hist(as.numeric(mapping_stats$mapping_rate), breaks=30, main="Mapping rate", xlab="mapping rate (%)", ylab = "# samples")
dev.off()

png("mapping_stats/mapping_quality.png",width=600, height=350)
hist(as.numeric(mapping_stats$mean_mapping_qual), breaks=30, main="Mean mapping quality", xlab="Mean mapping quality", ylab = "# samples")
dev.off()

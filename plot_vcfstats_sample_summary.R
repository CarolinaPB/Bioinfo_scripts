#!/usr/bin/env Rscript

library(data.table)

# Plot histogram of num ref homozygous, num non ref hom, num heterozygous, ts/tv, indels from the bcftools stats output
# Usage:  Rscript plot_vcfstats_sample_summary.r <stats>
# Will save plots in a directory called vcf_plots

args = commandArgs(trailingOnly=TRUE)


# sample_stats <- fread("concat/per_sample_stats.txt", drop = c("# PSC", "[2]id"), )
# create vcf_plots dir if it doesn't exist. Plots will be saved here
ifelse(!dir.exists(file.path("vcf_plots")), dir.create(file.path("vcf_plots")), "Saving plots to vcf_plots")


command <- paste0('zgrep PSC ',args[1], " > vcf_plots/per_sample_stats.txt")
system(command)

sample_stats <- fread("vcf_plots/per_sample_stats.txt", drop = c("# PSC", "[2]id"))


sample_stats[,"ts_tv":= `[7]nTransitions`/`[8]nTransversions`]



png("vcf_plots/num_ref_hom.png",width=600, height=350)
hist(sample_stats$`[4]nRefHom`, main="Num Ref Homozygous", xlab="# ref homozygous", ylab = "# samples", breaks=20)
dev.off()

png("vcf_plots/ts_tv.png",width=600, height=350)
hist(sample_stats$ts_tv, main="Ts/Tv", xlab="ts/tv", ylab = "# samples", breaks = 50)
dev.off()

png("vcf_plots/num_nonref_hom.png",width=600, height=350)
hist(sample_stats$`[5]nNonRefHom`, main="Num non Ref Homozygous", xlab="# non-ref homozygoys", ylab = "# samples", breaks=20)
dev.off()

png("vcf_plots/num_het.png",width=600, height=350)
hist(sample_stats$`[6]nHets`, main="Num Heterozygous", xlab = "# heterozygous", ylab = "# samples", breaks=20)
dev.off()

png("vcf_plots/num_indels.png",width=600, height=350)
hist(sample_stats$`[9]nIndels`, main="Num Indels", xlab = "# indels", ylab = "# samples", breaks=20)
dev.off()

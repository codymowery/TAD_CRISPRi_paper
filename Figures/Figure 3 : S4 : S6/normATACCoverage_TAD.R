library(tidyverse)
library(vroom)
library(DESeq2)
library(GenomicRanges)
library(rtracklayer)
library(edgeR)

## Files and code from Freimer et al 2022
## Subset to TAD coordinates via "for file in ../*; do grep chr2 ${file} | bedtools intersect -a stdin -b beds/tad_subsets/ref.bed > ${file}.tad; done"

# read in counts
counts <- vroom('count_mat_peaks_cluster150bp_peak_size_350bp.txt') ## From Freimer et al 2022
counts <- as.data.frame(counts)
rownames(counts) <- counts$name

# remove well position
colnames(counts) <- gsub("_[A-Z][0-9]{0,2}$", "", colnames(counts))

# Select relevant samples and generate count matrix
count_mat <- as.matrix(dplyr::select(counts, matches('Donor'), -matches('AAVS1_4')))

# Filter low count reads
min_cpm <- ceiling(10 / (min(colSums(count_mat)) / 1e6))
count_mat_filtered <- count_mat[rowSums(edgeR::cpm(count_mat) > min_cpm) >= 3, ]
print(dim(count_mat_filtered))

# Get DESeq2 size factors based on reads in peaks to normalize coverage
norm <- estimateSizeFactorsForMatrix(count_mat_filtered)
norm
colSums(count_mat_filtered)

# # Get control sample
control_samples <- grep("AAVS1", names(norm), value = T)

# Calculate coverage at TAD locus
locus_files <- dir('~/Dropbox (Gladstone)/Papers/CTLA4_TAD/code/JakeATAC/tadSubsets/', pattern = 'insertions', full.names = T)

# compute coverage using TAD_stop as width
calc_coverage <- function(x) {
    insertions <- import(grep(x, locus_files, value = T))
    cov <- IRanges::coverage(insertions + 50, width = 203967032)/norm[x] ## TAD_stop==203967032
}

# Control coverage at locus
control_coverage <- lapply(control_samples, calc_coverage)
control_coverage_average <- Reduce("+", control_coverage)/length(control_coverage)
export(GRanges(control_coverage_average), 'CTRL_coverage_at_tad.bw')


# samples <- grep('AAVS1', unique(gsub('Donor_[1-4]_', '', names(norm))), value = T, invert = T)
KO_coverage <- function(KO) {
  KO_samples <- grep(KO, names(norm), value = T)
  print(KO_samples)
  
  KO_coverage <- lapply(KO_samples, calc_coverage)
  KO_coverage_average <- Reduce("+", KO_coverage)/length(KO_coverage)
  export(GRanges(KO_coverage_average), paste0(KO, '_coverage_at_tad.bw'))
}
samples = c('ZNF217','IRF4','JAK3')
lapply(samples, KO_coverage)
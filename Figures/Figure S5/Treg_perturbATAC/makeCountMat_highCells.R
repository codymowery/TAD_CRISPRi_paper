library('tidyverse')
library('rtracklayer')
library('GenomicRanges')
library('GenomicAlignments')
list.files('output/merged_peaks/beds/')
final_peaks_file="output/merged_peaks/beds/peaks_cluster_150bp_peak_size_350bp.bed"
output_file="output/counts/count_mat_peaks_cluster150bp_peak_size_350bp.txt"
final_peaks <- import.bed(final_peaks_file)
names(final_peaks) <- paste0('peak', seq(1:length(final_peaks)))
final_peaks$name <- names(final_peaks)
get_count_matrix <- function(insertion_file) {
  insertions_df <- vroom::vroom(insertion_file, delim = '\t',
                                col_names = c('chr', 'start', 'stop', 'readid', 'score', 'strand'),
                                col_select = c(chr, start))
  insertions <- GRanges(seqnames = insertions_df$chr, IRanges(start = insertions_df$start+1, width = 1))
  rm(insertions_df)
  
  current_sample <- str_extract(insertion_file, 'Donor[A-Z]_[[:graph:]]')
  
  so <- summarizeOverlaps(final_peaks, insertions)
  count_mat <- assays(so)$counts
  colnames(count_mat) <- current_sample
  return(count_mat)
}

insertion_files <- list.files("output/beds/", pattern = 'highCells.insertions.bed$', full.names = T)
peak_counts <- lapply(insertion_files, get_count_matrix)


# Combine each count into 1 matrix and annotate peak information
combined_count_mat <- Reduce(cbind, peak_counts) %>% as_tibble(rownames = 'peakName')
annotated_count_mat <- final_peaks %>% as_tibble() %>% select(-score) %>% inner_join(., combined_count_mat, by = c('name' = 'peakName'))
annotated_count_mat
write_tsv(annotated_count_mat, output_file)

output_file
list.files()
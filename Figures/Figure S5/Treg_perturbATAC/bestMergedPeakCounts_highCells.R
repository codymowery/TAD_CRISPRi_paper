library('tidyverse')
library('GenomicRanges')
list.files('output/merged_peaks/beds/')
all_peaks <- read_rds('output/merged_peaks/beds/all_peaks_highCells.RDS')
list.files('output/beds/')
insertion_files <- list.files('output/beds/', pattern = 'highCells.insertions.bed$', full.names = T)

insertion_files
read_overlap_count <- function(reads, peaks) {
  read_count <- countOverlaps(reads, peaks)
  read_count_summary <- tibble(peak_size = unique(peaks$peak_size),
                               aggregate_distance = unique(peaks$aggregate_cluster_distance),
                               count_any_peak = mean(read_count > 0),
                               count_1_peak = mean(read_count == 1),
                               count_multi_peaks = mean(read_count > 1))
}
dir.create('output/merged_peaks/best_combo_count/')
for(insertion_file in insertion_files){
  # -----------------------------------------------------
  # Summarize best aggregation and peak size combination for all samples
  current_sample <- gsub('.insertions.bed', '', str_extract(insertion_file, 'Donor[A-Z]_[[:graph:]]'))
  
  insertions_df <- vroom::vroom(insertion_file, delim = '\t',
                                col_names = c('chr', 'start', 'stop', 'readid', 'score', 'strand'),
                                col_select = c(chr, start))
  insertions <- GRanges(seqnames = insertions_df$chr, IRanges(start = insertions_df$start+1, width = 1))
  rm(insertions_df)
  
  sample_count_summary <- lapply(all_peaks, function(x) read_overlap_count(reads = insertions, peaks = x)) %>%
    bind_rows() %>%
    mutate(sample = current_sample)
  
  write_tsv(x = sample_count_summary, file = paste0('output/merged_peaks/best_combo_count/', current_sample, '_peaks_count.tsv'))
}
list.files()
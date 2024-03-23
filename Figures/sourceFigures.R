## Everything in hg38 coordinates unless otherwise noted

library(tidyverse)
library(cowplot)
library(GenomicRanges)
library(ggh4x)
library(fuzzyjoin)
library(rtracklayer)
library(ComplexHeatmap)
library(svglite)
library(ggrepel)
library(ggpubr)
theme_set(theme_cowplot())

# Gene and sgRNA information
TAD_start=202527032
TAD_stop=203967032
CD28_start=203706475
CD28_stop=203738912
CTLA4_start=203867771
CTLA4_stop=203873965
ICOS_start=203936763
ICOS_stop=203961577
sgRNA_ref = read_csv('Library.tsv',col_names = F) %>% 
  dplyr::select(sgrna=X1,sequence=X2) %>% 
  tidyr::separate(sgrna, into = c('library','number','start')) %>% 
  dplyr::mutate(start=as.double(start))

## CRISPRi colors
orange = '#c18b33'
high_col = colorspace::darken(orange,0)
midHigh_col = colorspace::lighten(orange,0.75)
blue = '#2a5080'
low_col = colorspace::darken(blue,0)
midLow_col = colorspace::lighten(blue,0.75)
mid_col = 'grey90'

## Colors for trans RNA and ATAC
trans_lowCol='#2786a0'
trans_midCol='white'
trans_highCol='#e1be6a'
rna_lowCol='#001663'
rna_midCol='white'
rna_highCol='#7e4030'

## CRISPRi Datasets ---------
merged = read_rds('CRISPRi_merged.rds')

## References ----------
tad_genes=read_rds('tadGenes_ggplotRef.rds')
full_geneRef = function(genes=tad_genes,range=c(TAD_start,TAD_stop)){
  ggplot(genes, aes(start+(width/2),gene_num,width=width,label=sym_strand))+
    geom_tile(fill='black')+
    coord_cartesian(xlim = range, ylim=c(0,max(tad_genes$gene_num)+3), expand=FALSE)+
    geom_text(aes(x = start+(width/2), y = gene_num+1.5), size=rel(2))+
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.line.x = element_blank(), axis.line.y = element_line(color='white'), axis.ticks = element_blank())}

## Gene models
# Downloaded from UCSC Table Browser using Known Gene and coordinates 
# Run with pyGenomeTracks in /Users/codymowery/Dropbox (Gladstone)/Papers/CTLA4_TAD/code/ctla4_geneBody/
infile = import.bed('hgTables.bed')
sub = infile[gsub("\\..*","",infile$name) %in% read_rds('1tpm_donor1001_TregTeff_U_S_txList.rds')]
export.bed(sub,'hgTables_filtered1TPM.bed')
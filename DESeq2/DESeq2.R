library(tidyverse)
library(DESeq2)
sessionInfo()

## Script takes in sgRNA counts from Mageck and uses DESeq2 to identify sgRNAs enriched in high or low expressing cell populations

########################
## 210818-1 Treg Screens
########################
list.files('mageck/210818-1_mageckOuts/')
deseq2_treg = function(target,construct){
  d1_file = paste0('mageck/210818-1_mageckOuts/D1_',target,'.count.txt')
  d2_file = paste0('mageck/210818-1_mageckOuts/D2_',target,'.count.txt')
  cts = read_tsv(d1_file) %>% 
    full_join(., read_tsv(d2_file),by=c('sgRNA','Gene')) %>% 
    as.data.frame(.)
  rownames(cts)=paste0(cts$sgRNA)
  cts = subset(cts, select = -c(Gene,sgRNA))
  coldata = data.frame(row.names = colnames(cts), donor = c('D1','D1','D2','D2'), condition = c('LO','HI','LO','HI'), sample = colnames(cts))
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ donor + condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  ## Differential expression testing
  dds <- DESeq(dds)
  res <- results(dds)
  return(as_tibble(res,rownames = 'sgrna') %>% dplyr::mutate(target=target))
}
targets=c('CD28','CTLA4_0h','CTLA4_6h','ICOS')
merged_treg = tibble()
for(target in targets){
  merged_treg = bind_rows(
    merged_treg,
    deseq2_treg(target))
}
merged_treg_clean = merged_treg %>% tidyr::separate(sgrna, into=c('library','name','start'),'_') %>% dplyr::mutate(start=as.double(start))

range=c(CD28_start-20e3,ICOS_stop+20e3)
ggplot(merged_treg_clean,
       aes(start,-log2FoldChange,
           color=I(ifelse(padj<0.05,'black','grey50')),
           alpha=I(ifelse(padj<0.05,1,0.1))))+
  annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  geom_hline(yintercept = 0, color='gray88')+
  geom_point(data=merged_treg_clean %>% dplyr::filter(padj>0.05), shape=16, stroke=0)+
  geom_point(data=merged_treg_clean %>% dplyr::filter(padj<0.05), shape=16, stroke=0)+
  theme_classic()+
  scale_color_manual(values=c('black','grey50'))+
  theme(legend.position = "none")+
  scale_x_continuous(limits=range)+
  facet_grid(rows=vars(target))

table(merged_treg_clean$target)
merged_treg_clean %>% 
  dplyr::mutate(celltype='Treg',
                condition=ifelse((target=='CD28')|(target=='CTLA4_0h'),'0h',
                                 ifelse(target=='CTLA4_6h','6h','24h')),
                target2=ifelse((target=='CD28')|(target=='ICOS'),target,'CTLA4')) %>% 
  dplyr::select(-target) %>% dplyr::rename(target=target2) %>% 
  write_rds('DESeq2/210818-1_deseq2.rds')

########################
## 211028-1 Tconv Screens
########################
list.files('mageck/211028-1_mageckOuts/')
deseq2_tconv = function(target,construct){
  d1_file = paste0('mageck/211028-1_mageckOuts/D1_',target,'.count.txt')
  d2_file = paste0('mageck/211028-1_mageckOuts/D2_',target,'.count.txt')
  cts = read_tsv(d1_file) %>% 
    full_join(., read_tsv(d2_file),by=c('sgRNA','Gene')) %>% 
    as.data.frame(.)
  rownames(cts)=paste0(cts$sgRNA)
  cts = subset(cts, select = -c(Gene,sgRNA))
  coldata = data.frame(row.names = colnames(cts), donor = c('D1','D1','D2','D2'), condition = c('LO','HI','LO','HI'), sample = colnames(cts))
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ donor + condition)
  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  ## Differential expression testing
  dds <- DESeq(dds)
  res <- results(dds)
  return(as_tibble(res,rownames = 'sgrna') %>% dplyr::mutate(target=target))
}

targets=c('CD28','CTLA4','ICOS')
merged_tconv = tibble()
for(target in targets){
  merged_tconv = bind_rows(
    merged_tconv,
    deseq2_tconv(target))
}
merged_tconv_clean = merged_tconv %>% tidyr::separate(sgrna, into=c('library','name','start'),'_') %>% dplyr::mutate(start=as.double(start))

range=c(CD28_start-20e3,ICOS_stop+20e3)
ggplot(merged_tconv_clean,
       aes(start,-log2FoldChange,
           color=I(ifelse(padj<0.05,'black','grey50')),
           alpha=I(ifelse(padj<0.05,1,0.1))))+
  annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  geom_hline(yintercept = 0, color='gray88')+
  geom_point(data=merged_tconv_clean %>% dplyr::filter(padj>0.05), shape=16, stroke=0)+
  geom_point(data=merged_tconv_clean %>% dplyr::filter(padj<0.05), shape=16, stroke=0)+
  theme_classic()+
  scale_color_manual(values=c('black','grey50'))+
  theme(legend.position = "none")+
  scale_x_continuous(limits=range)+
  facet_grid(rows=vars(target))

merged_tconv_clean %>% 
  dplyr::mutate(celltype='Tconv',
                condition=ifelse(target=='CD28','0h',
                                 ifelse(target=='CTLA4','6h','24h'))) %>% 
  write_rds('DESeq2/211028-1_deseq2.rds')

########################
## Merge Screens
########################
bind_rows(read_rds('DESeq2/210818-1_deseq2.rds'),
          read_rds('DESeq2/211028-1_deseq2.rds')) %>% 
  write_rds('DESeq2/merged_deseq2.rds')

## merged_deseq2.rds provided as Supplementary Table 3

## Everything in hg38 coordinates 
source('sourceFigures.R')

## CRISPRi compare
library(tidyverse)
library(ComplexHeatmap)
library(BSgenome.Hsapiens.UCSC.hg38)
library(fuzzyjoin)
library(cowplot)
library(ggrepel)
library(ComplexHeatmap)
library(svglite)
theme_set(theme_cowplot())

# Look at concordance of ZIM3 and KRAB from limited 1-donor test
merged_compare=read_rds('210514-1_techrep_mergedResults.rds')
padj_thresh = 0.05
plot = merged_compare %>% 
  dplyr::select(-c(baseMean,lfcSE,pvalue)) %>% 
  pivot_wider(names_from = construct, values_from = c(log2FoldChange,padj)) %>% 
  dplyr::mutate(Significance = as.factor(ifelse((padj_KRAB<padj_thresh)&(padj_ZIM3<padj_thresh),'Both',
                                                ifelse((padj_KRAB<padj_thresh)&(padj_ZIM3>padj_thresh),'KRAB only',
                                                       ifelse((padj_KRAB>padj_thresh)&(padj_ZIM3<padj_thresh),'ZIM3 only','Neither')))))
ggplot(plot,aes(-log2FoldChange_ZIM3,-log2FoldChange_KRAB,color=Significance))+
  geom_point(data = . %>% dplyr::filter(Significance=='Neither'))+
  geom_point(data = . %>% dplyr::filter(Significance!='Neither'))+
  coord_fixed()+
  scale_color_manual(values=c('#792ab6','#edc948','#9e9e9e','#d95f03'))+
  theme(legend.position = 'bottom', legend.justification = 'center')

# Plot ZIM3 vs KRAB tracks
region=c(202500000,203980000)
merged_compare %>%
  dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
  dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
  dplyr::mutate(col=ifelse(log2FoldChange>0&padj<padj_thresh,low_col,
                           ifelse(log2FoldChange<0&padj<padj_thresh,high_col,mid_col))) %>% 
  dplyr::arrange(desc(pvalue)) %>% 
  ggplot()+
  annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
  geom_hline(yintercept = 0, linetype='solid', color='black')+
  geom_point(aes(start,-log10(pvalue),color=I(col)),shape=16,size=rel(3),alpha=0.6)+
  facet_grid(rows=vars(target,construct),switch = 'y',space='fixed', scales = 'free_y')+
  coord_cartesian(xlim = region, expand=FALSE, clip = "off")+
  theme(strip.background = element_blank(),
        panel.spacing = unit(0.7, "lines"))+
  NULL


## Donor:Donor Correlation
sgRNA_seqs = read_csv('UniqueLibraryRef.csv', col_names = F) %>% 
  tidyr::separate(X1,into=c('library','name','start')) %>% dplyr::select(name=X3,sequence=X2,start)
fnames = c(list.files('210818-1_mageckOuts/',pattern = 'normalized',full.names = T),
           list.files('211028-1_mageckOuts/',pattern = 'normalized',full.names = T))
norm = tibble()
for(file in fnames){
  norm=bind_rows(norm,
                 read_tsv(file) %>% 
                   dplyr::mutate(fname=str_replace_all(file,'Data/|mageck_outs\\/\\/|\\.count_normalized\\.txt|_\\.count_normalized\\.txt','')) %>% 
                   dplyr::rename(LO=3,HI=4))
} 
plot = norm %>% 
  dplyr::filter(!grepl('0h_CTLA4',fname)) %>% 
  dplyr::mutate(celltype=ifelse(grepl('210818-1',fname),'Treg','Tconv'),
                donor=ifelse(grepl('D1',fname),'D1','D2'),
                target=ifelse(grepl('CD28',fname),'CD28',
                              ifelse(grepl('ICOS',fname),'ICOS','CTLA4')),
                lfc=log2(HI/LO)) %>% 
  dplyr::filter(!grepl('CTLA4_0h',fname)) %>%  
  tidyr::separate(sgRNA, into=c('library','name','start'),'_') %>% dplyr::mutate(start=as.double(start)) %>%
  dplyr::select(library,name,start,celltype,donor,target,lfc) %>% 
  dplyr::filter(abs(lfc)<Inf) %>% 
  pivot_wider(names_from = donor,values_from = lfc) %>% 
  full_join(.,
          merged %>% dplyr::select(library,name,start,padj,celltype,target,log2FoldChange),
          by=c('library','name','start','celltype','target')) %>% 
  dplyr::mutate(color=ifelse(padj<0.05&log2FoldChange>0,low_col,
                             ifelse(padj<0.05&log2FoldChange<0,high_col,mid_col)))
ggplot(plot, aes(D1,D2,color=I(color), alpha=0.5))+
  geom_point(data=plot %>% dplyr::filter(color==mid_col), shape=16)+
  geom_smooth(method='lm', data=plot %>% dplyr::filter(color==mid_col))+
  geom_point(data=plot %>% dplyr::filter(color!=mid_col), shape=16)+
  geom_smooth(method='lm', data=plot %>% dplyr::filter(color!=mid_col))+
  stat_cor(method = "pearson")

## Distance from TSS Effect
plot = merged %>% 
  group_by(target) %>% rowwise() %>% 
  dplyr::mutate(category=ifelse(abs(start-get(paste0(target,'_start')))<500,'500bp',
                                ifelse(between(abs(start-get(paste0(target,'_start'))),500,1000),'1kb',
                                       ifelse(between(abs(start-get(paste0(target,'_start'))),1000,5000),'5kb','>5kb')))) %>% 
  dplyr::filter(!(celltype=='Treg'&target=='CTLA4'&condition=='0h')) %>% 
  dplyr::mutate(category=factor(category,levels=c('500bp','1kb','5kb','>5kb')))
plot %>% ggplot(aes(category,-log2FoldChange,fill=target))+
  geom_boxplot(outlier.shape = 16,outlier.size = rel(3),outlier.color = 'black',color='black')+
  geom_hline(yintercept = 0)+
  facet_grid(rows=vars(celltype))+
  scale_fill_brewer(palette = 'Greys')+
  theme(legend.position = 'none')+
  NULL
# ANOVA
compare_means(log2FoldChange ~ category,  
              group.by = c('celltype','target','condition'),
              method = 'anova',
              data = plot)
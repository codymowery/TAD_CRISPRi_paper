## Everything in hg38 coordinates 
source('sourceFigures.R')

# Load new Perturb-ATAC
trans = read_rds('tconv_transRegulators_geneSummary_merged.rds') ## Supplementary Table 4

# Load Perturb-ATAC published by Freimer et al 2022
trans_rna = read_csv('FreimerEtAl_Table_S5_RNA_Seq_results.csv')
trans_atac = 
  left_join(read_csv('ATAC_logFC.csv'), 
            read_tsv('ATAC_counts_forCoords.txt') %>% dplyr::select(c('Chr','start','end','width','strand','peakName')),
            by='peakName') %>% 
  dplyr::filter(Chr=='chr2'&(between(start,TAD_start,TAD_stop)|between(end,TAD_start,TAD_stop)))

## Trans Screen Heatmap
data = trans %>% 
  dplyr::filter(id %in% (trans %>% dplyr::filter(fdr<0.05) %>% pull(id))) %>% 
  group_by(screen,id) %>% 
  dplyr::mutate(sig=ifelse(min(fdr)<0.05,screen,NA)) %>% 
  dplyr::select(screen,lfc,id,sig) %>% 
  distinct() %>% 
  pivot_wider(names_from = screen, values_from = c(lfc,sig)) %>% 
  dplyr::mutate(from0 = abs((0-lfc_CD28)+(0-lfc_CTLA4)+(0-lfc_ICOS)),
                sigAll = ifelse(id=='Non-Targeting','Non-Targeting',paste(sig_CD28,sig_CTLA4,sig_ICOS,sep = '_'))) %>% 
  arrange(from0)
mat=as.matrix(data[2:4])
rownames(mat)=as.matrix(data[1])
Heatmap(mat,
        column_order = c('lfc_CD28','lfc_CTLA4','lfc_ICOS'),row_names_side = 'left',
        row_names_gp = gpar(fontsize=5 ),
        heatmap_legend_param = list(title = "Log2FC"),
        col = circlize::colorRamp2(c(-2, 0, 2), c(trans_highCol, trans_midCol, trans_lowCol)),
        row_split = factor(data$sigAll,levels=c("Non-Targeting","CD28_NA_NA","NA_CTLA4_NA","NA_NA_ICOS","CD28_CTLA4_NA","CD28_NA_ICOS","NA_CTLA4_ICOS","CD28_CTLA4_ICOS")),
        cluster_row_slices = FALSE,
        row_gap = unit(3, "mm"),
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 'gray'))),
        row_dend_side = "right",
        name='hm')

## Compare screen results vs arrayed bulk RNAseq validation
# RNAseq validation
plot = trans_atac %>%
  dplyr::filter(sample %in%
                  (trans_rna %>% dplyr::filter(adj.P.Val<0.05&(grepl('CD28|CTLA4|ICOS',gene_name))) %>% pull(sample))) %>%
  group_by(sample) %>%
  dplyr::mutate(mean = mean(log2FoldChange),
                ko=fct_reorder(sample,mean)) %>%
  arrange(mean)
trans_rna %>% dplyr::filter(gene_name %in% c('CD28','CTLA4','ICOS')) %>%
  dplyr::mutate(sample=factor(sample, levels=rev(unique(plot$ko)))) %>% dplyr::filter(!is.na(sample)) %>%
  ggplot(aes(gene_name,sample,fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(midpoint=0,low=rna_lowCol,high=rna_highCol,mid = 'white')

# Screens
genes = str_replace_all(rev(unique(plot$ko)),' KO','')
data = trans %>% 
  dplyr::filter(id %in% genes) %>% 
  dplyr::select(screen,lfc,id) %>% 
  distinct()
data %>% 
  dplyr::mutate(id=factor(id, levels=genes)) %>% dplyr::filter(!is.na(id)) %>%
  ggplot(aes(screen,id,fill=lfc))+
  geom_tile()+
  scale_fill_gradient2(midpoint=0,low=trans_highCol,high=trans_lowCol,mid =trans_midCol)


## Rug plot of individual sgRNAs
regs=c('FOXO1','BPTF','HINFP','MYC','SMARCB1','TFDP1','ZNF217')
guides = bind_rows(
  read_csv('FreimerEtAl_Table_S3_screen_sgrna_results.csv') %>% dplyr::filter(screen=='CTLA4'),
  read_tsv('CD28_low_high_pos_enrichment_2021-08-03.sgrna_summary.txt') %>% dplyr::mutate(screen='CD28'),
  read_tsv('ICOS_low_high_pos_enrichment_2021-11-09.sgrna_summary.txt') %>% dplyr::mutate(screen='ICOS')
) %>% 
  dplyr::mutate(isSig = ifelse(Gene %in% regs,Gene,NA)) %>% 
  dplyr::mutate(isSig = factor(isSig, levels=c('FOXO1','BPTF','HINFP','MYC','SMARCB1','TFDP1','ZNF217'))) %>% 
  dplyr::left_join(.,
                   trans %>% dplyr::select(id,screen,screen_lfc=lfc),
                   by=c('screen','isSig'='id')) %>% 
  dplyr::filter(control_mean>0&treat_mean>0) ## Remove a single guide in ICOS that had 0 HI counts and shifted the whole pop
for(gene in c('CD28','CTLA4','ICOS')){
  sub = guides %>% dplyr::filter(screen==gene)
  p = plot_grid(
    ggplot()+
      geom_rug(data=sub %>% dplyr::filter(is.na(isSig)) %>% dplyr::select(-isSig),
               aes(x=-LFC),alpha=0.01,color='black',length = unit(1, "npc"))+
      geom_rug(data=sub %>% dplyr::filter(!is.na(isSig)),
               aes(x=-LFC,color=I(ifelse(screen_lfc>0,trans_lowCol,trans_highCol))),alpha=1,length = unit(1, "npc"))+
      facet_grid(factor(isSig,c('FOXO1','BPTF','HINFP','MYC','SMARCB1','TFDP1','ZNF217'))~screen)+
      theme(panel.background = element_rect(fill = "grey95"))+
      geom_histogram()+
      NULL,
    ncol = 1, align = 'v', axis = 'lr',rel_heights = c(4))
  print(p)
}

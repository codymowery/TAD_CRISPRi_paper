source('sourceFigures.R')

trans = read_rds('tconv_transRegulators_geneSummary_merged.rds')
trans_rna = read_csv('FreimerEtAl_Table_S5_RNA_Seq_results.csv')
trans_atac = 
  left_join(read_csv('ATAC_logFC.csv'), 
            read_tsv('ATAC_counts_forCoords.txt') %>% 
              dplyr::select(c('Chr','start','end','width','strand','peakName')),
            by='peakName') %>% 
  dplyr::filter(Chr=='chr2'&(between(start,TAD_start,TAD_stop)|between(end,TAD_start,TAD_stop)))

## ZNF217 tracks
# See "ZNF217 Across Locus" in Figures3_S4_TransATAC.R

## ZNF217 KO RNAseq vs Trans Screens
inner_join(
  trans_rna %>% dplyr::filter(sample=='ZNF217 KO'&adj.P.Val<0.1),
  trans %>% dplyr::filter(fdr<0.1) %>% dplyr::select(id,lfc,screen) %>% distinct(),
  by=c('gene_name'='id')) %>% 
  ggplot(aes(logFC,lfc,label=gene_name))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(shape=16,size=rel(3))+
  geom_text_repel(max.overlaps = Inf)+
  facet_grid(cols=vars(screen))

## ZNF217 EnrichR
library(enrichR)
setEnrichrSite("Enrichr")
dbs = c('GO_Biological_Process_2018',
        'KEGG_2021_Human')
up = enrichr(trans_rna %>%
               dplyr::filter(sample=='ZNF217 KO'&adj.P.Val<0.05) %>%
               dplyr::filter(logFC>0) %>% pull(gene_name),
             dbs)
down = enrichr(trans_rna %>%
                 dplyr::filter(sample=='ZNF217 KO'&adj.P.Val<0.05) %>%
                 dplyr::filter(logFC<0) %>% pull(gene_name),
               dbs)
as_tibble(
  bind_rows(head(up$KEGG_2021_Human,3),
            head(up$GO_Biological_Process_2018,3),
            head(down$KEGG_2021_Human,3),
            head(down$GO_Biological_Process_2018,3))) %>% dplyr::mutate(dataset=c(rep('up',6),rep('down',6))) %>% 
  dplyr::arrange(Adjusted.P.value) %>% 
  dplyr::mutate(Term=fct_reorder(Term,-Adjusted.P.value)) %>% 
  ggplot(aes(x=Term,y=-log10(Adjusted.P.value)))+
  coord_flip()+
  geom_bar(stat='identity',fill='black')+
  facet_wrap(~dataset,scales = 'free',drop = T)+
  NULL

# Volcano Plot
sub = trans_rna %>%
  dplyr::filter(sample=='ZNF217 KO') 
up_diff=setdiff(down$GO_Biological_Process_2018$Term[1:3],up$GO_Biological_Process_2018$Term[1:3])
plot = sub %>% 
  dplyr::mutate(kegg=ifelse(gene_name %in% unlist(strsplit(up$KEGG_2021_Human[1,]$Genes,';')), paste0('up:',up$KEGG_2021_Human[1,]$Term),
                            ifelse(gene_name %in% unlist(strsplit(down$KEGG_2021_Human[1,]$Genes,';')), paste0('down:',down$KEGG_2021_Human[1,]$Term),
                                   ifelse(gene_name %in% unlist(strsplit(up$GO_Biological_Process_2018[1,]$Genes,';')), paste0('up:',up$GO_Biological_Process_2018[1,]$Term),
                                          ifelse(gene_name %in% unlist(strsplit(down$GO_Biological_Process_2018[1,]$Genes,';')), paste0('down:',down$GO_Biological_Process_2018[1,]$Term),1))))) %>% 
  dplyr::mutate(label=ifelse(gene_name %in% 
                               c('ZNF217',
                                 (sub %>% dplyr::filter(kegg!='1'&logFC<0) %>% slice_min(order_by=`adj.P.Val`,n=10) %>% pull(gene_name)),
                                 (sub %>% dplyr::filter(kegg!='1'&logFC>0) %>% slice_min(order_by=`adj.P.Val`,n=10) %>% pull(gene_name))),
                             gene_name,NA)) 
col_darken = c(0,0.5)
darker = function(col,x) {colorspace::darken(col,x)}
ggplot(data=plot %>% dplyr::arrange(kegg),
       aes(logFC,-log10(adj.P.Val),color=kegg,label=label,size=I(ifelse(!is.na(kegg),rel(3),rel(1)))),shape=16)+
  geom_point(shape=16,size=rel(3))+
  geom_hline(yintercept = -log10(0.05),linetype='dashed',color='black')+
  geom_text_repel(color='black',box.padding = 0.5,max.overlaps = Inf)+
  scale_color_manual(values=c('grey90',
                              sapply(col_darken, darker, col=wesanderson::wes_palette('Zissou1')[1]),
                              sapply(col_darken, darker, col=wesanderson::wes_palette('Zissou1')[3])))

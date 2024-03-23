source('sourceFigures.R')

trans = read_rds('tconv_transRegulators_geneSummary_merged.rds')
trans_rna = read_csv('FreimerEtAl_Table_S5_RNA_Seq_results.csv')
trans_atac = 
  left_join(read_csv('ATAC_logFC.csv'), 
            read_tsv('ATAC_counts_forCoords.txt') %>% 
              dplyr::select(c('Chr','start','end','width','strand','peakName')),
            by='peakName') %>% 
  dplyr::filter(Chr=='chr2'&(between(start,TAD_start,TAD_stop)|between(end,TAD_start,TAD_stop)))

## Trans ATAC of all 3 target genes --------
for(set in list(c('CD28',5e3,5e3,2),
                c('CTLA4',40e3,2e3,3),
                c('ICOS',5e3,0,2))){
  gene=set[1]
  upstream=get(paste0(gene,'_start'))-as.numeric(set[2])
  downstream=get(paste0(gene,'_stop'))+as.numeric(set[3])
  region=c(upstream,downstream)
  rel_height=set[4]
  plot = trans_atac %>% dplyr::filter(between(start,region[1],region[2])) %>% 
    dplyr::filter((sample %in% 
                     (trans_rna %>% dplyr::filter(adj.P.Val<0.05&(grepl('CD28|CTLA4|ICOS',gene_name))) %>% pull(sample)))|
                    (sample %in% (trans %>% dplyr::filter(fdr<0.05) %>% dplyr::mutate(sample = paste0(id,' KO')) %>% pull(sample)))) %>% 
    dplyr::filter(sample %in% (trans_atac %>% dplyr::filter(between(start,region[1],region[2])) %>% dplyr::filter(padj<0.05) %>% pull(sample))) %>% ## Removing samples wo any significant peaks
    group_by(sample) %>% 
    dplyr::mutate(mean = mean(log2FoldChange),
                  ko=fct_reorder(sample,mean)) %>% 
    arrange(mean)
  p = plot_grid(
    merged %>% dplyr::filter(target==gene&between(start,region[1],region[2])) %>%
      dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
      dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
      dplyr::mutate(celltype=factor(celltype,levels=c('Tconv','Treg'))) %>%
      ggplot()+
      annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
      annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
      annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
      geom_hline(yintercept = 0, linetype='solid', color='black')+
      geom_point(aes(start,-log10(pvalue),color=pval_ceiling_dir),shape=16,size=rel(3),alpha=0.75)+
      scale_color_gradient2(low=low_col, midpoint=0, mid=mid_col, high=high_col)+
      facet_grid(rows=vars(target,condition,celltype),switch = 'y',space='fixed', scales = 'free_y')+
      coord_cartesian(xlim = region, expand=FALSE, clip = "off")+
      theme_classic()+
      theme(strip.placement = 'outside', strip.background = element_rect(color='white',fill='white'),
            axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
            axis.title = element_blank(),
            legend.position = 'none',
            panel.spacing = unit(0.7, "lines"))+
      NULL,
    ggplot()+
      geom_rect(data=plot %>% arrange(abs(log2FoldChange)),
                aes(xmin=start,xmax=end,ymin=0,ymax=baseMean,fill=log2FoldChange,color=I(ifelse(padj<0.1,'black',NA))))+ ## Up/Down type plot
      facet_grid(factor(ko, levels=unique(plot$ko))~.,switch = 'y')+
      scale_fill_gradient2(low=trans_lowCol, midpoint=0, mid=trans_midCol, high=trans_highCol)+
      coord_cartesian(xlim = region, expand=F)+
      theme_void()+
      theme(strip.text.y.left = element_text(angle = 0),
            panel.spacing = unit(0, "lines")),
    ncol=1, align='v',axis='lr', rel_heights = c(as.numeric(rel_height),6,1,1))
  print(p)
}

# Tconv ATAC
for(set in list(c('CD28',5e3,5e3,2),
                c('CTLA4',40e3,2e3,3),
                c('ICOS',5e3,0,2))){
  gene=set[1]
  upstream=get(paste0(gene,'_start'))-as.numeric(set[2])
  downstream=get(paste0(gene,'_stop'))+as.numeric(set[3])
  region=c(upstream,downstream)
  max = max(as_tibble(import.bw('CTRL_coverage_at_tad.bw')) %>% dplyr::mutate(sample='AAVS1') %>% 
    dplyr::filter(between(start,region[1],region[2])) %>% pull(score))
  p = as_tibble(import.bw('CTRL_coverage_at_tad.bw')) %>% dplyr::mutate(sample='AAVS1') %>% 
    dplyr::filter(between(start,region[1],region[2])) %>% 
    ggplot()+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    geom_line(aes(x=start+(end-start)/2,y=score))+
    coord_cartesian(xlim = region, expand=FALSE)+
    theme_classic()+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.line.x = element_blank(), axis.ticks.x = element_blank())+
    labs(title=max)+
    NULL
  print(p)
}

## ID significant regulators to bold ATAC row labels -----------
full_join(trans_rna %>% dplyr::filter(gene_name=='CTLA4'&adj.P.Val<0.05),
          trans %>% dplyr::filter(screen=='CTLA4'&fdr<0.05) %>% dplyr::mutate(sample=paste0(id,' KO')),
          by='sample') %>% 
  dplyr::filter(sample%in%(trans_atac %>% pull(sample)))

full_join(trans_rna %>% dplyr::filter(gene_name=='CD28'&adj.P.Val<0.05),
          trans %>% dplyr::filter(screen=='CD28'&fdr<0.05) %>% dplyr::mutate(sample=paste0(id,' KO')),
          by='sample') %>% 
  dplyr::filter(sample%in%(trans_atac %>% pull(sample)))
full_join(trans_rna %>% dplyr::filter(gene_name=='ICOS'&adj.P.Val<0.05),
          trans %>% dplyr::filter(screen=='ICOS'&fdr<0.05) %>% dplyr::mutate(sample=paste0(id,' KO')),
          by='sample') %>% 
  dplyr::filter(sample%in%(trans_atac %>% pull(sample)))

## KO ATAC Tracks
# Get normalized ATAC files via normATACCoverage_TAD.R

## ZNF217 & IRF4 KO ATAC Track
region=c(CTLA4_start-38.5e3,CTLA4_start-35e3)
atac = tibble() 
samples=tibble(sample=c('CTRL','ZNF217','IRF4'),col=c('black',trans_highCol,trans_lowCol))
for(fname in list.files(path='tadCoverage/',pattern = paste(samples$sample,collapse = '|'))){
  atac = bind_rows(
    atac,
    as_tibble(import.bw(paste0('tadCoverage/',fname))) %>% 
      dplyr::mutate(sample=str_replace(fname,'_coverage_at_tad.bw','')) %>% 
      dplyr::filter(between(start,region[1],region[2])))
}
atac = bind_rows(atac %>% dplyr::select(-c(seqnames,width,strand)),
                 tibble(start=c(rep(region[1]-1,nrow(samples)),rep(region[2]+1,nrow(samples))),
                        end=c(rep(region[1]-1,nrow(samples)),rep(region[2]+1,nrow(samples))),
                        sample=rep(samples$sample,2),
                        score=rep(0,nrow(samples)*2)))
max(atac$score)
ggplot(atac)+
  geom_line(aes(x=start,y=score,color=sample))+
  annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
  annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
  annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
  coord_cartesian(xlim=region,expand = F)+
  theme_classic()+
  scale_color_manual(values=setNames(samples$col,samples$sample))+
  theme(strip.placement = 'outside', strip.background = element_rect(color='white'),legend.position = 'none',
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank())+
  NULL

## ZNF217 Across Locus
region=c(CD28_start-10e3,TAD_stop)
atac = tibble()
samples=tibble(sample=c('CTRL','ZNF217'),col=c('black',trans_highCol))
for(fname in list.files(path='tadCoverage/',pattern = paste(samples$sample,collapse = '|'))){
  atac = bind_rows(
    atac,
    as_tibble(import.bw(paste0('tadCoverage/',fname))) %>% 
      dplyr::mutate(sample=str_replace(fname,'_coverage_at_tad.bw','')) %>% 
      dplyr::filter(between(start,region[1],region[2])))
}
atac = bind_rows(atac %>% dplyr::select(-c(seqnames,width,strand)),
                 tibble(start=c(rep(region[1]-1,nrow(samples)),rep(region[2]+1,nrow(samples))),
                        end=c(rep(region[1]-1,nrow(samples)),rep(region[2]+1,nrow(samples))),
                        sample=rep(samples$sample,2),
                        score=rep(0,nrow(samples)*2)))
stim_atac = makeGRangesFromDataFrame(read_tsv('~/Downloads/41588_2019_505_MOESM6_ESM') %>% tidyr::separate(peak_id, into=c('CHR','START','END'), '_'),
                                     keep.extra.columns = T)
hg19_hg38_chain = import.chain('hg19ToHg38.over.chain')
stim_atac_hg38 = as_tibble(unlist(liftOver(stim_atac, hg19_hg38_chain)))

plot_grid(
  ggplot(atac)+
    geom_line(aes(x=start,y=score,color=sample))+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    coord_cartesian(xlim=region,expand = F)+
    theme_classic()+
    scale_color_manual(values=setNames(samples$col,samples$sample))+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white'),legend.position = 'none',
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.line.x = element_blank(), axis.ticks.x = element_blank())+
    NULL,
  trans_atac %>% dplyr::filter(sample=='ZNF217 KO'&padj<0.05&between(start,region[1],region[2])) %>% 
    ggplot()+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    ggplot2::geom_rect(aes(xmin=start,xmax=end,ymin=-1,ymax=1,fill=log2FoldChange))+
    coord_cartesian(xlim=region,expand = F)+
    theme_classic()+
    scale_fill_gradient2(low=trans_lowCol, midpoint=0, mid=trans_midCol, high=trans_highCol)+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white'),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.line.x = element_blank(), axis.ticks.x = element_blank())+
    NULL,
  stim_atac_hg38 %>% dplyr::filter(adj.P.Val<0.05&seqnames=='chr2'&between(start,region[1],region[2])) %>% 
    ggplot()+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color='black',linetype='dashed',fill=NA)+
    ggplot2::geom_rect(aes(xmin=start,xmax=end,ymin=-1,ymax=1,fill=logFC))+
    coord_cartesian(xlim=region,expand = F)+
    theme_classic()+
    scale_fill_gradient2(low='#2d004b', midpoint=0, mid='white', high='#7f3b08')+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white'),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.line.x = element_blank(), axis.ticks.x = element_blank())+
    NULL,
  ncol=1, align='v', axis='lr', rel_heights = c(3,0.5,0.5))

## ZNF217 KO Heatmap for schematic
samples = c('ZNF217')
plot = trans_rna %>% dplyr::filter(sample %in% paste0(samples,' KO') & gene_name %in% c('CTLA4','IRF4','ZNF217'))
ggplot(plot, aes(sample,gene_name,fill=logFC))+
  geom_tile(aes(color=I(ifelse(adj.P.Val<0.05,'black',NA))),size=2)+
  scale_fill_gradient2(low=rna_lowCol, midpoint=0, mid=rna_midCol, high=rna_highCol,
                       limits=c(
                         -max(abs(plot$logFC)),
                         max(abs(plot$logFC))))

## JASPAR 22 motifs
motifs = read_tsv('JASPAR22_motifs_chr2_202527032_203967032.txt', col_names = F) %>% 
  dplyr::filter(X5>400) %>% ## Default JASPAR 22 filter in UCSC
  dplyr::filter(grepl('IRF1$|MYB$',X7,ignore.case = T))

plot_region = c(203827771,203837771) ## CTLA4 enhancers, same coordinates as ChIP

ggplot(motifs %>% dplyr::filter(between(X2,plot_region[1],plot_region[2]) & !grepl('::',X7)))+
  annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
  annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
  annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
  geom_rect(aes(xmin=X2, xmax=X3, ymin=0, ymax=1), color='black')+
  coord_cartesian(xlim = plot_region, expand=FALSE)+
  facet_grid(rows=vars(X7))+
  theme_void()+
  NULL
## Everything in hg38 coordinates 
source('sourceFigures.R')

## CRISPRi tracks for CD28 and ICOS
for(set in list(c('CD28',5e3,5e3,2),
                c('ICOS',5e3,0,2))){
  gene=set[1]
  upstream=get(paste0(gene,'_start'))-as.numeric(set[2])
  downstream=get(paste0(gene,'_stop'))+as.numeric(set[3])
  region=c(upstream,downstream)
  rel_height=set[4]
  p = plot_grid(
    merged %>% dplyr::filter(target==gene&between(start,region[1],region[2])) %>%
      dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
      dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
      dplyr::mutate(celltype=factor(celltype,levels=c('Treg','Tconv'))) %>%
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
    ncol=1, align='v',axis='lr', rel_heights = c(as.numeric(rel_height),6,1,1))
  print(p)
}

## All signals at Stim-Responsive CiRE
region=c(203829159,203831178)
p = plot_grid(
  merged %>% dplyr::filter(between(start,region[1],region[2])) %>%
    dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
    dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
    dplyr::mutate(celltype=factor(celltype,levels=c('Treg','Tconv'))) %>%
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
          axis.title = element_blank(),
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  ncol=1, align='v',axis='lr', rel_heights = c(as.numeric(rel_height),6,1,1))
print(p)

## All signals at Treg-Dominant CiRE
region=c(203835010,203836851)
p = plot_grid(
  merged %>% dplyr::filter(between(start,region[1],region[2])) %>%
    dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
    dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
    dplyr::mutate(celltype=factor(celltype,levels=c('Treg','Tconv'))) %>%
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
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  ncol=1, align='v',axis='lr', rel_heights = c(as.numeric(rel_height),6,1,1))
print(p)

## Histogram plots for sgRNA binning in Figure2_CTLA4Enhancers.R
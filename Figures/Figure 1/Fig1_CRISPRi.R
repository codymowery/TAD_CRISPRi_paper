## Everything in hg38 coordinates 

source('sourceFigures.R')

# HiC downloaded from http://3dgenome.fsm.northwestern.edu/ using K562, Rao_2014-raw, 10kb res, approx coordinates chr2:201000000–205500000 

## Gene Model for assay schematic in Fig 1a
ggplot()+
  ggplot2::geom_segment(aes(x=202520000, xend=203970000, y=0, yend=0))+
  ggplot2::geom_rect(aes(xmin=CD28_start, xmax=CD28_stop, ymin=-1, ymax=1))+
  ggplot2::geom_rect(aes(xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-1, ymax=1))+
  ggplot2::geom_rect(aes(xmin=ICOS_start, xmax=ICOS_stop, ymin=-1, ymax=1))+
  theme_void()

## Plot CRISPRi
region=c(202500000,203980000) ## HiC plot coordinates
padj_thresh=0.05
plot_grid(
  full_geneRef(range=region),
  NULL,
  merged %>% dplyr::filter(!(target=='CTLA4'&condition=='0h')&between(start,region[1],region[2])) %>%
    dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
    dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
    dplyr::mutate(col=ifelse(log2FoldChange>0&padj<padj_thresh,low_col,
                             ifelse(log2FoldChange<0&padj<padj_thresh,high_col,mid_col)),
                  celltype=factor(celltype,levels=c('Tconv','Treg'))) %>% 
    ggplot()+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    geom_hline(yintercept = 0, linetype='solid', color='black')+
    geom_point(aes(start,-log10(pvalue),color=I(col)),shape=16,size=rel(2),alpha=0.6)+
    facet_grid(rows=vars(target,condition,celltype),switch = 'y',space='fixed', scales = 'free_y')+
    coord_cartesian(xlim = region, expand=FALSE, clip = "off")+
    theme(strip.background = element_blank(),
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  ncol=1, align='v',axis='lr', rel_heights = c(1,0.25,6))

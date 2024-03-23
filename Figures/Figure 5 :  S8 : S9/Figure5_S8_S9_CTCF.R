## Everything in hg38 coordinates 

source('sourceFigures.R')
library(BSgenome.Hsapiens.UCSC.hg38)

###
## Figure 5
###
## Load files -------
# GSM5668095 in hg38 coords
loops = read_tsv('GSM5668095_ENCFF116POO_long_range_chromatin_interactions_GRCh38.bedpe',col_names = F) %>%
  rowwise() %>%
  dplyr::mutate(min=min(X2+(X3-X2),X5+(X6-X5)),max=max(X2+(X3-X2),X5+(X6-X5))) %>%
  dplyr::filter(X1=='chr2'&X4=='chr2')
tracks = as_tibble(import.bw('GSM5668095_ENCFF371IIJ_signal_of_unique_reads_GRCh38.bigWig')) %>%
  dplyr::filter(seqnames=='chr2'&between(start,TAD_start,TAD_stop)&between(end,TAD_start,TAD_stop))
peaks = read_tsv('GSM5668095_ENCFF629TKQ_peaks_GRCh38.bed',col_names = F) %>%
  dplyr::filter(X1=='chr2')
write_rds(list(loops,tracks,peaks),'chiapet_cd4.rds')
# GSM5668096 in hg38 coords
loops = read_tsv('GSM5668096_ENCFF085MVS_long_range_chromatin_interactions_GRCh38.bedpe.gz',col_names = F) %>%
  rowwise() %>%
  dplyr::mutate(min=min(X2+(X3-X2),X5+(X6-X5)),max=max(X2+(X3-X2),X5+(X6-X5))) %>%
  dplyr::filter(X1=='chr2'&X4=='chr2')
tracks = as_tibble(import.bw('GSM5668096_ENCFF678SVN_signal_of_unique_reads_GRCh38.bigWig')) %>%
  dplyr::filter(seqnames=='chr2'&between(start,TAD_start,TAD_stop)&between(end,TAD_start,TAD_stop))
peaks = read_tsv('GSM5668096_ENCFF776EHF_peaks_GRCh38.bed.gz',col_names = F) %>%
  dplyr::filter(X1=='chr2')
write_rds(list(loops,tracks,peaks),'chiapet_GSM5668096_cd4.rds')

########## 2 Donors ----------
loops=bind_rows(
  read_rds('chiapet_cd4.rds')[[1]] %>% dplyr::mutate(donor='GSM5668095'),
  read_rds('chiapet_GSM5668096_cd4.rds')[[1]] %>% dplyr::mutate(donor='GSM5668096'))
tracks=bind_rows(
  read_rds('chiapet_cd4.rds')[[2]] %>% dplyr::mutate(donor='GSM5668095'),
  read_rds('chiapet_GSM5668096_cd4.rds')[[2]] %>% dplyr::mutate(donor='GSM5668096'))
peaks=bind_rows(
  read_rds('chiapet_cd4.rds')[[3]] %>% dplyr::mutate(donor='GSM5668095'),
  read_rds('chiapet_GSM5668096_cd4.rds')[[3]] %>% dplyr::mutate(donor='GSM5668096'))

target_region=c(CD28_start-10e3,TAD_stop)
plot_region=c(CD28_start-10e3,TAD_stop)

## Filter for target region
loops2 = loops %>% 
  dplyr::filter(between(max,target_region[1],target_region[2])&between(min,target_region[1],target_region[2])) %>% 
  dplyr::mutate(midpoint=min+(max-min)/2,
                width=max-min) %>% 
  group_by(donor) %>% 
  dplyr::mutate(rangeFraction=width/max(.$width)) 
## Filter by loops shared between donors within 5k of one another
intersect = 
  distance_inner_join(
    loops2 %>% dplyr::filter(donor=='GSM5668095'),
    loops2 %>% dplyr::filter(donor=='GSM5668096'),
    by=c('min','max'),
    max_dist=5000) %>% 
  pivot_longer(cols = contains('.'), names_to = c('.value','donor_var'), names_sep = '\\.')
## Plot over Tconv and Treg screen
plot_grid(
  merged %>% dplyr::filter(between(start,plot_region[1],plot_region[2])) %>%
    dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
    dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
    dplyr::mutate(celltype=factor(celltype,levels=c('Tconv','Treg'))) %>% 
    ggplot()+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    geom_hline(yintercept = 0, linetype='solid', color='black')+
    geom_point(aes(start,-log10(pvalue),color=pval_ceiling_dir),shape=16,size=rel(2),alpha=0.75)+
    scale_color_gradient2(low=low_col, midpoint=0, mid=mid_col, high=high_col)+
    facet_grid(rows=vars(target,celltype,condition),switch = 'y',space='fixed', scales = 'free_y')+
    coord_cartesian(xlim = plot_region, expand=FALSE, clip = "off")+
    theme_classic()+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white',fill='white'),
          axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  NULL,
  ggplot(intersect)+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    geom_curve(aes(x = min, y = 0, xend = max, yend = 0), curvature = 0.35)+
    coord_cartesian(xlim = plot_region,ylim=c(-2,0), expand=FALSE)+
    theme_void()+
    NULL,
  ncol = 1, align = 'v', axis = 'lr', rel_heights = c(6,1,3))
plot_grid(
  ggplot(tracks %>% dplyr::filter(between(start,plot_region[1],plot_region[2])))+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    geom_line(aes(start+(end-start)/2,score))+
    coord_cartesian(xlim = plot_region, expand=FALSE)+
    theme_void()+
    NULL,
  ncol = 1, align = 'v', axis = 'lr', rel_heights = c(2))
max(tracks %>% dplyr::filter(between(start,plot_region[1],plot_region[2])) %>% pull(score))

## 1D Guide Plot for 4C track
plot_region=c(203494157,203972174)
plot = merged %>% dplyr::filter(celltype=='Tconv') %>% 
  dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
  dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling))

plot_grid(
  full_geneRef(range=plot_region),
  ggplot()+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5,)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    coord_cartesian(xlim = plot_region, ylim=c(0,1), expand=FALSE, clip = "off")+
    theme_void(),
  plot %>% dplyr::filter(between(start,plot_region[1],plot_region[2])) %>% arrange(desc(padj),start) %>% 
    ggplot()+
    geom_hline(yintercept = 0, linetype='solid', color='black')+
    geom_point(aes(start,1,color=pval_ceiling_dir),shape=15,size=rel(1),alpha=0.75)+
    scale_color_gradient2(low=low_col, midpoint=0, mid='grey95', high=high_col)+
    facet_grid(rows=vars(target))+
    coord_cartesian(xlim = plot_region, ylim=c(1,1), expand=FALSE, clip = "off")+
    theme_classic(),
  ncol = 1, align = 'v', axis = 'lr',rel_heights = c(1,0.25,3))

####
## Figure S8
#### 

# ## CRISPRi Validation
combine = read_rds('CTCF_crispri_combinedExperiments_unprocessed.rds') %>% 
  dplyr::filter(Count>500) %>%
  group_by(exp,Donor,flowTarget,celltype,Condition,Plate) %>% 
  dplyr::mutate(ctrl_median_recalc = median(MFI[sgRNA=='NTC'], na.rm=T)) %>%
  dplyr::mutate(fc_recalc = MFI/ctrl_median_recalc)
# Plot
combine %>% 
  dplyr::mutate(sgRNA = factor(sgRNA, levels=c('NTC','5p_CTCF'))) %>% 
  ggplot(aes(sgRNA,fc_recalc,fill=sgRNA))+
  geom_boxplot(color='black', outlier.shape = NA)+
  geom_point(color='black', shape = 16, alpha = 0.5)+
  facet_grid(rows=vars(flowTarget),
             cols=vars(celltype,Condition),
             scales='free_y')+
  scale_fill_manual(values = c('#000000','#00a087'))+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle=60,hjust=1),
        panel.spacing = unit(0.5, "lines"))
combine %>% group_by(celltype,Condition,sgRNA,flowTarget) %>% dplyr::count() %>% view()
compare_means(data=combine %>% dplyr::filter(sgRNA %in% c('5p_CTCF','NTC') & Condition!='24h') %>% ungroup(),
              formula = fc_recalc ~ sgRNA,
              ref.group = 'NTC',
              group.by = c('flowTarget','celltype','Condition'),
              method='t.test')

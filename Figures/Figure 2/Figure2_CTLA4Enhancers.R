## Everything in hg38 coordinates 

source('sourceFigures.R')

## File processing
## eQTL data from DOI 10.1126/science.abf1970  ------
eqtl = read_tsv('Table6_lupusPaper_gracieGordon.txt', skip = 1) %>%
  dplyr::filter((grepl('CD28|CTLA4|ICOS',RSID))) %>% dplyr::filter(!grepl('ICOSLG',RSID)) %>%
  tidyr::separate(RSID, into=c('coord','gene'), '_') %>% tidyr::separate(coord, into=c('chr','start'), ':') %>%
  dplyr::mutate(start=as.double(start))
pre = eqtl %>%
  dplyr::mutate(seqnames='chr2',end=start) %>%
  dplyr::relocate(seqnames,start,end) %>% dplyr::select(-chr)
lifted = as_tibble(unlist(liftOver(makeGRangesFromDataFrame(pre,keep.extra.columns = T), hg19_hg38_chain)))
write_rds(lifted, 'Table6_lupusPaper_gracieGordon_TAD_lifted.txt')

## Download RA GWAS summary stats from http://plaza.umin.ac.jp/~yokada/datasource/software.htm ---------
ra = bind_rows(
  read_tsv('~/Downloads/RA_GWASmeta_TransEthnic_v2.txt.gz') %>% dplyr::filter(Chr==2&between(`Position(hg19)`,203391755,204831755)) %>% dplyr::mutate(study='Transethnic'),
  read_tsv('~/Downloads/RA_GWASmeta_European_v2.txt.gz') %>% dplyr::filter(Chr==2&between(`Position(hg19)`,203391755,204831755)) %>% dplyr::mutate(study='European'),
  read_tsv('~/Downloads/RA_GWASmeta_Asian_v2.txt.gz') %>% dplyr::filter(Chr==2&between(`Position(hg19)`,203391755,204831755)) %>% dplyr::mutate(study='Asian'))  ## TAD hg19 coords
ra_hg38 = as_tibble(unlist(liftOver(makeGRangesFromDataFrame(ra %>% dplyr::mutate(Chr='chr2'),keep.extra.columns = T,start.field = 'Position(hg19)',end.field='Position(hg19)'),chain = hg19_hg38_chain)))
write_rds(ra_hg38,'RA_GWASmeta_hg38lifted.rds')
ra_hg38 = read_rds('RA_GWASmeta_hg38lifted.rds')
region=c(CTLA4_start-42e3,CTLA4_stop+5e3)
top = ra_hg38 %>% dplyr::filter(between(start,region[1],region[2])) %>% slice_min(n=1,P.val) %>% pull(SNPID)
top_ld = as_tibble(LDproxy(snp = top,
                           pop = "ALL",
                           r2d = "r2",
                           token = 'adae89590394')) %>%
  tidyr::separate(Coord, into=c('chr','start'), ':') %>%
  dplyr::mutate(ref_snp = top,start=as.double(start))
lifted_topRA_ld = as_tibble(unlist(liftOver(makeGRangesFromDataFrame(top_ld,keep.extra.columns = T,seqnames.field = 'chr',start.field = 'start',end.field='start'),chain = hg19_hg38_chain)))
write_rds(left_join(ra_hg38,lifted_topRA_ld, by=c('start')),'RA_GWASmeta_hg38lifted_wLD.rds')

## Lupus LD ----------
lupus_eqtl = read_rds('Table6_lupusPaper_gracieGordon_TAD_lifted.txt')
top_lupus = lupus_eqtl %>% dplyr::filter(between(start,region[1],region[2])) %>% slice_min(n = 1,PVALUE_FE) %>% pull(start)
top_lupus_rs = 'rs3087243' ## searched for 2:start on https://www.snpedia.com/
top_lupus_ld = as_tibble(LDproxy(snp = top_lupus_rs,
                                 pop = "ALL",
                                 r2d = "r2",
                                 token = 'adae89590394')) %>%
  tidyr::separate(Coord, into=c('chr','start'), ':') %>%
  dplyr::mutate(ref_snp = top_lupus_rs,start=as.double(start))
lifted_top_lupus_ld = as_tibble(unlist(liftOver(makeGRangesFromDataFrame(top_lupus_ld,keep.extra.columns = T,seqnames.field = 'chr',start.field = 'start',end.field='start'),chain = hg19_hg38_chain)))
write_rds(lifted_top_lupus_ld,'Table6_lupusPaper_gracieGordon_TAD_lifted_wLD.rds')

## Plot -------
lupus_eqtl = read_rds('Table6_lupusPaper_gracieGordon_TAD_lifted.txt')
lifted_top_lupus_ld = read_rds('Table6_lupusPaper_gracieGordon_TAD_lifted_wLD.rds')
ra_hg38_ld = read_rds('RA_GWASmeta_hg38lifted_wLD.rds')

## Identify peaks of CRISPRi signal
region=c(CTLA4_start-42e3,CTLA4_stop+5e3)
sig_df = merged %>% dplyr::filter(target=='CTLA4'&padj<0.05) %>% arrange(start) %>% dplyr::select(name, start) %>% distinct()
sig_df$distanceToNextSig = 0
i = 1 
while(i<nrow(sig_df)){
  sig_df$distanceToNextSig[i]=sig_df$start[i+1]-sig_df$start[i]
  i = i+1
}
sig_df
ggplot(sig_df, aes(distanceToNextSig))+geom_histogram(fill='black')+scale_x_log10()+labs(y='Number of sgRNAs',x='Bp to Next Significant sgRNA')
ggplot(sig_df %>% dplyr::filter(distanceToNextSig<1000), aes(distanceToNextSig))+geom_histogram(fill='black')+scale_x_continuous(limits=c(0,1000))+labs(y='Number of sgRNAs',x='Bp to Next Significant sgRNA')
# Add new "peak" annotation 
i = 1
j = 1
sig_df$peak = NA
while(i<nrow(sig_df)){
  sig_df$peak[i]=paste0('peak_',as.character(j))
  if(sig_df$distanceToNextSig[i]>500){j=j+1}
  i = i+1
}
peakPlot = sig_df %>% group_by(peak) %>% add_count() %>% dplyr::filter(n>1) %>% dplyr::mutate(peakMin=min(start),peakMax=max(start)) %>% 
  dplyr::select(-c(name,start,distanceToNextSig)) %>% distinct() %>% ungroup()
ggplot(peakPlot %>% mutate(x=NA), aes(x,n))+geom_boxplot()+geom_point()+labs(x=NULL,y='Number of Significant sgRNAs per Peak')+scale_y_continuous(n.breaks=8)
peakPlot %>% dplyr::mutate(chr=2) %>% dplyr::select(chr,start=peakMin,stop=peakMax) %>% 
  dplyr::filter(between(start,region[1],region[2])&between(stop,region[1],region[2])) %>% 
  write_tsv('ctla4_tilingPeaks_chr2_203825771_203878965.tsv')

## Plot CRISPRi tracks
plot = merged %>% dplyr::filter(target=='CTLA4') %>% 
  dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
  dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
  dplyr::mutate(celltype=factor(celltype,levels=c('Tconv','Treg')))
plot_grid(
  plot %>% dplyr::filter(celltype=='Tconv'&condition=='6h'&between(start,region[1],region[2])) %>% 
    ggplot()+
    geom_rect(data=peakPlot,aes(xmin=peakMin,xmax=peakMax,ymin=-Inf,ymax=Inf),fill='#d6c4a8',alpha=0.3)+
    geom_hline(yintercept = 0, linetype='solid', color='black')+
    geom_point(aes(start,-log10(pvalue),color=pval_ceiling_dir),shape=16,size=rel(3),alpha=0.75)+
    scale_color_gradient2(low=low_col, midpoint=0, mid=mid_col, high=high_col)+
    facet_grid(rows=vars(target,condition,celltype),switch = 'y',space='fixed', scales = 'free_y')+
    coord_cartesian(xlim = region, expand=FALSE, clip = "off")+
    theme_classic()+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white',fill='white'),
          axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
          axis.title = element_blank(),
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  plot %>% dplyr::filter(celltype=='Treg'&condition=='0h'&between(start,region[1],region[2])) %>% 
    ggplot()+
    geom_rect(data=peakPlot,aes(xmin=peakMin,xmax=peakMax,ymin=-Inf,ymax=Inf),fill='#d6c4a8',alpha=0.3)+
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
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  plot %>% dplyr::filter(celltype=='Treg'&condition=='6h'&between(start,region[1],region[2])) %>% 
    ggplot()+
    geom_rect(data=peakPlot,aes(xmin=peakMin,xmax=peakMax,ymin=-Inf,ymax=Inf),fill='#d6c4a8',alpha=0.3)+
    geom_hline(yintercept = 0, linetype='solid', color='black')+
    geom_point(aes(start,-log10(pvalue),color=pval_ceiling_dir),shape=16,size=rel(3),alpha=0.75)+
    scale_color_gradient2(low=low_col, midpoint=0, mid=mid_col, high=high_col)+
    facet_grid(rows=vars(target,condition,celltype),switch = 'y',space='fixed', scales = 'free_y')+
    coord_cartesian(xlim = region, expand=FALSE, clip = "off")+
    theme_classic()+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white',fill='white'),
          axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
          axis.title = element_blank(),
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  ra_hg38_ld %>% dplyr::filter(between(start,region[1],region[2])&study=='Transethnic') %>% 
    ggplot()+
    geom_rect(data=peakPlot,aes(xmin=peakMin,xmax=peakMax,ymin=-Inf,ymax=Inf),fill='#d6c4a8',alpha=0.3)+
    geom_hline(yintercept = -log10(5*10^(-8)), linetype='dashed', color='black')+
    geom_point(aes(start,-log10(P.val),color=R2,shape=I(ifelse(R2==1,17,16))),size=rel(2),alpha=1)+
    theme_classic()+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),strip.placement = 'none')+
    facet_grid(rows=vars(study))+
    scale_color_gradient(low='#d6e9f3',high='black')+
    coord_cartesian(xlim = region, expand=FALSE, clip = 'off')+
    labs(y='RA GWAS\nLog10(Transethnic P-Value)'),
  left_join(lupus_eqtl,lifted_top_lupus_ld, by=c('start')) %>% dplyr::filter(between(start,region[1],region[2]) & cell=='t4') %>% 
    ggplot()+
    geom_rect(data=peakPlot,aes(xmin=peakMin,xmax=peakMax,ymin=-Inf,ymax=Inf),fill='#d6c4a8',alpha=0.3)+
    geom_hline(yintercept = -log10(5*10^(-8)), linetype='dashed', color='black')+
    geom_point(aes(start,-log10(PVALUE_FE),color=R2,shape=I(ifelse(R2==1,17,16))),size=rel(2),alpha=1)+
    theme_classic()+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),strip.placement = 'none')+
    scale_color_gradient(low='#d6e9f3',high='black')+
    coord_cartesian(xlim = region, expand=FALSE, clip = "off")+
    labs(y='Lupus sc-eQTL\nLog10(P-Value)'),
  ncol = 1, align = 'v', axis='lr', rel_heights = c(1,1,1,1,1))

## CRISPRi Validation
plot = read_rds('220319-1_merged_unprocessed.rds') %>% 
  dplyr::filter(count>500) %>% 
  group_by(celltype,donor,flowTarget) %>% 
  dplyr::mutate(ctrl_median = median(MFI[grepl('NTC',sgrnaTarget)], na.rm=T)) %>%
  dplyr::mutate(fc = MFI/ctrl_median)
plot %>% 
  ggplot(aes(sgrnaTarget,fc,fill=I(ifelse(sgrnaTarget=='NTC',NA,low_col))))+
  geom_boxplot(color='black',outlier.shape = NA)+
  geom_point(color='black',alpha=0.5,shape=16,size=rel(2))+
  facet_grid(rows=vars(celltype))+
  geom_hline(yintercept = 1)+
  scale_y_continuous(limits = c(0,1.75))+
  theme(axis.text.x = element_text(angle=60,hjust=1))
plot %>% group_by(celltype, sgrnaTarget,flowTarget,donor) %>% count()
compare_means(fc ~ sgrnaTarget,
              data=plot,
              group.by = 'celltype',
              method='t.test',
              ref.group = "NTC")
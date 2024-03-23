## Everything in hg38 coordinates 
source('sourceFigures.R')

trans = read_rds('tconv_transRegulators_geneSummary_merged.rds')

## Cis co-regulation
region1=c(CD28_start-5e3,CD28_start+5e3)
region2=c(CTLA4_start-5e3,CTLA4_start+5e3)
region3=c(ICOS_start-5e3,ICOS_start+5e3)
plot_grid(
  merged %>% dplyr::filter(!(target=='CTLA4'&condition=='0h')&between(start,region1[1],region1[2])) %>%
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
    coord_cartesian(xlim = region1, expand=FALSE, clip = "off")+
    theme_classic()+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white',fill='white'),
          axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  merged %>% dplyr::filter(!(target=='CTLA4'&condition=='0h')&between(start,region2[1],region2[2])) %>%
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
    coord_cartesian(xlim = region2, expand=FALSE, clip = "off")+
    theme_classic()+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white',fill='white'),
          axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  merged %>% dplyr::filter(!(target=='CTLA4'&condition=='0h')&between(start,region3[1],region3[2])) %>%
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
    coord_cartesian(xlim = region3, expand=FALSE, clip = "off")+
    theme_classic()+
    theme(strip.placement = 'outside', strip.background = element_rect(color='white',fill='white'),
          axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.spacing = unit(0.7, "lines"))+
    NULL,
  nrow=1,align='h',axis='tb', rel_widths = c(1,1,1))

## Arrayed Validation
# CRISPRi
plot = read_rds('220319-1_merged_unprocessed.rds') %>% 
  dplyr::filter(count>500) %>% 
  group_by(celltype,donor,flowTarget) %>% 
  dplyr::mutate(ctrl_median = median(MFI[grepl('NTC',sgrnaTarget)], na.rm=T)) %>%
  dplyr::mutate(fc = MFI/ctrl_median)
plot %>% dplyr::filter(sgrnaTarget!='NTC') %>% 
  ggplot(aes(sgrnaTarget,fc,fill=flowTarget))+
  geom_boxplot(color='black',outlier.shape = NA)+
  geom_point(color='black',alpha=0.5,shape=16,size=rel(2),position=position_dodge(width=0.75))+
  facet_grid(rows=vars(celltype))+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  scale_fill_manual(values=c('black','grey60','grey90'))+
  scale_y_continuous(breaks = c(0.0,0.5,1.0,1.5),limits = c(0,1.8))
plot %>% group_by(celltype, sgrnaTarget,flowTarget,donor) %>% dplyr::count() %>% view()
compare_means(data=plot,
              formula = fc ~ sgrnaTarget,
              method='t.test',
              group.by = c('celltype','flowTarget'),
              comparisons = list(c('NTC','CTLA4 TSS')))

# KO
ko=read_rds('220105-1_merged_unprocessed.rds') %>% 
  dplyr::filter(count>500) %>% 
  group_by(celltype,donor,flowTarget,condition) %>% 
  dplyr::mutate(ctrl_median = median(MFI[grepl('AAVS1',sgrnaTarget)], na.rm=T)) %>%
  dplyr::mutate(fc = MFI/ctrl_median)
ko %>% dplyr::filter(sgrnaTarget!='AAVS1') %>% 
  ggplot(aes(sgrnaTarget,fc,fill=flowTarget))+
  geom_boxplot(color='black',outlier.shape = NA)+
  geom_point(color='black',alpha=0.5,shape=16,size=rel(2),position=position_dodge(width=0.75))+
  facet_grid(rows=vars(celltype))+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  scale_fill_manual(values=c('black','grey60','grey90'))+
  scale_y_continuous(breaks = c(0.0,0.5,1.0,1.5),limits = c(0,1.8))
ko %>% group_by(celltype, sgrnaTarget,flowTarget,donor) %>% dplyr::count() %>% view()
compare_means(data=ko %>% dplyr::filter(sgrnaTarget %in% c('AAVS1','CTLA4')),
              formula = fc ~ sgrnaTarget,
              method='t.test',
              group.by = c('celltype','flowTarget'),
              comparisons = list(c('AAVS1','CTLA4')))


#####
## Figure S7
#####
## Promoter-Capture-C from Su et al., Nat. Commun. 2020
##  Processing
ibedToBedpe = function(file){
  converted = read_tsv(paste0(file,'.ibed')) %>%
    dplyr::mutate(name=paste0(bait_name,'___',otherEnd_name)) %>%
    dplyr::select(chrom1=bait_chr,start1=bait_start,end1=bait_end,
                  chrom2=otherEnd_chr,start2=otherEnd_start,end2=otherEnd_end,
                  name,score) %>%
    dplyr::mutate(start1=start1-1,start2=start2-1)
  write_tsv(converted, file = paste0(file,'.bedpe'), col_names = F)
}
ibedToBedpe_targetSubset = function(file){
  converted = read_tsv(paste0(file,'.ibed')) %>%
    dplyr::mutate(name=paste0(bait_name,'___',otherEnd_name)) %>%
    dplyr::select(chrom1=bait_chr,start1=bait_start,end1=bait_end,
                  chrom2=otherEnd_chr,start2=otherEnd_start,end2=otherEnd_end,
                  name,score) %>%
    dplyr::mutate(start1=start1-1,start2=start2-1)
  for(target in c('CD28','CTLA4','ICOS')){
    sub = converted %>% dplyr::filter(grepl(target,name))
    write_tsv(sub, file = paste0(file,"_",target,'.bedpe'), col_names = F)
  }
}

ibeds=list.files(pattern = 'ibed')
files=gsub("\\.ibed","", ibeds)

###### Loading
file_type = '1frag_'
files=list.files(pattern = file_type)
combined = files %>%
  map(function(x) read_tsv(x, col_names=F) %>% dplyr::mutate(file=x)) %>%
  do.call("rbind",.) %>%
  dplyr::rename(chrom1=X1,start1=X2,end1=X3,chrom2=X4,start2=X5,end2=X6,name=X7,score=X8) %>%
  tidyr::separate(file, into = c('celltype','data','type'),sep = '\\.') %>% tidyr::separate(data, into = c('frag','target'),sep = '_') %>%
  rowwise() %>% dplyr::mutate(midstart=min(mean(start1,end1),mean(start2,end2)),midend=max(mean(start1,end1),mean(start2,end2)))
combined

## Liftover
lifted = as_tibble(unlist(
  makeGRangesFromDataFrame(combined %>% dplyr::rename(seqnames=chrom1,start=midstart,end=midend),keep.extra.columns = T) %>%
    liftOver(., chain=import.chain('~/Dropbox (Gladstone)/Data/TADscreen/reference_data/hg19ToHg38.over.chain'))))
write_rds(lifted, paste0(file_type,'merged_hg38Lifted.rds'))

pcc = read_rds('1frag_merged_hg38Lifted.rds')
RAPH1_start=203433682 ## RAPH1 NM_213589.3 dominant transcript, hg38
RAPH1_stop=203535301 ## RAPH1 NM_213589.3 dominant transcript, hg38
plot_region = c(TAD_stop-500e3,TAD_stop+100e3) # RAPH1 to ICOS zoom
plot_grid(
  ggplot(pcc)+
    annotate("rect", xmin=RAPH1_start, xmax=RAPH1_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
    geom_curve(aes(x = start, y = 0, xend = end, yend = 0, color=celltype), curvature = -0.35)+
    coord_cartesian(xlim = plot_region,ylim=c(0,2), expand=FALSE)+
    scale_color_manual(values=c('grey60','black'))+
    facet_grid(rows=vars(target))+
    theme_void()+
    NULL,
  merged %>% dplyr::filter(between(start,plot_region[1],plot_region[2])) %>%
    dplyr::mutate(pval_ceiling = ifelse(padj<0.05,log10(0.05),log10(padj))) %>%
    dplyr::mutate(pval_ceiling_dir = ifelse(log2FoldChange>0,pval_ceiling,-pval_ceiling)) %>%
    dplyr::mutate(celltype=factor(celltype,levels=c('Tconv','Treg'))) %>% 
    ggplot()+
    annotate("rect", xmin=RAPH1_start, xmax=RAPH1_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.5)+
    geom_hline(yintercept = 0, linetype='solid', color='black')+
    geom_vline(xintercept = TAD_stop, linetype='solid', color='black')+
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
  ncol = 1, align = 'v', axis = 'lr', rel_heights = c(1,1))


##  Co-regulation of Trans factors
data = trans %>% 
  dplyr::select(screen,lfc,id) %>% 
  distinct() %>% 
  pivot_wider(names_from = screen, values_from = lfc) %>% 
  dplyr::mutate(from0 = abs((0-CD28)+(0-CTLA4)+(0-ICOS))) %>% 
  arrange(desc(from0))
## Limit to just TFs significiant in any 1 gene
data2 = trans %>% 
  dplyr::filter(id %in% (trans %>% dplyr::filter(fdr<0.05) %>% pull(id))) %>% 
  dplyr::select(screen,lfc,id) %>% 
  distinct() %>% 
  pivot_wider(names_from = screen, values_from = lfc) %>% 
  dplyr::mutate(from0 = abs((0-CD28)+(0-CTLA4)+(0-ICOS))) %>% 
  arrange(desc(from0))
plot_grid(
  ggplot(data=data2, 
         aes(CD28, CTLA4))+
    geom_density_2d(data=data, aes(CD28, CTLA4), color='grey70')+
    geom_point()+
    geom_smooth(method = 'glm', color='black'),
  ggplot(data=data2, aes(CD28, ICOS))+
    geom_density_2d(data=data, aes(CD28, ICOS), color='grey70')+
    geom_point()+
    geom_smooth(method = 'glm', color='black'),
  ggplot(data=data2, 
         aes(CTLA4,ICOS))+
    geom_density_2d(data=data, aes(CTLA4,ICOS), color='grey70')+
    geom_point()+
    geom_smooth(method = 'glm', color='black'),
  nrow=1,align='h',axis='tb', rel_widths = c(1,1,1))
cor.test(data2$CD28,data2$CTLA4,
         method = c("pearson"))
cor.test(data2$CD28,data2$ICOS,
         method = c("pearson"))
cor.test(data2$ICOS,data2$CTLA4,
         method = c("pearson"))
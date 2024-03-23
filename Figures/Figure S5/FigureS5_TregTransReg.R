## Cody Mowery
source('sourceFigures.R')

## Arrayed Trans KO in Human Tregs
plot=read_rds('230313-1_merged_unprocessed.rds') %>%
  dplyr::filter(Count>500) %>%
  group_by(Donor,flowTarget) %>%
  dplyr::mutate(ctrl_median = median(MFI[grepl('AAVS1',Gene)], na.rm=T)) %>%
  dplyr::mutate(fc = MFI/ctrl_median)

plot %>% ggplot(aes(Gene,fc))+
  geom_hline(yintercept = 1)+
  geom_boxplot(color='black',outlier.shape = NA)+
  geom_point(color='black', shape=16, alpha=0.5)+
  facet_grid(cols=vars(flowTarget))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+ 
  NULL
table(plot$Donor, plot$Gene)
compare_means(data=plot,
              formula = fc ~ Gene,
              ref.group = 'AAVS1',
              group.by = c('flowTarget'),
              method='t.test')

## Motifs
motifs = read_tsv('JASPAR22_motifs_chr2_202527032_203967032.txt', col_names = F) %>% 
  dplyr::filter(X5>400) %>% ## Default JASPAR 22 filter in UCSC
  dplyr::filter(grepl('FOXP3$|STAT5',X7,ignore.case = T))
plot_region = c(203833771,203837771) ## Treg-Dom Enhancer
ggplot(motifs %>% dplyr::filter(between(X2,plot_region[1],plot_region[2]) & !grepl('::',X7)))+
  annotate("rect", xmin=CD28_start, xmax=CD28_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
  annotate("rect", xmin=CTLA4_start, xmax=CTLA4_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
  annotate("rect", xmin=ICOS_start, xmax=ICOS_stop, ymin=-Inf, ymax=Inf, color=NA, fill='grey90', alpha=0.25)+
  geom_rect(aes(xmin=X2, xmax=X3, ymin=0, ymax=1), color='black')+
  coord_cartesian(xlim = plot_region, expand=FALSE)+
  facet_grid(rows=vars(X7))+
  theme_void()+
  NULL

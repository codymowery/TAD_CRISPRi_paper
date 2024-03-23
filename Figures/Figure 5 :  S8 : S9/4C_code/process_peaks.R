###########################################################
## Extract the normalized 4C values and represent the average
## value in a boxplot 


library(peakC)

## Load datasets

setwd('./pipe4Cout/RDS/')

rds1 <- readRDS('CTLA4StimEn_HandD_001.rds')
rds2 <- readRDS('CTLA4StimEn_HandD_002.rds')
rds3 <- readRDS('CTLA4StimEn_HandD_003.rds')
rds4 <- readRDS('CTLA4StimEn_HandD_004.rds')


## Plot Genes

library(scales)

##CTLA4

pdf(file = 'boxplot_CTLA4_donor1_blue.pdf',width = 4,height = 7)
boxplot(log(subsetByOverlaps(rds2$reads, GRanges('chr2',IRanges( 203867771, 203873965)))$norm4C), 
        log(subsetByOverlaps(rds1$reads, GRanges('chr2',IRanges( 203867771, 203873965)))$norm4C),
        frame.plot=F, main = 'CTLA4', ylab='log(norm4C)', names=c('HandD_002','HandD_001'),
        col=c(alpha('blue',0.75),alpha('red3',0.75)))
dev.off()

t.test(log(subsetByOverlaps(rds2$reads, GRanges('chr2',IRanges( 203867771, 203873965)))$norm4C), 
       log(subsetByOverlaps(rds1$reads, GRanges('chr2',IRanges( 203867771, 203873965)))$norm4C))$p.value

pdf(file = 'boxplot_CTLA4_donor2_blue.pdf',width = 4,height = 7)
boxplot(log(subsetByOverlaps(rds4$reads, GRanges('chr2',IRanges( 203867771, 203873965)))$norm4C), 
        log(subsetByOverlaps(rds3$reads, GRanges('chr2',IRanges( 203867771, 203873965)))$norm4C),
        frame.plot=F, main = 'CTLA4', ylab='log(norm4C)', names=c('HandD_004','HandD_003'),
        col=c(alpha('blue',0.75),alpha('red3',0.75)))
dev.off()

t.test(log(subsetByOverlaps(rds4$reads, GRanges('chr2',IRanges( 203867771, 203873965)))$norm4C), 
       log(subsetByOverlaps(rds3$reads, GRanges('chr2',IRanges( 203867771, 203873965)))$norm4C))$p.value



##CD28

pdf(file = 'boxplot_CD28_blue_donor1.pdf',width = 4,height = 7)
boxplot(log(subsetByOverlaps(rds2$reads, GRanges('chr2',IRanges( 203706475, 203738912)))$norm4C),
        log(subsetByOverlaps(rds1$reads, GRanges('chr2',IRanges( 203706475, 203738912)))$norm4C),
        frame.plot=F, main = 'CD28', ylab='log(norm4C)', names=c('HandD_002','HandD_001'),
        col=c(alpha('blue',0.75),alpha('red3',0.75)))
dev.off()

t.test(log(subsetByOverlaps(rds2$reads, GRanges('chr2',IRanges( 203706475, 203738912)))$norm4C),
       log(subsetByOverlaps(rds1$reads, GRanges('chr2',IRanges( 203706475, 203738912)))$norm4C))$p.value

pdf(file = 'boxplot_CD28_blue_donor2.pdf',width = 4,height = 7)
boxplot(log(subsetByOverlaps(rds4$reads, GRanges('chr2',IRanges( 203706475, 203738912)))$norm4C), 
        log(subsetByOverlaps(rds3$reads, GRanges('chr2',IRanges( 203706475, 203738912)))$norm4C),
        frame.plot=F, main = 'CD28', ylab='log(norm4C)', names=c('HandD_004','HandD_003'),
        col=c(alpha('blue',0.75),alpha('red3',0.75)))
dev.off()

t.test(log(subsetByOverlaps(rds4$reads, GRanges('chr2',IRanges( 203706475, 203738912)))$norm4C), 
       log(subsetByOverlaps(rds3$reads, GRanges('chr2',IRanges( 203706475, 203738912)))$norm4C))$p.value


##ICOS

pdf(file = 'boxplot_ICOS_blue_donor1.pdf',width = 4,height = 7)
boxplot(log(subsetByOverlaps(rds2$reads, GRanges('chr2',IRanges( 203936763, 203961577)))$norm4C),
        log(subsetByOverlaps(rds1$reads, GRanges('chr2',IRanges( 203936763, 203961577)))$norm4C),
        frame.plot=F, main = 'ICOS', ylab='log(norm4C)', names=c('HandD_002','HandD_001'),
        col=c(alpha('blue',0.75),alpha('red3',0.75)))
dev.off()

t.test(log(subsetByOverlaps(rds2$reads, GRanges('chr2',IRanges( 203936763, 203961577)))$norm4C),
       log(subsetByOverlaps(rds1$reads, GRanges('chr2',IRanges( 203936763, 203961577)))$norm4C))$p.value

pdf(file = 'boxplot_ICOS_blue_donor2.pdf',width = 4,height = 7)
boxplot(log(subsetByOverlaps(rds4$reads, GRanges('chr2',IRanges( 203936763, 203961577)))$norm4C), 
        log(subsetByOverlaps(rds3$reads, GRanges('chr2',IRanges( 203936763, 203961577)))$norm4C),
        frame.plot=F, main = 'ICOS', ylab='log(norm4C)', names=c('HandD_004','HandD_003'),
        col=c(alpha('blue',0.75),alpha('red3',0.75)))
dev.off()

t.test(log(subsetByOverlaps(rds4$reads, GRanges('chr2',IRanges( 203936763, 203961577)))$norm4C), 
       log(subsetByOverlaps(rds3$reads, GRanges('chr2',IRanges( 203936763, 203961577)))$norm4C))$p.value


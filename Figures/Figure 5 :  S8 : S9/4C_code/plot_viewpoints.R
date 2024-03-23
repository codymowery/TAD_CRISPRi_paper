############################################################
## Create coverage plots using the wig files and the beds

library('rtracklayer')
library('Sushi')

## Gene tracks

gene.file <- read.table(gzfile('~/Downloads/hg38.refGene.gtf.gz'), sep = '\t', header = F)
gene.file <- gene.file[grep('exon',gene.file$V3),]
gene.name <- unlist(lapply(strsplit(gene.file$V9,';'), function(x){gsub(" ","",gsub(" gene_name","",x[grep('gene_name',x)]))}))
gene.bed <- cbind(gene.file[,c(1,4,5)], gene.name, gene.file[,c(6,7)])
colnames(gene.bed) <- c('chrom','start','stop','gene','score','strand')
gene.bed <- gene.bed[-c(grep(pattern = 'alt', gene.bed$chrom),grep(pattern = 'fix', gene.bed$chrom),grep(pattern = 'random', gene.bed$chrom)),]


## Load wig files and process to include chromosome

CTLA4StimEn_HandD_001_WIN21.wig <- read.table('pipe4Cout/WIG/CTLA4StimEn_HandD_001_WIN21.wig.gz', sep = '\t', header = F, skip = 3)
CTLA4StimEn_HandD_001_WIN21.bed <- read.table('pipe4Cout/CTLA4StimEn_HandD_001.bed', sep = '\t', header = F, skip = 2)
CTLA4StimEn_HandD_001_WIN21.wig <- cbind('chr2',CTLA4StimEn_HandD_001_WIN21.wig)
new.col <- c()
for (n in CTLA4StimEn_HandD_001_WIN21.wig[,2]){
  new.col <- c(new.col,n+100)
}
#new.col <- c(new.col[2:length(new.col)],new.col[length(new.col)]+1)
CTLA4StimEn_HandD_001_WIN21.wig <- cbind(CTLA4StimEn_HandD_001_WIN21.wig[,c(1,2)],new.col,CTLA4StimEn_HandD_001_WIN21.wig[,3])
rownames(CTLA4StimEn_HandD_001_WIN21.wig) <- CTLA4StimEn_HandD_001_WIN21.wig[,2]

CTLA4StimEn_HandD_002_WIN21.wig <- read.table('pipe4Cout/WIG/CTLA4StimEn_HandD_002_WIN21.wig.gz', sep = '\t', header = F, skip = 3)
CTLA4StimEn_HandD_002_WIN21.bed <- read.table('pipe4Cout/CTLA4StimEn_HandD_002.bed', sep = '\t', header = F, skip = 2)
CTLA4StimEn_HandD_002_WIN21.wig <- cbind('chr2',CTLA4StimEn_HandD_002_WIN21.wig)
new.col <- c()
for (n in CTLA4StimEn_HandD_002_WIN21.wig[,2]){
  new.col <- c(new.col,n+100)
}
#new.col <- c(new.col[2:length(new.col)],new.col[length(new.col)]+1)
CTLA4StimEn_HandD_002_WIN21.wig <- cbind(CTLA4StimEn_HandD_002_WIN21.wig[,c(1,2)],new.col,CTLA4StimEn_HandD_002_WIN21.wig[,3])
rownames(CTLA4StimEn_HandD_002_WIN21.wig) <- CTLA4StimEn_HandD_002_WIN21.wig[,2]

CTLA4StimEn_HandD_003_WIN21.wig <- read.table('pipe4Cout/WIG/CTLA4StimEn_HandD_003_WIN21.wig.gz', sep = '\t', header = F, skip = 3)
CTLA4StimEn_HandD_003_WIN21.bed <- read.table('pipe4Cout/CTLA4StimEn_HandD_003.bed', sep = '\t', header = F, skip = 2)
CTLA4StimEn_HandD_003_WIN21.wig <- cbind('chr2',CTLA4StimEn_HandD_003_WIN21.wig)
new.col <- c()
for (n in CTLA4StimEn_HandD_003_WIN21.wig[,2]){
  new.col <- c(new.col,n+100)
}
#new.col <- c(new.col[2:length(new.col)],new.col[length(new.col)]+1)
CTLA4StimEn_HandD_003_WIN21.wig <- cbind(CTLA4StimEn_HandD_003_WIN21.wig[,c(1,2)],new.col,CTLA4StimEn_HandD_003_WIN21.wig[,3])
rownames(CTLA4StimEn_HandD_003_WIN21.wig) <- CTLA4StimEn_HandD_003_WIN21.wig[,2]

CTLA4StimEn_HandD_004_WIN21.wig <- read.table('pipe4Cout/WIG/CTLA4StimEn_HandD_004_WIN21.wig.gz', sep = '\t', header = F, skip = 3)
CTLA4StimEn_HandD_004_WIN21.bed <- read.table('pipe4Cout/CTLA4StimEn_HandD_004.bed', sep = '\t', header = F, skip = 2)
CTLA4StimEn_HandD_004_WIN21.wig <- cbind('chr2',CTLA4StimEn_HandD_004_WIN21.wig)
new.col <- c()
for (n in CTLA4StimEn_HandD_004_WIN21.wig[,2]){
  new.col <- c(new.col,n+100)
}
#new.col <- c(new.col[2:length(new.col)],new.col[length(new.col)]+1)
CTLA4StimEn_HandD_004_WIN21.wig <- cbind(CTLA4StimEn_HandD_004_WIN21.wig[,c(1,2)],new.col,CTLA4StimEn_HandD_004_WIN21.wig[,3])
rownames(CTLA4StimEn_HandD_004_WIN21.wig) <- CTLA4StimEn_HandD_004_WIN21.wig[,2]


## Calculate Splines for smoothing and Plot region using overlaps between wt and ko

library(scales)


pdf(file = 'plot4C_region_sushi_continuous.pdf',width = 4,height = 6)
par(mfrow=c(3,1),mar = c(0.7, 6.5, 1.4, 6.5))
pg = plotGenes(geneinfo = gene.bed[which(gene.bed$chrom=='chr2'),], chrom = 'chr2', chromstart = 203494157 ,chromend = 203972174 ,
               bheight=0.2,plotgenetype="arrow",bentline=FALSE,
               labeloffset=.4,fontsize=0.9,arrowlength = 0.025,
               labeltext=TRUE, labelat = 'middle')
labelgenome('chr2',203494157,203972174,n=3,scale="Mb")

par(mar = c(0.7, 5.5, 1.4, 5.5))

CTLA4StimEn_HandD_001_WIN21_red.wig <- CTLA4StimEn_HandD_001_WIN21.wig[CTLA4StimEn_HandD_001_WIN21.wig$V1>203494150&CTLA4StimEn_HandD_001_WIN21.wig$V1<203972181,]
spline.data = smooth.spline(x=CTLA4StimEn_HandD_001_WIN21_red.wig[,1], y=CTLA4StimEn_HandD_001_WIN21_red.wig[,2], spar = 0.75, all.knots=T)
plot(spline.data, type='l', ylim=c(0,3000), col = 'red3', frame.plot=F, xaxt='n')
xy=predict(spline.data, CTLA4StimEn_HandD_001_WIN21.wig[CTLA4StimEn_HandD_001_WIN21.wig$V1>203494150&CTLA4StimEn_HandD_001_WIN21.wig$V1<203972181,1])
m <- length(xy$x)
x.poly <- c(xy$x, xy$x[m], xy$x[1])
y.poly <- c(xy$y, 0, 0)
polygon(x.poly, y.poly, col = alpha('red3', 0.6), border = NA)

CTLA4StimEn_HandD_002_WIN21_red.wig <- CTLA4StimEn_HandD_002_WIN21.wig[CTLA4StimEn_HandD_002_WIN21.wig$V1>203494150&CTLA4StimEn_HandD_002_WIN21.wig$V1<203972181,]
spline.data = smooth.spline(x=CTLA4StimEn_HandD_002_WIN21_red.wig[,1], y=CTLA4StimEn_HandD_002_WIN21_red.wig[,2], spar = 0.75, all.knots=T)
lines(spline.data, type='l', ylim=c(0,3000), col = 'darkblue')
xy=predict(spline.data, CTLA4StimEn_HandD_002_WIN21.wig[CTLA4StimEn_HandD_002_WIN21.wig$V1>203494150&CTLA4StimEn_HandD_002_WIN21.wig$V1<203972181,1])
m <- length(xy$x)
x.poly <- c(xy$x, xy$x[m], xy$x[1])
y.poly <- c(xy$y, 0, 0)
polygon(x.poly, y.poly, col = alpha('blue', 0.3), border = NA)

CTLA4StimEn_HandD_003_WIN21_red.wig <- CTLA4StimEn_HandD_003_WIN21.wig[CTLA4StimEn_HandD_003_WIN21.wig$V1>203494150&CTLA4StimEn_HandD_003_WIN21.wig$V1<203972181,]
spline.data = smooth.spline(x=CTLA4StimEn_HandD_003_WIN21_red.wig[,1], y=CTLA4StimEn_HandD_003_WIN21_red.wig[,2], spar = 0.75, all.knots=T)
plot(spline.data, type='l', ylim=c(0,3000), col = 'red3', frame.plot=F, xaxt='n')
xy=predict(spline.data, CTLA4StimEn_HandD_003_WIN21.wig[CTLA4StimEn_HandD_003_WIN21.wig$V1>203494150&CTLA4StimEn_HandD_003_WIN21.wig$V1<203972181,1])
m <- length(xy$x)
x.poly <- c(xy$x, xy$x[m], xy$x[1])         
y.poly <- c(xy$y, 0, 0)                     
polygon(x.poly, y.poly, col = alpha('red3', 0.6), border = NA)

CTLA4StimEn_HandD_004_WIN21_red.wig <- CTLA4StimEn_HandD_004_WIN21.wig[CTLA4StimEn_HandD_004_WIN21.wig$V1>203494150&CTLA4StimEn_HandD_004_WIN21.wig$V1<203972181,]
spline.data = smooth.spline(x=CTLA4StimEn_HandD_004_WIN21_red.wig[,1], y=CTLA4StimEn_HandD_004_WIN21_red.wig[,2], spar = 0.75, all.knots=T)
lines(spline.data, type='l', ylim=c(0,3000), col = 'darkblue')
xy=predict(spline.data, CTLA4StimEn_HandD_004_WIN21.wig[CTLA4StimEn_HandD_004_WIN21.wig$V1>203494150&CTLA4StimEn_HandD_004_WIN21.wig$V1<203972181,1])
m <- length(xy$x)
x.poly <- c(xy$x, xy$x[m], xy$x[1])
y.poly <- c(xy$y, 0, 0)
polygon(x.poly, y.poly, col = alpha('blue', 0.3), border = NA)

dev.off()


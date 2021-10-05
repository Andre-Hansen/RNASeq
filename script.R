read.csv("~/RNAseq1/SraRunTable.txt")
runtable<-read.csv("~/RNAseq1/SraRunTable.txt")
library("ggpubr")
library("stringr")


runtable<-data.frame(runtable)
runtab<-data.frame(runtable$Run)
runtab$Lib<-runtable$Library.Name
runtab$ID<-runtable$submitted_subject_id
runtab$time<-NA
runtab$time<-rep("post", 20)
runtab$time<-rep("post", 19)
runtab[21:nrow(runtab),]$time<-rep("pre", 18)
runtab$filename<-NA
runtab$filename<-paste0(runtab$runtable.Run, "*.bam")
runtab$IDtime<-NA
runtab$IDtime<-paste0(runtab$ID, "_", runtab$time)
#filenames<-list.files("RNAseq1/",pattern="*bam")
filenames<-paste0(runtab$runtable.Run, ".bam")
library("Rsamtools")
bamfiles<-BamFileList(filenames, yieldSize = 100000)
library("GenomicFeatures")
gfffile<-("/Volumes/scratch/cancergeneticslab/ConorM/ref_genomes/human_GRCh38/Homo_sapiens.GRCh38.97.chr.gff3")
txdb<-makeTxDbFromGFF(gfffile, format = "gff3", organism = "Homo sapiens")
saveDb(txdb, file="HomoSapiens.sqlite")
txdb<-loadDb("HomoSapiens.sqlite")
ebg <- exonsBy(txdb, by="gene")
ebg
library("GenomicAlignments")
library("BiocParallel")
setwd("RNAseq1/")
se<-summarizeOverlaps(ebg, bamfiles, singleEnd=FALSE, ignore.strand=TRUE)
eMatrix<-assays(se)$counts
library(utils)
library(DESeq2)
conds=data.frame(conds=factor(runtab$time), row.names = rownames(runtab))
#dds <- DESeqDataSetFromMatrix(countData = as.matrix(eMatrix),colData = conds,design = formula(~conds))

##RPKM values
returnRPKM <- function(counts, gffsub) {
  # Length of exon union per gene in kbp
  geneLengthsInKB <- sum(width(reduce(gffsub)))/1000
  # Factor for converting to million of mapped reads. 
  millionsMapped <- sum(counts)/1e+06
  # RPK: reads per kilobase of exon model.
  rpm <- counts/millionsMapped
  # RPKM: reads per kilobase of exon model per million mapped reads. 
  rpkm <- rpm/geneLengthsInKB
  return(rpkm)
}

#countRrpkm <- apply(eMatrix, 2, function(x) returnRPKM(counts=x, gffsub=ebg))
#library(tidyverse)
#counts<-data.frame(countRrpkm, row.names = FALSE)
#counts<-rownames_to_column(counts, quote = FALSE)

?rownames_to_column
colnames(counts)[1]<-"Gene symbol"
write.table(counts,"~/millerfiles/RNAseq/RPKM.txt", sep = "\t", row.names = FALSE, quote = FALSE)
rd<-read.delim("~/millerfiles/RNAseq/RPKM.txt")
read.delim("~/faimos data/SymbolsOnly.txt")
##deseq differential between pre/post
Deseq<-DESeq(dds)
res<-results(Deseq)
res<-na.omit(res)
sum(res$padj <= 0.05)
sum(p.adjust(res$pvalue, nrow(res), method = "fdr") <= 0.05)
sigres<-res[p.adjust(res$pvalue, nrow(res), method = "fdr") <= 0.05,]
library(lattice)
library(gplots)
rownames(sigres)[1:20]
##limma
library(limma)
library(edgeR)
conds$conds
design<-model.matrix(~conds$conds)
dge<-DGEList(counts=eMatrix)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge<-calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
v <- voom(dge, design, plot=TRUE, normalize.method = "none")
for(i in 1:3){
  plot(logCPM[,i], v$E[,i], xlab="LogCPM", ylab="Voom", main=colnames(logCPM)[i])
  abline(0,1) }
fit <- lmFit(v, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef=ncol(design), n=nrow(eMatrix)) 
sum(p.adjust(tt$P.Value, n= nrow(eMatrix), method = "fdr")<=0.05)
sig<-tt[p.adjust(tt$P.Value, n= nrow(eMatrix), method = "fdr")<=0.05,]

##gene list conversion
#1. get your list of genes as a vector
#2. convert to ENSG symbols using bioMart
#3.extract these genes from v$E or normalised expression mtrix, and create a heatmap


library(lattice)
library(gplots)
rownames(sig)

y<-norm[rownames(sig),]
y<-t(scale(t(as.matrix(y))))
y <- y[order(y[,1]),]

levelplot(t(y), height=0.2, col.regions=colorpanel(40, "blue", "white", "red"), colorkey=list(space="top"), xlab="", ylab="Gene ID")

#####CIBERSORT
csO<-read.csv("Miller RNAseq CIBERSORT.csv")
csO$Timepoint<-NA
csO$Timepoint<-runtab$time
csO$ID<-runtab$ID
prepost<-csO
norm<-logCPM
genelist<-c("PRF1", "ARHGAP15", "GZMB", "ARHGAP25", "CXCL13", "CCL5",
            'IRF1', "CCR2", "IKZF1", "CCR7", "HLA-E", "CD2",
            "CD200","CXCL13", "FBLN7", "ICOS", "SGPP2", "SH2D1A",
            "FOXP3", "GZMA", "C15orf53", "PRF1", "IL5",
            "CTLA4",
            "IL32", "GPR15","IL4",
            "HLA-DMA", "HLA-DQB1", "HLA-DRA","HLA-DRB4",
            "KLRF1","KLRC1",
            "FUCA1","MMP9","LGMN","HS3ST2", "TM4SF19", "CLEC5A","GPNMB", "C11orf45", "CD68", "CYBB",
            "HLA-A", "HLA-B", "HLA-C","HLA-F", "HLA-G", "HLA-J",
            "CXCL10", "CXCL11", "GBP1", "STAT1",
            "DDX58", "CD2", "HERC6", "CD226", "IFI44", "CD27", "IFI44L", "CD28", "IFIT1", "CD40", "IFIT2", "CD40LG", "MX1", "CD58","OAS1", "CD70", "OAS3", "ICOS", "RSAD2", "ICOSLG",
            "BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","LGALS9",
            "CD247","TIGIT","CD27","PDCD1","CD3D",
            "CD48",
            "CD53","CORO1A","CSF2RB","EVI2B","FGL2","GIMAP4","GIMAP5","GMFG","GZMA","GZMK","HCLS1","IL10RA","IL2RG","IL7R","INPP5D","IRF8",
            "ITK","KLRK1","LCK","LCP2","LPXN",
            "LTB","PIK3CD","PLAC8","PRG1","PRKCB1","PTPRC","RAC2","SAMSSN1","SCYA5","SELL","SD2D1A","SLA","SLAMF1","STAT4","TNFRSF7","TRBC1",
            "SLAMF1","PDCD1LG2","TNFRSF18","PVRL3","TNFRSF25","TIGIT","TNFRSF4",
            "TNFRSF8","TNFRSF9","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9")
library(GenomicRanges)
#normRPKM <- apply(norm, 2, function(x) returnRPKM(counts=x, gffsub=ebg))
#norm<-data.frame(norm)
#norm
#norm[which(rownames(norm) =="RAB27B"),]

GS<-NULL
for(i in genelist){
  GS<-rbind(norm[which(rownames(norm)==i),], GS)
  print(i)
}

GS

heatmap.2(data.matrix(GS), trace="none")
immgenes<-hclust(dist(t(GS)))

hiclust<-cutree(immgenes, 3)
plot(as.dendrogram(immgenes), leaflab = "none")
prepost$cluster<- hiclust[match(prepost$Input.Sample,names(hiclust))]
immpops<-colnames(prepost[c(2:23, 29)])
ggpaired(prepostwide, cond1 = "T.cells.CD8.pre", cond2 = "T.cells.CD8.post", width = 0)
table(prepost$cluster)
for (i in immpops){
  g<-ggplot(prepost, aes(as.factor(cluster), eval(parse(text=paste0(i))))) + 
    ylab(i)+
    geom_boxplot(alpha = 0.5) +
    geom_point(position = "jitter")+
    geom_signif(comparisons = list(c(1,2), c(1,3), c(2,3)),
                y_position = c(max(eval(parse(text=paste0("prepost$",i)))),
                               max(eval(parse(text=paste0("prepost$",i))))+0.1*max(eval(parse(text=paste0("prepost$",i)))),
                               max(eval(parse(text=paste0("prepost$",i))))+0.2*max(eval(parse(text=paste0("prepost$",i)))),
                               max(eval(parse(text=paste0("prepost$",i))))+0.3*max(eval(parse(text=paste0("prepost$",i)))), map_signif_level = TRUE))
  print(g)
}
library(forcats)
as.factor(prepost$Timepoint)
prepost$Timepoint<-fct_relevel(prepost$Timepoint, c("pre", "post"))
for (i in immpops){
  g<-ggplot(prepost,aes(Timepoint, eval(parse(text=paste0(i))))) + geom_boxplot(alpha = 0.5) + 
    geom_point(position = "jitter") +
    ylab(i) + 
    geom_signif(comparisons = list(c("pre", "post")), map_signif_level = TRUE,
                y_position = max(eval(parse(text=paste0("prepost$",i))))+0.1*max(eval(parse(text=paste0("prepost$",i)))))
  print(g)
}
t.test(prepost$T.cells.CD8[which(prepost$Timepoint == "pre")])
prepostwide<-reshape(prepost, 
                     v.names= immpops, 
                     timevar= "Timepoint", 
                     idvar = "ID", 
                     direction = "wide")
immpopchanges<-NULL
for (i in immpops){
  pre<-eval(parse(text=paste0("prepostwide$", i, ".pre")))
  post<-eval(parse(text=paste0("prepostwide$", i, ".post")))
  doop<-post-pre
  immpopchanges<-cbind(doop, immpopchanges)
}
rev(immpops)
colnames(immpopchanges)<-rev(immpops)
immpopchanges<-data.frame(immpopchanges)
immpopchanges$ID<-prepostwide$ID
immpopchanges<-immpopchanges[!is.na(immpopchanges$Neutrophils),]
pvals<-list()
for (i in immpops){
  o<-order(eval(parse(text=paste0("immpopchanges$",i))))
  immpopchanges<-immpopchanges[o,]
  #use a cutoff of 50 as used in Dunbier 2013 paper
  #col<- ifelse(immpopchangeski67$ki67 > 50, "red", "green")
  barplot(eval(parse(text=paste0("immpopchanges$",i)))[!is.na(eval(parse(text=paste0("immpopchanges$",i))))],
          main = paste("Waterfall plot for changes in", i, "(post-pre)"),
          xlab= "Patient",
          ylab=paste(i, "gene signature shift (post-pre treatment)"))
  #plot(eval(parse(text=paste0("immpopchangeski67$",i))), immpopchangeski67$ki67,
  #xlab=paste(i, "gene signature shift (post-pre treatment)"),
  #lab="Ki67 change")
  g0<-sum(eval(parse(text=paste0("immpopchanges$",i)))>0, na.rm = TRUE)
  l0<-sum(eval(parse(text=paste0("immpopchanges$",i)))<0, na.rm = TRUE)
  ##binomial test-do the values go up or down?
  binom<-binom.test(c(g0,l0))
  pvals[i]<-binom$p.value
}
pvals
names(counts)
sum(tt$)
head(tt)


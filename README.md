# ChIPATAC

Tools for analyzing genomic regions. 

### readBed()
Read in a bed file as a GenomicRanges object:
```
x<-readBed("path/to/bedfile.bed")
```

### writeGTF()
Export a GenomicRanges object as a GTF file:
```
writeGTF(x, source="ATAC-Peakset", peakID="peak", filename="path/to/gtffile.gtf")
```

### writeBETA()
Export a GenomicRanges object as a bed file in the format compatible with Cistrome BETA
```
writeBETA(x, peakID="peak", filename="path/to/betabed.bed", scoreCol="V4")
```

### getAnnoGenes()
Return genes to which regions are annoated
```
#First annotate regions to genome using ChIPseeker
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
anno<-annotatePeak(x, tssRegion=c(-3000,3000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")

#Now retrieve genes from annotation
allGenes<-getAnnoGenes(anno, feature="all", unique=T, geneCol="SYMBOL")
promoterGenes<-getAnnoGenes(anno, feature="promoter", unique=T, geneCol="SYMBOL")
```
### retRegion()
Return genomic regions associated with specific genes of interest
```
retRegion(anno, symbol="GAPDH")
```

### FeatureEnrichment()
Calculate empirical p-values for enrichment of genomic features in given dataset over random (works on Windows with WSL and bedtools installed only) - This uses ChIPseeker to annotate bedfile of interest followed by N iterations of random genomic regions to calculate empirical p-value for each genomic feature.
```
FeatureEnrichment(bedFile="path/to/bedfile.bed", species="hsapiens", N=1000, tss=c(-1000,1000))
```

### plotVenn_gr()
Plot a Venn diagram of peak overlaps
```
x<-readBed("path/to/group1.bed")
y<-readBed("path/to/group2.bed")
plotVenn_gr(x=list(group1=x, group2=y), title="", scale=T)
```

### aggrAcc()
Aggregate peak counts to gene features
```
#Example workflow:
#Read in bedfile of peaks of interest
x<-readBed("path/to/bedfile.bed") #Fourth column in bedfile contains peakIDs
#Annotate peaks using ChIPseeker's annotatePeak
anno<-annotatePeak(x, tssRegion=c(-1000,1000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db") #There will be a column in anno named 'V4' with peakIDs
#Write peaks of interest to GTF file using peakID column V4
writeGTF(x, source="Peaks_of_interest", peakID=x$V4, filename="peaksOfInterest.gtf")

#now count bam files to regions in GTF file using preferred method (recommended Rsubread featureCounts)
library(Rsubread)
bamfiles<-c("path/to/indexed_bam1.bam","path/to/indexed_bam2.bam",...,"path/to/indexed_bamN.bam")
counts<-featureCounts(bamfiles, annot.ext="peaksOfInterest.gtf", isGTFAnnotationFile=T, isPairedEnd=T)

#Now we have counts in each sample for each peak of interest. We want to combine all counts for peaks annotated to each gene (and we can subset them based on genomic feature annotations)
allCounts<-aggrAcc(counts$counts, anno, annoCol="V4", feature="all", method="sum")
promoterCounts<-aggrAcc(counts$counts, anno, annoCol="V4", feature="promoter", method="sum")
```

### customGREAT()
Run a GREAT-like analysis on genomic regions of interest using a custom gene set gmt. 
At the moment, parallelization works only on Windows
```
gmt<-"path/to/gmt.gmt" #Can be .html address or R named list object with gene sets of interest

#If using same gmt for multiple regions of interst, start by calculating coverage of gene sets, then re-use for each analysis to save time

genomicCoverage<-customGREAT(x, gmt, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db", promoter=F, tssRegion=c(-1000,1000), parallel=T, returnCov=T, genCov=NULL, gsName="customGMT", FDR=T, enr="pos", significant=T)

#Now that we have the genomic coverage for the gene sets, we can run this multiple times more quickly

x_great<-customGREAT(x, gmt, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db", promoter=F, tssRegion=c(-1000,1000), parallel=T, returnCov=F, genCov=genomicCoverage, gsName="customGMT", FDR=T, enr="pos", significant=T)

y_great<-customGREAT(y, gmt, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db", promoter=F, tssRegion=c(-1000,1000), parallel=T, returnCov=F, genCov=genomicCoverage, gsName="customGMT", FDR=T, enr="pos", significant=T)
```

### compHOMER()
Compare HOMER known motif enrichment results from 2 different analyses
```
compHomer(xfile="path/to/x/HOMER/knownResults.txt", yfile="path/to/y/HOMER/knownResults.txt", title="Comparing x and y", xlab="sample x", ylab="sample y", lab=T, spec="CTCF", numLab=5, enrichment=F, returnRes=F)
```

### Metagene Plots
Make customized metagene plots using output from deeptools computeMatrix
#### First run deeptools computeMatrix:
```
#!/bin/bash

computeMatrix scale-regions -S x.bw y.bw -R ROI.bed -b 1000 -a 3000 -o scaled_regions.mat.gz

computeMatrix reference-point -S x.bw y.bw -R ROI.bed -b 3000 -a 3000 -o reference_point.mat.gz
```

#### Now Make metagene plots
```
#Reads in the deeptools matrices and annotates each region to its corresponding gene
scaledCov<-makeCovList("path/to/scaled_regions.mat.gz", sampleTable=data.frame(sampleName=c("group1","group2"), sampleID=c("x","y")), annoBed=T, species="human")
referenceCov<-makeCovList("path/to/reference_point.mat.gz", sampleTable=data.frame(sampleName=c("group1","group2"), sampleID=c("x","y")), annoBed=T, sepecies="human")

#Optional: using gene annotations, get one region for each annotated gene, using the region with the strongest signal in the provided control sample
scaledCov<-addSym(scaledCov, rmZeroes=T, control="group1")
referenceCov<-addSym(referenceCov, rmZeroes=T, control="group1")

#Now plot the metagene profile:
plotMetaGene(scaledCov, samples=NULL, subgenes=NULL, title="Metagene Plot", cols=NULL, marks=c("-1000","TSS","TTS","+3000"), xlab="Position (bp)", ylab="Occupancy", pal="Dark2", a=100, b=100, scale=T, retDat=F)
plotMetaGene(referenceCov, sample=NULL, subgenes=NULL, title="Metagene Plot", cols=NULL, marks=c("-3000","TSS","+3000"), xlab="Position (bp)", ylab="Occupancy", pal="Dark2", a=100, b=300, scale=F, retDat=F)

#There is another function: ChIPATAC:::plotMetaGene2(), which uses the same arguments as above, but, in addition takes an argument 'eb' (must = "sd" or "se") to plot the metagene plot with error (sd=standard deviation, se= standard error)
```


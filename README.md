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

###getAnnoGenes()
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
###retRegion()
Return genomic regions associated with specific genes of interest
```
retRegion(anno, symbol="GAPDH")
```

###FeatureEnrichment()
Calculate empirical p-values for enrichment of genomic features in given dataset over random (works on Windows with WSL and bedtools installed only) - This uses ChIPseeker to annotate bedfile of interest followed by N iterations of random genomic regions to calculate empirical p-value for each genomic feature.
```
FeatureEnrichment(bedFile="path/to/bedfile.bed", species="hsapiens", N=1000, tss=c(-1000,1000))
```

###aggrAcc()
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



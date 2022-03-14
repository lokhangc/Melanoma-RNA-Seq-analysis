library(openxlsx)

# 1. Differential splicing analysis

## 1.1. Module set-up
rMATS analysis were performed on the spartan HPC at The University of Melbourne. The foss, python, samtools and bedtools modules are required in prior of rMATS.
```{r rMATS_modules, echo=FALSE, results = 'asis'}
module <- readImage('C:/Users/Rex/Desktop/rMATS set up.PNG')
display(module, method = 'raster')
```

## 1.2. Input format
#Sample's file input for rMATS can be BAM or SAM format. BAM files were used for the following analysis. rMATS is able to analyse samples with replicates. Samples replicates of the same condition are stored in a text file. An example below showing that the SA-SPLICE19-5_sorted.bam and SA-SPLICE19-6_sorted.bam are replicates that stored in the DMSO_72hr.txt file.
```{r rMATS_input, echo=FALSE,results = 'asis'}
input <- readImage('C:/Users/Rex/Desktop/rMATS file input.PNG')
display(input, method = 'raster')
```

## 1.3. Differential splicing analysis
rMATS analysis was performed under an interactive section. The function for rMATS analysis was stored in the RNASeq-MATS.py. To perform differential splicing analysis, the function is called and parameters are set.

The '--gtf' options is the annotation file, which the hg19 GTF file is used. \
The '--b1' and '--b2' option is the BAM files input, where b1 and b2 are files to be compared between two conditions. \ 
The '--od' options is to set the output directory which contains the rMATS output.\
The '--tmp' options is to set the temporary directory which contains tracking information of the rMATS analysis.\
The '--t' options is to specify if the samples are paired-end or single-end.\
The '--readLength' options is the length of each read. Samtools can be used to find the read length of each reads.\
The '--variable-read-length' options is to allow reads with length different to --readLength also be processed.\
The '--anchorLength' options sets the minimum number of mapped bases (of each side of the splice junction) required for a novel mapped splice junction.\
The '--nthread' options assign the number of threads of the CPU to run rMATS.\
The '--cstat' options set the cutoff of splicing difference.A 0.0001 for 0.01% difference.\
The '--novelSS' options enable to detection of novel splice sites.
```{r Differential_splicing_analysis, echo=FALSE,results = 'asis'}
rMATS <- readImage('C:/Users/Rex/Desktop/rMATS code.PNG')
display(rMATS, method = 'raster')
```

## 1.4. Output files
rMATS detects five different alternative splicing events (AS_events), including alternative 3' splice site (A3SS), alternative 3' splice site (A5SS),  mutually exclusive exon (MXE), retained intron (RI) and skipped exon (SE). Seven files are generated for each splicing events. There are thirty five files in total.

The '[AS_Event].MATS.JC.txt' contains output including only reads that span junctions defined by rMATS (Junction Counts).\
The '[AS_Event].MATS.JCEC.txt' contains output including both reads that span junctions defined by rMATS (Junction Counts) and reads that do not cross an exon boundary (Exon Counts).\
The 'fromGTF.[AS_Event].txt' contains alternative splicing (AS) events derived from GTF and RNA.\
The 'fromGTF.novelJunction.[AS_Event].txt' contains alternative splicing (AS) events which were identified only after considering the RNA.\ This does not include events with an unannotated splice site.\
The 'fromGTF.novelSpliceSite.[AS_Event].txt' contains only those events which include an unannotated splice site. Only relevant if --novelSS is enabled.\
The 'JC.raw.input.[AS_Event].txt' contains event counts including only reads that span junctions defined by rMATS (Junction Counts).\
The 'JCEC.raw.input.[AS_Event].txt' contains event counts including both reads that span junctions defined by rMATS (Junction Counts) and reads that do not cross an exon boundary (Exon Counts).\
```{r rMATS_output_files, echo=FALSE,results = 'asis'}
output <- readImage('C:/Users/Rex/Desktop/rMATS output files.PNG')
display(output, method = 'raster')
```

## 1.5. Output explanation

### Headers
For detected splicing events, some of the file headers are shared.

'ID' is the the rMATS event id.\
'GeneID' is the gene id.\
'geneSymbol' is the gene name.\
'chr' is the chromosome location.\
'strand' is the strand of the gene.\
'IJC_SAMPLE_1' is the inclusion counts for sample 1. Replicates are comma separated.\
'SJC_SAMPLE_1' is the skipping counts for sample 1.\
'IJC_SAMPLE_2' is the inclusion counts for sample 2. \
'SJC_SAMPLE_2' is the skipping counts for sample 2. \
'IncFormLen' is the length of inclusion form, used for normalization.\
'SkipFormLen' is the length of skipping form, used for normalization.\
'PValue' is the significance of splicing difference between the two sample groups. (Only available if the statistical model is on).\
'FDR' is the False Discovery Rate calculated from p-value. (Only available if statistical model is on).\
'IncLevel1' is the inclusion level for sample 1. Calculated from normalized counts.\
'IncLevel2' is the Inclusion level for sample 2. Calculated from normalized counts.\
'IncLevelDifference' is the average(IncLevel1) - average(IncLevel2).\

### A3SS
Headers distinct under A3SS detection are explained in the following screeshot.

'longExonStart_0base' = longExonStart\
'longExonEnd'	= longExonEnd\
'shortES' = shortExonStart\
'shortEE'	= shortExonEnd\
'flankingES' = flankingExonStart\
'flankingEE' = flankingExonEnd\
```{r A3SS,echo=FALSE,results = 'asis'}
A3SS <- readImage('C:/Users/Rex/Desktop/A3SS.PNG')
display(A3SS, method = 'raster')
```

### A5SS
Headers distinct under A5SS detection are explained in the following screeshot.It is similar to A3SS.
```{r A5SS,echo=FALSE,results = 'asis'}
A5SS <- readImage('C:/Users/Rex/Desktop/A5SS.PNG')
display(A5SS, method = 'raster')
```

### MXE
Headers distinct under MXE detection are explained in the following screeshot.

'1stExonStart_0base' = firstExonStart\
'1stExonEnd' = firstExonEnd\
'2ndExonStart_0base' = secondExonStart\
'2ndExonEnd' = secondExonEnd\
'upstreamES' = upstreamExonStart\
'upstreamEE' = upstreamExonEnd\
'downstreamES' = downstreamExonStart\
'downstreamEE' = downstreatExonEnd\
```{r MXE,echo=FALSE,results = 'asis'}
MXE <- readImage('C:/Users/Rex/Desktop/MXE.PNG')
display(MXE, method = 'raster')
```

### RI
Headers distinct under RI detection are explained in the following screeshot.

'riExonStart_0base' = riExonStart\
'riExonEnd' = riExonEnd\
'upstreamES' = upstreamExonStart\
'upstreamEE' = upstreamExonEnd\
'downstreamES' = downstreamExonStart\
'downstreamEE' = downstreamExonEnd\
```{r RI,echo=FALSE,results = 'asis'}
RI <- readImage('C:/Users/Rex/Desktop/RI.PNG')
display(RI, method = 'raster')
```

### SE
Headers distinct under SE detection are explained in the following screeshot.

'exonStart_0base' = exonStart\
'exonEnd'	= exonEnd
'upstreamES' = firstFlankingExonStart\
'upstreamEE' = firstFlankingExonEnd\
'downstreamES' = secondFlankingExonStart\
'downstreamEE' = secondFlankingExonEnd\
```{r SE,echo=FALSE,results = 'asis'}
SE <- readImage('C:/Users/Rex/Desktop/SE.PNG')
display(SE, method = 'raster')
```

# 2. Analysing rMATS results
```{r, message=FALSE,warning=FALSE,results = 'asis'}
library(ggplot2)
library(rtracklayer)
library(enrichR)
library(limma)
```

## 2.1. Novel Junction
Input the detected novel junction by splicing events and by samples comparisons. Below shows the A3SS novel junction of different sample comparisons.
```{r novel_junction_input,results = 'asis'}
#A3SS novel junction
CDK4i_72hr_DMSO_72hr_novelJ_A3SS <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelJunction.A3SS.txt')
PRMT5i_72hr_DMSO_72hr_novelJ_A3SS <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelJunction.A3SS.txt')
CDK4i_6d_DMSO_72hr_novelJ_A3SS <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelJunction.A3SS.txt')

```
```{r echo=FALSE}
#A5SS novel junction
CDK4i_72hr_DMSO_72hr_novelJ_A5SS <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelJunction.A5SS.txt')
PRMT5i_72hr_DMSO_72hr_novelJ_A5SS <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelJunction.A5SS.txt')
CDK4i_6d_DMSO_72hr_novelJ_A5SS <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelJunction.A5SS.txt')

#MXE novel junction
CDK4i_72hr_DMSO_72hr_novelJ_MXE <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelJunction.MXE.txt')
PRMT5i_72hr_DMSO_72hr_novelJ_MXE <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelJunction.MXE.txt')
CDK4i_6d_DMSO_72hr_novelJ_MXE <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelJunction.MXE.txt')

#RI novel junction
CDK4i_72hr_DMSO_72hr_novelJ_RI <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelJunction.RI.txt')
PRMT5i_72hr_DMSO_72hr_novelJ_RI <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelJunction.RI.txt')
CDK4i_6d_DMSO_72hr_novelJ_RI <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelJunction.RI.txt')

#SE novel junction
CDK4i_72hr_DMSO_72hr_novelJ_SE <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelJunction.SE.txt')
PRMT5i_72hr_DMSO_72hr_novelJ_SE <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelJunction.SE.txt')
CDK4i_6d_DMSO_72hr_novelJ_SE <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelJunction.SE.txt')
```

The novel junctions are then summerised by the gene. 
```{r novel_junction_count_by_gene,results = 'asis'}
#A3SS novel junction count by gene
CDK4i_72hr_DMSO_72hr_novelJ_A3SS_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelJ_A3SS[, 3]))
PRMT5i_72hr_DMSO_72hr_novelJ_A3SS_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelJ_A3SS[, 3]))
CDK4i_6d_DMSO_72hr_novelJ_A3SS_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelJ_A3SS[, 3]))

```
```{r echo=FALSE}
#A5SS novel junction count by gene
CDK4i_72hr_DMSO_72hr_novelJ_A5SS_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelJ_A5SS[, 3]))
PRMT5i_72hr_DMSO_72hr_novelJ_A5SS_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelJ_A5SS[, 3]))
CDK4i_6d_DMSO_72hr_novelJ_A5SS_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelJ_A5SS[, 3]))

#MXE novel junction
CDK4i_72hr_DMSO_72hr_novelJ_MXE_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelJ_MXE[, 3]))
PRMT5i_72hr_DMSO_72hr_novelJ_MXE_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelJ_MXE[, 3]))
CDK4i_6d_DMSO_72hr_novelJ_MXE_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelJ_MXE[, 3]))

#RI novel junction
CDK4i_72hr_DMSO_72hr_novelJ_RI_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelJ_RI[, 3]))
PRMT5i_72hr_DMSO_72hr_novelJ_RI_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelJ_RI[, 3]))
CDK4i_6d_DMSO_72hr_novelJ_RI_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelJ_RI[, 3]))

#SE novel junction
CDK4i_72hr_DMSO_72hr_novelJ_SE_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelJ_SE[, 3]))
PRMT5i_72hr_DMSO_72hr_novelJ_SE_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelJ_SE[, 3]))
CDK4i_6d_DMSO_72hr_novelJ_SE_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelJ_SE[, 3]))
```

Calculate the ratio between the number of novel junction and novel junction genes. The ratio shows helps deciding the number of novel junction detected per gene in average.
```{r norvelJ_ratio,results = 'asis'}
novelJ_ratio_A3SS <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelJ_A3SS) / nrow(CDK4i_72hr_DMSO_72hr_novelJ_A3SS_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelJ_A3SS) / nrow(PRMT5i_72hr_DMSO_72hr_novelJ_A3SS_gene),
nrow(CDK4i_6d_DMSO_72hr_novelJ_A3SS) / nrow(CDK4i_6d_DMSO_72hr_novelJ_A3SS_gene)
)

```
```{r echo=FALSE}
novelJ_ratio_A5SS <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelJ_A5SS) / nrow(CDK4i_72hr_DMSO_72hr_novelJ_A5SS_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelJ_A5SS) / nrow(PRMT5i_72hr_DMSO_72hr_novelJ_A5SS_gene),
nrow(CDK4i_6d_DMSO_72hr_novelJ_A5SS) / nrow(CDK4i_6d_DMSO_72hr_novelJ_A5SS_gene)
)
novelJ_ratio_MXE <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelJ_MXE) / nrow(CDK4i_72hr_DMSO_72hr_novelJ_MXE_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelJ_MXE) / nrow(PRMT5i_72hr_DMSO_72hr_novelJ_MXE_gene),
nrow(CDK4i_6d_DMSO_72hr_novelJ_MXE) / nrow(CDK4i_6d_DMSO_72hr_novelJ_MXE_gene)
)
novelJ_ratio_RI <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelJ_RI) / nrow(CDK4i_72hr_DMSO_72hr_novelJ_RI_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelJ_RI) / nrow(PRMT5i_72hr_DMSO_72hr_novelJ_RI_gene),
nrow(CDK4i_6d_DMSO_72hr_novelJ_RI) / nrow(CDK4i_6d_DMSO_72hr_novelJ_RI_gene)
)
novelJ_ratio_SE <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelJ_SE) / nrow(CDK4i_72hr_DMSO_72hr_novelJ_SE_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelJ_SE) / nrow(PRMT5i_72hr_DMSO_72hr_novelJ_SE_gene),
nrow(CDK4i_6d_DMSO_72hr_novelJ_SE) / nrow(CDK4i_6d_DMSO_72hr_novelJ_SE_gene)
)
```
```{r novelJ_ratio_matrix,results = 'asis'}
novelJ_ratio <- data.frame(
'A3SS' = novelJ_ratio_A3SS,
'A5SS' = novelJ_ratio_A5SS,
'MXE' = novelJ_ratio_MXE,
'RI' = novelJ_ratio_RI,
'SE' = novelJ_ratio_SE
)
row.names(novelJ_ratio) <-
c('CDK4i_72hr_DMSO_72hr',
'PRMT5i_72hr_DMSO_72hr',
'CDK4i_6d_DMSO_72hr')
kable(novelJ_ratio)
```

Creating a novel junction count matrix.
```{r novel_junction_count_matrix,results = 'asis'}
comparison <- c(
rep("CDK4i_72hr_DMSO_72hr" , 5),
rep("PRMT5i_72hr_DMSO_72hr" , 5),
rep("CDK4i_6d_DMSO_72hr" , 5)
)
splice_event <- c(rep(c('A3SS', 'A5SS', 'MXE', 'RI', 'SE'), 3))
Number_of_noveljunction <-
c(
dim(CDK4i_72hr_DMSO_72hr_novelJ_A3SS)[1],
dim(CDK4i_72hr_DMSO_72hr_novelJ_A5SS)[1],
dim(CDK4i_72hr_DMSO_72hr_novelJ_MXE)[1],
dim(CDK4i_72hr_DMSO_72hr_novelJ_RI)[1],
dim(CDK4i_72hr_DMSO_72hr_novelJ_SE)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelJ_A3SS)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelJ_A5SS)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelJ_MXE)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelJ_RI)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelJ_SE)[1],
dim(CDK4i_6d_DMSO_72hr_novelJ_A3SS)[1],
dim(CDK4i_6d_DMSO_72hr_novelJ_A5SS)[1],
dim(CDK4i_6d_DMSO_72hr_novelJ_MXE)[1],
dim(CDK4i_6d_DMSO_72hr_novelJ_RI)[1],
dim(CDK4i_6d_DMSO_72hr_novelJ_SE)[1]
)

novelJunction <-
data.frame(comparison, splice_event, Number_of_noveljunction)

kable(novelJunction)
```

Visualize the novel junction count with barplot.
```{r novel_junction_count_plot,results = 'asis'}
ggplot(novelJunction,aes(fill=splice_event,y=Number_of_noveljunction,x=comparison))+geom_bar(position="dodge", stat="identity")+scale_x_discrete(guide = guide_axis(n.dodge=3))+ggtitle('Number of novel junction')
```

## 2.2. Novel splice site
Input the detected novel splice site by splicing events and by samples comparisons. Below shows the A3SS novel splice site of different sample comparisons.
```{r novel_splice_site_input,results = 'asis'}
#A3SS novel splice site
CDK4i_72hr_DMSO_72hr_novelSpliceSite_A3SS <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.A3SS.txt')
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A3SS <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.A3SS.txt')
CDK4i_6d_DMSO_72hr_novelSpliceSite_A3SS <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelSpliceSite.A3SS.txt')

```
```{r echo=FALSE}
#A5SS novel splice site
CDK4i_72hr_DMSO_72hr_novelSpliceSite_A5SS <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.A5SS.txt')
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A5SS <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.A5SS.txt')
CDK4i_6d_DMSO_72hr_novelSpliceSite_A5SS <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelSpliceSite.A5SS.txt')

#MXE novel splice site
CDK4i_72hr_DMSO_72hr_novelSpliceSite_MXE <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.MXE.txt')
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_MXE <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.MXE.txt')
CDK4i_6d_DMSO_72hr_novelSpliceSite_MXE <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelSpliceSite.MXE.txt')

#RI novel splice site
CDK4i_72hr_DMSO_72hr_novelSpliceSite_RI <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.RI.txt')
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_RI <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.RI.txt')
CDK4i_6d_DMSO_72hr_novelSpliceSite_RI <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelSpliceSite.RI.txt')

#SE novel splice site
CDK4i_72hr_DMSO_72hr_novelSpliceSite_SE <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.SE.txt')
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_SE <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/fromGTF.novelSpliceSite.SE.txt')
CDK4i_6d_DMSO_72hr_novelSpliceSite_SE <-
read.delim('Results/CDK4i_6d_DMSO_72hr/fromGTF.novelSpliceSite.SE.txt')
```

The novel splice site are then summerised by the gene. 
```{r novel_splice_site_count_by_gene,results = 'asis'}
#A3SS novel splice site count by gene
CDK4i_72hr_DMSO_72hr_novelSpliceSite_A3SS_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelSpliceSite_A3SS[, 3]))
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A3SS_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A3SS[, 3]))
CDK4i_6d_DMSO_72hr_novelSpliceSite_A3SS_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelSpliceSite_A3SS[, 3]))

```
```{r echo=FALSE}
#A5SS novel splice site count by gene
CDK4i_72hr_DMSO_72hr_novelSpliceSite_A5SS_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelSpliceSite_A5SS[, 3]))
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A5SS_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A5SS[, 3]))
CDK4i_6d_DMSO_72hr_novelSpliceSite_A5SS_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelSpliceSite_A5SS[, 3]))

#MXE novel splice site count by gene
CDK4i_72hr_DMSO_72hr_novelSpliceSite_MXE_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelSpliceSite_MXE[, 3]))
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_MXE_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_MXE[, 3]))
CDK4i_6d_DMSO_72hr_novelSpliceSite_MXE_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelSpliceSite_MXE[, 3]))

#RI novel splice site count by gene
CDK4i_72hr_DMSO_72hr_novelSpliceSite_RI_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelSpliceSite_RI[, 3]))
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_RI_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_RI[, 3]))
CDK4i_6d_DMSO_72hr_novelSpliceSite_RI_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelSpliceSite_RI[, 3]))

#SE novel splice site count by gene
CDK4i_72hr_DMSO_72hr_novelSpliceSite_SE_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_novelSpliceSite_SE[, 3]))
PRMT5i_72hr_DMSO_72hr_novelSpliceSite_SE_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_SE[, 3]))
CDK4i_6d_DMSO_72hr_novelSpliceSite_SE_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_novelSpliceSite_SE[, 3]))
```

Calculate the ratio between the number of novel splice site and novel splice site genes. The ratio shows helps deciding the number of novel splice site detected per gene in average.
```{r novelSpliceSite_ratio,results = 'asis'}
novelSpliceSite_ratio_A3SS <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_A3SS) / nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_A3SS_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A3SS) / nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A3SS_gene),
nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_A3SS) / nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_A3SS_gene)
)

```
```{r echo=FALSE}
novelSpliceSite_ratio_A5SS <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_A5SS) / nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_A5SS_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A5SS) / nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A5SS_gene),
nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_A5SS) / nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_A5SS_gene)
)
novelSpliceSite_ratio_MXE <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_MXE) / nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_MXE_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_MXE) / nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_MXE_gene),
nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_MXE) / nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_MXE_gene)
)
novelSpliceSite_ratio_RI <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_RI) / nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_RI_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_RI) / nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_RI_gene),
nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_RI) / nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_RI_gene)
)
novelSpliceSite_ratio_SE <-
c(
nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_SE) / nrow(CDK4i_72hr_DMSO_72hr_novelSpliceSite_SE_gene),
nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_SE) / nrow(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_SE_gene),
nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_SE) / nrow(CDK4i_6d_DMSO_72hr_novelSpliceSite_SE_gene)
)
```
```{r novelSpliceSite_ratio_matrix,results = 'asis'}
novelSpliceSite_ratio <- data.frame(
'A3SS' = novelSpliceSite_ratio_A3SS,
'A5SS' = novelSpliceSite_ratio_A5SS,
'MXE' = novelSpliceSite_ratio_MXE,
'RI' = novelSpliceSite_ratio_RI,
'SE' = novelSpliceSite_ratio_SE
)
row.names(novelSpliceSite_ratio) <-
c('CDK4i_72hr_DMSO_72hr',
'PRMT5i_72hr_DMSO_72hr',
'CDK4i_6d_DMSO_72hr')
kable(novelSpliceSite_ratio)
```
Creating a novel splice site count matrix.
```{r novel_splice_site_count_matrix,results = 'asis'}
Number_of_novelSpliceSite <-
c(
dim(CDK4i_72hr_DMSO_72hr_novelSpliceSite_A3SS)[1],
dim(CDK4i_72hr_DMSO_72hr_novelSpliceSite_A5SS)[1],
dim(CDK4i_72hr_DMSO_72hr_novelSpliceSite_MXE)[1],
dim(CDK4i_72hr_DMSO_72hr_novelSpliceSite_RI)[1],
dim(CDK4i_72hr_DMSO_72hr_novelSpliceSite_SE)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A3SS)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_A5SS)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_MXE)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_RI)[1],
dim(PRMT5i_72hr_DMSO_72hr_novelSpliceSite_SE)[1],
dim(CDK4i_6d_DMSO_72hr_novelSpliceSite_A3SS)[1],
dim(CDK4i_6d_DMSO_72hr_novelSpliceSite_A5SS)[1],
dim(CDK4i_6d_DMSO_72hr_novelSpliceSite_MXE)[1],
dim(CDK4i_6d_DMSO_72hr_novelSpliceSite_RI)[1],
dim(CDK4i_6d_DMSO_72hr_novelSpliceSite_SE)[1]
)

novelSpliceSite <-
data.frame(comparison, splice_event, Number_of_novelSpliceSite)

kable(novelSpliceSite)

```

Visualize the novel splice site count with barplot.
```{r novel_splice_site_count_plot,results = 'asis'}
ggplot(novelSpliceSite,aes(fill=splice_event,y=Number_of_novelSpliceSite,x=comparison))+geom_bar(position="dodge",stat="identity")+scale_x_discrete(guide = guide_axis(n.dodge=3))+ggtitle('Number of novel splice site')
```

## 2.3. Splicing events
Input the detected splicing events by samples comparisons. Below read the A3SS result of different sample comparisons.
```{r Splicing_events,results = 'asis'}
#A3SS
CDK4i_72hr_DMSO_72hr_A3SS <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/A3SS.MATS.JC.txt')
PRMT5i_72hr_DMSO_72hr_A3SS <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/A3SS.MATS.JC.txt')
CDK4i_6d_DMSO_72hr_A3SS <-
read.delim('Results/CDK4i_6d_DMSO_72hr/A3SS.MATS.JC.txt')

```
```{r echo=FALSE}
#A5SS
CDK4i_72hr_DMSO_72hr_A5SS <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/A5SS.MATS.JC.txt')
PRMT5i_72hr_DMSO_72hr_A5SS <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/A5SS.MATS.JC.txt')
CDK4i_6d_DMSO_72hr_A5SS <-
read.delim('Results/CDK4i_6d_DMSO_72hr/A5SS.MATS.JC.txt')

#MXE
CDK4i_72hr_DMSO_72hr_MXE <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/MXE.MATS.JC.txt')
PRMT5i_72hr_DMSO_72hr_MXE <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/MXE.MATS.JC.txt')
CDK4i_6d_DMSO_72hr_MXE <-
read.delim('Results/CDK4i_6d_DMSO_72hr/MXE.MATS.JC.txt')

#RI
CDK4i_72hr_DMSO_72hr_RI <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/RI.MATS.JC.txt')
PRMT5i_72hr_DMSO_72hr_RI <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/RI.MATS.JC.txt')
CDK4i_6d_DMSO_72hr_RI <-
read.delim('Results/CDK4i_6d_DMSO_72hr/RI.MATS.JC.txt')

#SE
CDK4i_72hr_DMSO_72hr_SE <-
read.delim('Results/CDK4i_72hr_DMSO_72hr/SE.MATS.JC.txt')
PRMT5i_72hr_DMSO_72hr_SE <-
read.delim('Results/PRMT5i_72hr_DMSO_72hr/SE.MATS.JC.txt')
CDK4i_6d_DMSO_72hr_SE <-
read.delim('Results/CDK4i_6d_DMSO_72hr/SE.MATS.JC.txt')
```

Filter significant splicing events with false discovery rate < 0.05.
```{r Significant_splicing_events,results = 'asis'}
#significant A3SS
CDK4i_72hr_DMSO_72hr_A3SS_sig <-
CDK4i_72hr_DMSO_72hr_A3SS[CDK4i_72hr_DMSO_72hr_A3SS$FDR < 0.05, ]
PRMT5i_72hr_DMSO_72hr_A3SS_sig <-
PRMT5i_72hr_DMSO_72hr_A3SS[PRMT5i_72hr_DMSO_72hr_A3SS$FDR < 0.05, ]
CDK4i_6d_DMSO_72hr_A3SS_sig <-
CDK4i_6d_DMSO_72hr_A3SS[CDK4i_6d_DMSO_72hr_A3SS$FDR < 0.05, ]

```
```{r echo=FALSE}
#significant A5SS
CDK4i_72hr_DMSO_72hr_A5SS_sig <-
CDK4i_72hr_DMSO_72hr_A5SS[CDK4i_72hr_DMSO_72hr_A5SS$FDR < 0.05, ]
PRMT5i_72hr_DMSO_72hr_A5SS_sig <-
PRMT5i_72hr_DMSO_72hr_A5SS[PRMT5i_72hr_DMSO_72hr_A5SS$FDR < 0.05, ]
CDK4i_6d_DMSO_72hr_A5SS_sig <-
CDK4i_6d_DMSO_72hr_A5SS[CDK4i_6d_DMSO_72hr_A5SS$FDR < 0.05, ]

#significant MXE
CDK4i_72hr_DMSO_72hr_MXE_sig <-
CDK4i_72hr_DMSO_72hr_MXE[CDK4i_72hr_DMSO_72hr_MXE$FDR < 0.05, ]
PRMT5i_72hr_DMSO_72hr_MXE_sig <-
PRMT5i_72hr_DMSO_72hr_MXE[PRMT5i_72hr_DMSO_72hr_MXE$FDR < 0.05, ]
CDK4i_6d_DMSO_72hr_MXE_sig <-
CDK4i_6d_DMSO_72hr_MXE[CDK4i_6d_DMSO_72hr_MXE$FDR < 0.05, ]

#significant RI
CDK4i_72hr_DMSO_72hr_RI_sig <-
CDK4i_72hr_DMSO_72hr_RI[CDK4i_72hr_DMSO_72hr_RI$FDR < 0.05, ]
PRMT5i_72hr_DMSO_72hr_RI_sig <-
PRMT5i_72hr_DMSO_72hr_RI[PRMT5i_72hr_DMSO_72hr_RI$FDR < 0.05, ]
CDK4i_6d_DMSO_72hr_RI_sig <-
CDK4i_6d_DMSO_72hr_RI[CDK4i_6d_DMSO_72hr_RI$FDR < 0.05, ]

#significant SE
CDK4i_72hr_DMSO_72hr_SE_sig <-
CDK4i_72hr_DMSO_72hr_SE[CDK4i_72hr_DMSO_72hr_SE$FDR < 0.05, ]
PRMT5i_72hr_DMSO_72hr_SE_sig <-
PRMT5i_72hr_DMSO_72hr_SE[PRMT5i_72hr_DMSO_72hr_SE$FDR < 0.05, ]
CDK4i_6d_DMSO_72hr_SE_sig <-
CDK4i_6d_DMSO_72hr_SE[CDK4i_6d_DMSO_72hr_SE$FDR < 0.05, ]
```

The splicing events are then summerised by the gene.
```{r Splicing_event_count_by_gene,results = 'asis'}
#A3SS count by gene
CDK4i_72hr_DMSO_72hr_A3SS_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_A3SS_sig[, 3]))
PRMT5i_72hr_DMSO_72hr_A3SS_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_A3SS_sig[, 3]))
CDK4i_6d_DMSO_72hr_A3SS_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_A3SS_sig[, 3]))

```
```{r echo=FALSE}
#A3SS count by gene
CDK4i_72hr_DMSO_72hr_A5SS_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_A5SS_sig[, 3]))
PRMT5i_72hr_DMSO_72hr_A5SS_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_A5SS_sig[, 3]))
CDK4i_6d_DMSO_72hr_A5SS_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_A5SS_sig[, 3]))

#MXE count by gene
CDK4i_72hr_DMSO_72hr_MXE_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_MXE_sig[, 3]))
PRMT5i_72hr_DMSO_72hr_MXE_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_MXE_sig[, 3]))
CDK4i_6d_DMSO_72hr_MXE_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_MXE_sig[, 3]))

#RI count by gene
CDK4i_72hr_DMSO_72hr_RI_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_RI_sig[, 3]))
PRMT5i_72hr_DMSO_72hr_RI_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_RI_sig[, 3]))
CDK4i_6d_DMSO_72hr_RI_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_RI_sig[, 3]))

#SE
CDK4i_72hr_DMSO_72hr_SE_gene <-
data.frame(table(CDK4i_72hr_DMSO_72hr_SE_sig[, 3]))
PRMT5i_72hr_DMSO_72hr_SE_gene <-
data.frame(table(PRMT5i_72hr_DMSO_72hr_SE_sig[, 3]))
CDK4i_6d_DMSO_72hr_SE_gene <-
data.frame(table(CDK4i_6d_DMSO_72hr_SE_sig[, 3]))
```

Calculate the ratio between the number of splice events and spliced genes. The ratio shows helps deciding the number of splice events detected per gene in average.
```{r splice_event_ratio,results = 'asis'}
event_ratio_A3SS <-
c(
nrow(CDK4i_72hr_DMSO_72hr_A3SS_sig) / nrow(CDK4i_72hr_DMSO_72hr_A3SS_gene),
nrow(PRMT5i_72hr_DMSO_72hr_A3SS_sig) / nrow(PRMT5i_72hr_DMSO_72hr_A3SS_gene),
nrow(CDK4i_6d_DMSO_72hr_A3SS_sig) / nrow(CDK4i_6d_DMSO_72hr_A3SS_gene)
)

```
```{r echo=FALSE}
event_ratio_A5SS <-
c(
nrow(CDK4i_72hr_DMSO_72hr_A5SS_sig) / nrow(CDK4i_72hr_DMSO_72hr_A5SS_gene),
nrow(PRMT5i_72hr_DMSO_72hr_A5SS_sig) / nrow(PRMT5i_72hr_DMSO_72hr_A5SS_gene),
nrow(CDK4i_6d_DMSO_72hr_A5SS_sig) / nrow(CDK4i_6d_DMSO_72hr_A5SS_gene)
)
event_ratio_MXE <-
c(
nrow(CDK4i_72hr_DMSO_72hr_MXE_sig) / nrow(CDK4i_72hr_DMSO_72hr_MXE_gene),
nrow(PRMT5i_72hr_DMSO_72hr_MXE_sig) / nrow(PRMT5i_72hr_DMSO_72hr_MXE_gene),
nrow(CDK4i_6d_DMSO_72hr_MXE_sig) / nrow(CDK4i_6d_DMSO_72hr_MXE_gene)
)
event_ratio_RI <-
c(
nrow(CDK4i_72hr_DMSO_72hr_RI_sig) / nrow(CDK4i_72hr_DMSO_72hr_RI_gene),
nrow(PRMT5i_72hr_DMSO_72hr_RI_sig) / nrow(PRMT5i_72hr_DMSO_72hr_RI_gene),
nrow(CDK4i_6d_DMSO_72hr_RI_sig) / nrow(CDK4i_6d_DMSO_72hr_RI_gene)
)
event_ratio_SE <-
c(
nrow(CDK4i_72hr_DMSO_72hr_SE_sig) / nrow(CDK4i_72hr_DMSO_72hr_SE_gene),
nrow(PRMT5i_72hr_DMSO_72hr_SE_sig) / nrow(PRMT5i_72hr_DMSO_72hr_SE_gene),
nrow(CDK4i_6d_DMSO_72hr_SE_sig) / nrow(CDK4i_6d_DMSO_72hr_SE_gene)
)
```
```{r event_ratio_matrix,results = 'asis'}
event_ratio <- data.frame(
'A3SS' = event_ratio_A3SS,
'A5SS' = event_ratio_A5SS,
'MXE' = event_ratio_MXE,
'RI' = event_ratio_RI,
'SE' = event_ratio_SE
)
row.names(event_ratio) <-
c('CDK4i_72hr_DMSO_72hr',
'PRMT5i_72hr_DMSO_72hr',
'CDK4i_6d_DMSO_72hr')
kable(event_ratio)
```
Creating a splice events count matrix.
```{r Splice_event_count_matrix,results = 'asis'}
Number_of_spliced_event <- c(
dim(CDK4i_72hr_DMSO_72hr_A3SS_sig)[1],
dim(CDK4i_72hr_DMSO_72hr_A5SS_sig)[1],
dim(CDK4i_72hr_DMSO_72hr_MXE_sig)[1],
dim(CDK4i_72hr_DMSO_72hr_RI_sig)[1],
dim(CDK4i_72hr_DMSO_72hr_SE_sig)[1],
dim(PRMT5i_72hr_DMSO_72hr_A3SS_sig)[1],
dim(PRMT5i_72hr_DMSO_72hr_A5SS_sig)[1],
dim(PRMT5i_72hr_DMSO_72hr_MXE_sig)[1],
dim(PRMT5i_72hr_DMSO_72hr_RI_sig)[1],
dim(PRMT5i_72hr_DMSO_72hr_SE_sig)[1],
dim(CDK4i_6d_DMSO_72hr_A3SS_sig)[1],
dim(CDK4i_6d_DMSO_72hr_A5SS_sig)[1],
dim(CDK4i_6d_DMSO_72hr_MXE_sig)[1],
dim(CDK4i_6d_DMSO_72hr_RI_sig)[1],
dim(CDK4i_6d_DMSO_72hr_SE_sig)[1]
)

SpliceEvent <-
data.frame(comparison, splice_event, Number_of_spliced_event)

kable(SpliceEvent)
```

Visualize the splice events count with barplot.
```{r Splice_event_count_plot,results = 'asis'}
ggplot(SpliceEvent,aes(fill=splice_event,y=Number_of_spliced_event,x=comparison))+geom_bar(position="dodge", stat="identity")+
scale_x_discrete(guide = guide_axis(n.dodge=3))+ggtitle('Number of splicing events')
```

## 2.4. Compare splicing events
Since there are more splice events detected under CDK4i_6d than CDK4i_72hr, thus, the splicing events under CDK4i_6d are compared to PRMT5i_72hr. Below compares the A3SS result of PRMT5i_72hr and CDK4i_6d.
```{r compare_splicing_events,results = 'asis'}
#A3SS
DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A3SS <-
merge(
PRMT5i_72hr_DMSO_72hr_A3SS_sig[, -c(1:2, 12:22)],
CDK4i_6d_DMSO_72hr_A3SS_sig[, -c(1:2, 12:22)],
by = 1:9,
all = TRUE
)
colnames(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A3SS)[10:11] <-
c('PRMT5i', 'CDK4i_6d')

#MXE
DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_MXE <-
merge(
PRMT5i_72hr_DMSO_72hr_MXE_sig[, -c(1:2, 14:24)],
CDK4i_6d_DMSO_72hr_MXE_sig[, -c(1:2, 14:24)],
by = 1:11,
all = TRUE
)
colnames(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_MXE)[12:13] <-
c('PRMT5i', 'CDK4i_6d')


```
```{r echo=FALSE}
#A5SS
DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A5SS <-
merge(
PRMT5i_72hr_DMSO_72hr_A5SS_sig[, -c(1:2, 12:22)],
CDK4i_6d_DMSO_72hr_A5SS_sig[, -c(1:2, 12:22)],
by = 1:9,
all = TRUE
)
colnames(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A5SS)[10:11] <-
c('PRMT5i', 'CDK4i_6d')

#RI
DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_RI <-
merge(
PRMT5i_72hr_DMSO_72hr_RI_sig[, -c(1:2, 12:22)],
CDK4i_6d_DMSO_72hr_RI_sig[, -c(1:2, 12:22)],
by = 1:9,
all = TRUE
)
colnames(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_RI)[10:11] <-
c('PRMT5i', 'CDK4i_6d')

#SE
DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_SE <-
merge(
PRMT5i_72hr_DMSO_72hr_SE_sig[, -c(1:2, 12:22)],
CDK4i_6d_DMSO_72hr_SE_sig[, -c(1:2, 12:22)],
by = 1:9,
all = TRUE
)
colnames(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_SE)[10:11] <-
c('PRMT5i', 'CDK4i_6d')
```

Visualize the splice events count between PRMT5i_72hr and CDK4i_6d with venn diagrams.
```{r plot_splicig_events_between_CDK4i_&_PRMT5i,results = 'asis'}
par(mfrow = c(2, 3))
vennDiagram(ifelse(is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A3SS[, 10:11]), 0, 1), cex = 1.0, main =
'A3SS')
vennDiagram(ifelse(is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A5SS[, 10:11]), 0, 1), cex = 1.0, main =
'A5SS')
vennDiagram(ifelse(is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_MXE[, 12:13]), 0, 1), cex = 1.0, main =
'MXE')
vennDiagram(ifelse(is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_RI[, 10:11]), 0, 1), cex = 1.0, main =
'RI')
vennDiagram(ifelse(is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_SE[, 10:11]), 0, 1), cex = 1.0, main =
'SE')
```

There are much more SE detected than the other splicing events. The rMATS output dose not include the information of exon number and transcript ID of a particular splicing events. Below the SE events are annotated with the transcript ID and exon number.
These information are stored in the GTF file.
```{r annotate_splicind_exon,results = 'asis'}
GTF = 'GTF/Homo_sapiens.GRCh37.87.gtf'
GTF_read <- readGFF(GTF)
GTF_exon <- GTF_read[GTF_read$type == 'exon', ]
GTF_exon <- GTF_exon[, c(4, 5, 11, 16, 22)]
```

A function that annotate SE.
```{r SE_annotation_function,results = 'asis'}
SE_annotation <- function(rMATS_result) {
rMATS_result$exonStart_0base <- rMATS_result$exonStart_0base + 1
rMATS_result <-
merge(
rMATS_result,
GTF_exon,
by.x = c('geneSymbol', 'exonStart_0base', 'exonEnd'),
by.y = c('gene_name', 'start', 'end'),
all.x = TRUE,
all.y = FALSE
)
rMATS_result$exonStart_0base <- rMATS_result$exonStart_0base - 1
return(rMATS_result)
}
```

The compared splicing events are stored in three files. The splicing_under_PRMT5i_&_CDK4i.xlsx contains splicing events that are commonly detected in both inhibitions; the splicing_under_CDK4i_not_PRMT5i.xlsx contains splicing events that are only detected under CDK4i; and the splicing_under_PRMT5i_not_CDK4i.xlsx contains splicing events that are only detected under PRMT5i. The SE events are annotated with transcript ID and exon number.
```{r write_distinct_&_common_splicing_events_SE_annotated,echo=FALSE}
DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_splicing <-
list(
'A3SS' = na.omit(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A3SS),
'A5SS' = na.omit(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A5SS),
'MXE' = na.omit(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_MXE),
'RI' = na.omit(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_RI),
'SE' = na.omit(SE_annotation(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_SE))
)

splicing_under_CDK4i_not_PRMT5i <-
list(
'A3SS' = DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A3SS[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A3SS$PRMT5i), ],
'A5SS' = DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A5SS[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A5SS$PRMT5i), ],
'MXE' = DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_MXE[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_MXE$PRMT5i), ],
'RI' = DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_RI[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_RI$PRMT5i), ],
'SE' = SE_annotation(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_SE[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_SE$PRMT5i), ])
)

splicing_under_PRMT5i_not_CDK4i <-
list(
'A3SS' = DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A3SS[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A3SS$CDK4i_6d), ],
'A5SS' = DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A5SS[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_A5SS$CDK4i_6d), ],
'MXE' = DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_MXE[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_MXE$CDK4i_6d), ],
'RI' = DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_RI[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_RI$CDK4i_6d), ],
'SE' = SE_annotation(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_SE[is.na(DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_SE$CDK4i_6d), ])
)

write.xlsx(
DMSO_PRMT5i_72hr_DMSO_CDK4i_6d_splicing,
"Analysis/splice_event_comparison/splicing_under_PRMT5i_&_CDK4i.xlsx"
)
write.xlsx(
splicing_under_CDK4i_not_PRMT5i,
"Analysis/splice_event_comparison/splicing_under_CDK4i_not_PRMT5i.xlsx"
)
write.xlsx(
splicing_under_PRMT5i_not_CDK4i,
"Analysis/splice_event_comparison/splicing_under_PRMT5i_not_CDK4i.xlsx"
)

```

## 2.5. Compare splicing events and differentially expressed genes
The SE are then compared to the differential gene expression result to study the relationship of alternative splicing and gene expression.The SE are annotated with transcript ID and exon number.
```{r SE_anotation,results = 'asis'}
CDK4i_72hr_DMSO_72hr_SE_sig_ann <-
SE_annotation(CDK4i_72hr_DMSO_72hr_SE_sig)
PRMT5i_72hr_DMSO_72hr_SE_sig_ann <-
SE_annotation(PRMT5i_72hr_DMSO_72hr_SE_sig)
CDK4i_6d_DMSO_72hr_SE_sig_ann <-
SE_annotation(CDK4i_6d_DMSO_72hr_SE_sig)
```

Input the limma differential gene expression result by samples comparisons.
```{r limma_DGE_result,message=FALSE, warning=FALSE,results = 'asis'}
setwd('../Differential_gene_expression/limma/')
CDK4i_72hr_DMSO_72hr_DEG <-
read.xlsx('limma_sigchange_dge.xlsx',
sheet = 1)
PRMT5i_72hr_DMSO_72hr_DEG <-
read.xlsx('limma_sigchange_dge.xlsx',
sheet = 2)
CDK4i_6d_DMSO_72hr_DEG <-
read.xlsx('limma_sigchange_dge.xlsx',
sheet = 4)
```

Annotate SE by the deferentially expressed genes.
```{r limma_DGE_result_annotation,results = 'asis'}
CDK4i_72hr_DMSO_72hr_SE_DEG <- merge(CDK4i_72hr_DMSO_72hr_SE_sig_ann,CDK4i_72hr_DMSO_72hr_DEG,by.x='geneSymbol',by.y='',all=TRUE)
PRMT5i_72hr_DMSO_72hr_SE_DEG <- merge(PRMT5i_72hr_DMSO_72hr_SE_sig_ann,PRMT5i_72hr_DMSO_72hr_DEG,by.x='geneSymbol',by.y='',all=TRUE)
CDK4i_6d_DMSO_72hr_SE_DEG <- merge(CDK4i_6d_DMSO_72hr_SE_sig_ann,CDK4i_6d_DMSO_72hr_DEG,by.x='geneSymbol',by.y='',all=TRUE)
```

The DGE_annotated_SE.xlsx contains the SE with annotated transcript ID, exon number and differentially expressed genes.
```{r write_DGE_annotated_SE, echo=FALSE}
DGE_annotated_SE <- list('CDK4i_72hr_DMSO_72hr' = CDK4i_72hr_DMSO_72hr_SE_DEG,
'PRMT5i_72hr_DMSO_72hr' = PRMT5i_72hr_DMSO_72hr_SE_DEG,
'CDK4i_6d_DMSO_72hr' = CDK4i_6d_DMSO_72hr_SE_DEG)

write.xlsx(DGE_annotated_SE,"Analysis/Annotated_SE/DGE_annotated_SE.xlsx")
```

Compare differentially spliced genes and differentially expressed genes.
```{r,results = 'asis'}
CDK4i_72hr_DMSO_72hr_DSPG_DEG <- merge(CDK4i_72hr_DMSO_72hr_SE_gene,CDK4i_72hr_DMSO_72hr_DEG,by.x='Var1',by.y='',all=TRUE)
PRMT5i_72hr_DMSO_72hr_DSPG_DEG <- merge(PRMT5i_72hr_DMSO_72hr_SE_gene,PRMT5i_72hr_DMSO_72hr_DEG,by.x='Var1',by.y='',all=TRUE)
CDK4i_6d_DMSO_72hr_DSPG_DEG <- merge(CDK4i_6d_DMSO_72hr_SE_gene,CDK4i_6d_DMSO_72hr_DEG,by.x='Var1',by.y='',all=TRUE)
```
```{r, echo=FALSE}
CDK4i_72hr_DSPG_vs_DEG <- CDK4i_72hr_DMSO_72hr_DSPG_DEG[, 2:3]
colnames(CDK4i_72hr_DSPG_vs_DEG) <- c('DSPG', 'DEG')

PRMT5i_72hr_DSPG_vs_DEG <- PRMT5i_72hr_DMSO_72hr_DSPG_DEG[, 2:3]
colnames(PRMT5i_72hr_DSPG_vs_DEG) <- c('DSPG', 'DEG')

CDK4i_6d_DSPG_vs_DEG <- CDK4i_6d_DMSO_72hr_DSPG_DEG[, 2:3]
colnames(CDK4i_6d_DSPG_vs_DEG) <- c('DSPG', 'DEG')
```

Visualize the differentially spliced genes and differentially expressed genes under CDK4i_72hr, PRMT5i_72hr and CDK4i_6d with venn diagrams.
```{r,results = 'asis'}
par(mfrow = c(2, 2))
vennDiagram(ifelse(is.na(CDK4i_72hr_DSPG_vs_DEG), 0, 1), cex = 1.0, main =
'CDK4i_72hr_DSPG_vs_DEG')
vennDiagram(ifelse(is.na(PRMT5i_72hr_DSPG_vs_DEG), 0, 1), cex = 1.0, main =
'PRMT5i_72hr_DSPG_vs_DEG')
vennDiagram(ifelse(is.na(CDK4i_6d_DSPG_vs_DEG), 0, 1), cex = 1.0, main =
'CDK4i_6d_DSPG_vs_DEG')
```

# 3. Gene sets analysis
## 3.1. Set-up
The listEnrichrDbs() function shows all available gene sets collection for the analysis. The GO_Biological_Process_2018 and GO_Molecular_Function_2018 were selected for the analysis.
```{r select_gene_sets_database,results = 'asis'}
dbs <- listEnrichrDbs()
selected_dbs <-
c('GO_Biological_Process_2018', 'GO_Molecular_Function_2018')
```

## 3.2. Gene set analysis on differentially spliced genes
Gene sets analysis were performed on a list of differentially spliced genes. Below shows the gene sets analysis for CDK4i_72hr.
```{r gene_sets_analysis,message=FALSE,results = 'asis'}
CDK4i_72hr_DMSO_72hr_spliced_gene <-
unique(
rbind(
CDK4i_72hr_DMSO_72hr_A3SS_gene,
CDK4i_72hr_DMSO_72hr_A5SS_gene,
CDK4i_72hr_DMSO_72hr_MXE_gene,
CDK4i_72hr_DMSO_72hr_RI_gene,
CDK4i_72hr_DMSO_72hr_SE_gene
)[, 1]
)

CDK4i_72hr_DMSO_72hr_Spliced_gene_enriched <-
enrichr(as.character(CDK4i_72hr_DMSO_72hr_spliced_gene),
selected_dbs)

```
```{r echo=FALSE,message=FALSE}
PRMT5i_72hr_DMSO_72hr_spliced_gene <-
unique(
rbind(
PRMT5i_72hr_DMSO_72hr_A3SS_gene,
PRMT5i_72hr_DMSO_72hr_A5SS_gene,
PRMT5i_72hr_DMSO_72hr_MXE_gene,
PRMT5i_72hr_DMSO_72hr_RI_gene,
PRMT5i_72hr_DMSO_72hr_SE_gene
)[, 1]
)

CDK4i_6d_DMSO_72hr_spliced_gene <-
unique(
rbind(
CDK4i_6d_DMSO_72hr_A3SS_gene,
CDK4i_6d_DMSO_72hr_A5SS_gene,
CDK4i_6d_DMSO_72hr_MXE_gene,
CDK4i_6d_DMSO_72hr_RI_gene,
CDK4i_6d_DMSO_72hr_SE_gene
)[, 1]
)

PRMT5i_72hr_DMSO_72hr_Spliced_gene_enriched <-
enrichr(as.character(PRMT5i_72hr_DMSO_72hr_spliced_gene),
selected_dbs)
CDK4i_6d_DMSO_72hr_Spliced_gene_enriched <-
enrichr(as.character(CDK4i_6d_DMSO_72hr_spliced_gene), selected_dbs)
```

For GO_Biological_Process_2018, filter significant gene sets by p-value < 0.05 and sort the significant gene sets by p-value. 
```{r GO_Biological_Process_2018_sig_gene_sets,results = 'asis'}
sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018 <-
CDK4i_72hr_DMSO_72hr_Spliced_gene_enriched$GO_Biological_Process_2018[CDK4i_72hr_DMSO_72hr_Spliced_gene_enriched$GO_Biological_Process_2018$Adjusted.P.value < 0.05,]

sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018 <-
sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018[order(sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018[, "Adjusted.P.value"], decreasing = FALSE),]

```
```{r echo=FALSE}
sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018 <-
PRMT5i_72hr_DMSO_72hr_Spliced_gene_enriched$GO_Biological_Process_2018[PRMT5i_72hr_DMSO_72hr_Spliced_gene_enriched$GO_Biological_Process_2018$Adjusted.P.value < 0.05,]

sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018 <-
sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018[order(sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018[, "Adjusted.P.value"], decreasing = FALSE),]

sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Biological_Process_2018 <-
CDK4i_6d_DMSO_72hr_Spliced_gene_enriched$GO_Biological_Process_2018[CDK4i_6d_DMSO_72hr_Spliced_gene_enriched$GO_Biological_Process_2018$Adjusted.P.value < 0.05,]

sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Biological_Process_2018 <-
sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Biological_Process_2018[order(sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Biological_Process_2018[, "Adjusted.P.value"], decreasing = FALSE),]
```

For GO_Molecular_Function_2018, filter significant gene sets by p-value < 0.05 and sort the significant gene sets by p-value.
```{r GO_Molecular_Function_2018_sig_gene_sets,results = 'asis'}
sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018 <-
CDK4i_72hr_DMSO_72hr_Spliced_gene_enriched$GO_Molecular_Function_2018[CDK4i_72hr_DMSO_72hr_Spliced_gene_enriched$GO_Molecular_Function_2018$Adjusted.P.value < 0.05, ]

sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018 <-
sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018[order(sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018[, "Adjusted.P.value"], decreasing = FALSE), ]

```
```{r echo=FALSE}
sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018 <-
PRMT5i_72hr_DMSO_72hr_Spliced_gene_enriched$GO_Molecular_Function_2018[PRMT5i_72hr_DMSO_72hr_Spliced_gene_enriched$GO_Molecular_Function_2018$Adjusted.P.value <  0.05,]

sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018 <-
sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018[order(sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018[, "Adjusted.P.value"], decreasing = FALSE),]

sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Molecular_Function_2018 <-
CDK4i_6d_DMSO_72hr_Spliced_gene_enriched$GO_Molecular_Function_2018[CDK4i_6d_DMSO_72hr_Spliced_gene_enriched$GO_Molecular_Function_2018$Adjusted.P.value < 0.05,]

sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Molecular_Function_2018 <-
sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Molecular_Function_2018[order(sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Molecular_Function_2018[, "Adjusted.P.value"], decreasing = FALSE),]
```

The results are stored in two files. The GO_Biological_Process_2018_gene_set_DSPG.xlsx contains significant GO_Biological_Process_2018 gene sets, while the GO_Molecular_Function_2018_gene_set_DSPG.xlsx contains significant GO_Molecular_Function_2018 gene sets.
```{r write_gene_set_analysis_result,echo=FALSE}
GO_Biological_Process_2018_DSPG <-
list(
'CDK4i_72hr_DMSO_72hr' = sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018,
'PRMT5i_72hr_DMSO_72hr' = sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018,
'CDK4i_6d_DMSO_72hr' = sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Biological_Process_2018
)
GO_Molecular_Function_2018_DSPG <-
list(
'CDK4i_72hr_DMSO_72hr' = sig_CDK4i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018,
'PRMT5i_72hr_DMSO_72hr' = sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018,
'CDK4i_6d_DMSO_72hr' = sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Molecular_Function_2018
)

write.xlsx(
GO_Biological_Process_2018_DSPG,
"Analysis/Gene_sets_analysis/DSPG/GO_Biological_Process_2018_gene_set_DSPG.xlsx"
)
write.xlsx(
GO_Molecular_Function_2018_DSPG,
"Analysis/Gene_sets_analysis/DSPG/GO_Molecular_Function_2018_gene_set_DSPG.xlsx"
)

```

Compare the gene sets between CDK4i_6d and PRMT5i_72hr.
```{r compare_splicing_gene_sets,results = 'asis'}
com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Biological_Process_2018 <-
merge(
sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018[, c(1, 3)],
sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Biological_Process_2018[, c(1, 3)],
by = 'Term',
all = TRUE
)
colnames(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Biological_Process_2018)[2:3] <-
c('PRMT5i', 'CDK4i_6d')
com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Molecular_Function_2018 <-
merge(
sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018[, c(1, 3)],
sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Molecular_Function_2018[, c(1, 3)],
by = 'Term',
all = TRUE
)
colnames(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Molecular_Function_2018)[2:3] <-
c('PRMT5i', 'CDK4i_6d')
```

Visualise the gene sets between CDK4i_6d and PRMT5i_72hr with venn diagrams.
```{r plot_splicig_gene_sets_between_CDK4i_&_PRMT5i,results = 'asis' }
par(mfrow = c(1, 2))
vennDiagram(ifelse(
is.na(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Biological_Process_2018[, 2:3]),
0,
1
), cex = 1.0, main = 'GO_Biological_Process_2018_DSPG')
vennDiagram(ifelse(
is.na(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Molecular_Function_2018[, 2:3]),
0,
1
), cex = 1.0, main = 'GO_Molecular_Function_2018_DSPG')
```

The compared results are stored in two files. The GO_Biological_Process_2018_gene_set_DSPG_compared.xlsx and the GO_Molecular_Function_2018_gene_set_DSPG_compared.xlsx contains gene sets that are commonly detected in both inhibitions and distinct in both inhibitions.
```{r write_distinct_&_common_splicing_gene_sets,echo=FALSE}
GO_Biological_Process_2018_DSPG_compared <-
list(
'common' = na.omit(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Biological_Process_2018),
'Distinct under PRMT5i' = com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Biological_Process_2018[is.na(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Biological_Process_2018$CDK4i_6d), ],
'Distinct under CDK4i' = com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Biological_Process_2018[is.na(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Biological_Process_2018$PRMT5i_6d), ]
)

GO_Molecular_Function_2018_DSPG_compared <-
list(
'common' = na.omit(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Molecular_Function_2018),
'Distinct under PRMT5i' = com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Molecular_Function_2018[is.na(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Molecular_Function_2018$CDK4i_6d), ],
'Distinct under CDK4i' = com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Molecular_Function_2018[is.na(com_PRMT5i_72hr_CDK4i_6d_DSPG_GO_Molecular_Function_2018$PRMT5i_6d), ]
)

write.xlsx(
GO_Biological_Process_2018_DSPG_compared,
"Analysis/Gene_sets_analysis/DSPG/GO_Biological_Process_2018_gene_set_DSPG_compared.xlsx"
)
write.xlsx(
GO_Molecular_Function_2018_DSPG_compared,
"Analysis/Gene_sets_analysis/DSPG/GO_Molecular_Function_2018_gene_set_DSPG_compared.xlsx"
)
```

## 3.3. Gene sets analysis on differentially expressed genes
Gene sets analysis were performed on a list of differentially expressed genes. Below shows the gene sets analysis for CDK4i_72hr.
```{r gene_set_analysis_on_DEGs,message=FALSE,results = 'asis'}
CDK4i_72hr_DMSO_72hr_DE_gene_enriched <-
enrichr(CDK4i_72hr_DMSO_72hr_DEG[, 1], selected_dbs)
```
```{r echo=FALSE}
PRMT5i_72hr_DMSO_72hr_DE_gene_enriched <-
enrichr(PRMT5i_72hr_DMSO_72hr_DEG[, 1], selected_dbs)
CDK4i_6d_DMSO_72hr_DE_gene_enriched <-
enrichr(CDK4i_6d_DMSO_72hr_DEG[, 1], selected_dbs)
```

Like the GO_Biological_Process_2018  differentially spliced gene sets, significant gene sets are filtered by p-value < 0.05 and sorted by p-value.
```{r GO_Biological_Process_2018_sig_gene_sets_DGEs,results = 'asis'}
sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018 <-
CDK4i_72hr_DMSO_72hr_DE_gene_enriched$GO_Biological_Process_2018[CDK4i_72hr_DMSO_72hr_DE_gene_enriched$GO_Biological_Process_2018$Adjusted.P.value < 0.05, ]

sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018 <-
sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018[order(sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018[, "Adjusted.P.value"], decreasing = FALSE), ]

```
```{r echo=FALSE}
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018 <-
PRMT5i_72hr_DMSO_72hr_DE_gene_enriched$GO_Biological_Process_2018[PRMT5i_72hr_DMSO_72hr_DE_gene_enriched$GO_Biological_Process_2018$Adjusted.P.value < 0.05, ]

sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018 <-
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018[order(sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018[, "Adjusted.P.value"], decreasing = FALSE), ]

sig_CDK4i_6d_DMSO_72hr_DEG_GO_Biological_Process_2018 <-
CDK4i_6d_DMSO_72hr_DE_gene_enriched$GO_Biological_Process_2018[CDK4i_6d_DMSO_72hr_DE_gene_enriched$GO_Biological_Process_2018$Adjusted.P.value < 0.05, ]

sig_CDK4i_6d_DMSO_72hr_DEG_GO_Biological_Process_2018 <-
sig_CDK4i_6d_DMSO_72hr_DEG_GO_Biological_Process_2018[order(sig_CDK4i_6d_DMSO_72hr_DEG_GO_Biological_Process_2018[, "Adjusted.P.value"], decreasing = FALSE), ]
```

Like the GO_Molecular_Function_2018 differentially spliced gene sets, significant gene sets are filtered by p-value < 0.05 and sorted by p-value.
```{r GO_Molecular_Function_2018_sig_gene_sets_DGEs,results = 'asis'}
sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018 <-
CDK4i_72hr_DMSO_72hr_DE_gene_enriched$GO_Molecular_Function_2018[CDK4i_72hr_DMSO_72hr_DE_gene_enriched$GO_Molecular_Function_2018$Adjusted.P.value <
0.05, ]
sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018 <-
sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018[order(sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018[, "Adjusted.P.value"],
decreasing = FALSE), ]

```
```{r echo=FALSE}
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018 <-
PRMT5i_72hr_DMSO_72hr_DE_gene_enriched$GO_Molecular_Function_2018[PRMT5i_72hr_DMSO_72hr_DE_gene_enriched$GO_Molecular_Function_2018$Adjusted.P.value <
0.05, ]
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018 <-
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018[order(sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018[, "Adjusted.P.value"],
decreasing = FALSE), ]

sig_CDK4i_6d_DMSO_72hr_DEG_GO_Molecular_Function_2018 <-
CDK4i_6d_DMSO_72hr_DE_gene_enriched$GO_Molecular_Function_2018[CDK4i_6d_DMSO_72hr_DE_gene_enriched$GO_Molecular_Function_2018$Adjusted.P.value <
0.05, ]
sig_CDK4i_6d_DMSO_72hr_DEG_GO_Molecular_Function_2018 <-
sig_CDK4i_6d_DMSO_72hr_DEG_GO_Molecular_Function_2018[order(sig_CDK4i_6d_DMSO_72hr_DEG_GO_Molecular_Function_2018[, "Adjusted.P.value"],
decreasing = FALSE), ]
```

The results are stored in two files. The GO_Biological_Process_2018_gene_set_DEG.xlsx contains significant GO_Biological_Process_2018 gene sets, while the GO_Molecular_Function_2018_gene_set_DEG.xlsx contains significant GO_Molecular_Function_2018 gene sets.
```{r write_gene_set_analysis_DEG_result,echo=FALSE}
GO_Biological_Process_2018_DEG <-
list(
'CDK4i_72hr_DMSO_72hr' = sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018,
'PRMT5i_72hr_DMSO_72hr' = sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018,
'CDK4i_6d_DMSO_72hr' = sig_CDK4i_6d_DMSO_72hr_DEG_GO_Biological_Process_2018
)
GO_Molecular_Function_2018_DEG <-
list(
'CDK4i_72hr_DMSO_72hr' = sig_CDK4i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018,
'PRMT5i_72hr_DMSO_72hr' = sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018,
'CDK4i_6d_DMSO_72hr' = sig_CDK4i_6d_DMSO_72hr_DEG_GO_Molecular_Function_2018
)

write.xlsx(
GO_Biological_Process_2018_DSPG,
"Analysis/Gene_sets_analysis/DEG/GO_Biological_Process_2018_gene_set_DEG.xlsx"
)
write.xlsx(
GO_Molecular_Function_2018_DSPG,
"Analysis/Gene_sets_analysis/DEG/GO_Molecular_Function_2018_gene_set_DEG.xlsx"
)

```

Compare the gene sets between CDK4i_6d and PRMT5i_72hr.
```{r compare_DE_gene_sets,results = 'asis'}
com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Biological_Process_2018 <-
merge(
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018[, c(1, 3)],
sig_CDK4i_6d_DMSO_72hr_DEG_GO_Biological_Process_2018[, c(1, 3)],
by = 'Term',
all = TRUE
)
colnames(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Biological_Process_2018)[2:3] <-
c('PRMT5i', 'CDK4i_6d')

com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Molecular_Function_2018 <-
merge(
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018[, c(1, 3)],
sig_CDK4i_6d_DMSO_72hr_DEG_GO_Molecular_Function_2018[, c(1, 3)],
by = 'Term',
all = TRUE
)
colnames(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Molecular_Function_2018)[2:3] <-
c('PRMT5i', 'CDK4i_6d')
```

Visualise the gene sets between CDK4i_6d and PRMT5i_72hr with venn diagrams.
```{r plot_DE_gene_sets_between_CDK4i_&_PRMT5i,results = 'asis' }
par(mfrow = c(1, 2))
vennDiagram(ifelse(
is.na(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Biological_Process_2018[, 2:3]),
0,
1
), cex = 1.0, main = 'GO_Biological_Process_2018_DEG')
vennDiagram(ifelse(
is.na(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Molecular_Function_2018[, 2:3]),
0,
1
), cex = 1.0, main = 'GO_Molecular_Function_2018_DEG')
```

The compared results are stored in two files. The GO_Biological_Process_2018_gene_set_DEG_compared.xlsx and the GO_Molecular_Function_2018_gene_set_DEG_compared.xlsx contains gene sets that are commonly detected in both inhibitions and distinct in both inhibitions.
```{r write_distinct_&_common_DE_gene_sets, echo=FALSE}
GO_Biological_Process_2018_DEG_compared <-
list(
'common' = na.omit(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Biological_Process_2018),
'Distinct under PRMT5i' = com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Biological_Process_2018[is.na(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Biological_Process_2018$CDK4i_6d),],
'Distinct under CDK4i' = com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Biological_Process_2018[is.na(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Biological_Process_2018$PRMT5i_6d),]
)

GO_Molecular_Function_2018_DEG_compared <-
list(
'common' = na.omit(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Molecular_Function_2018),
'Distinct under PRMT5i' = com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Molecular_Function_2018[is.na(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Molecular_Function_2018$CDK4i_6d),],
'Distinct under CDK4i' = com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Molecular_Function_2018[is.na(com_PRMT5i_72hr_CDK4i_6d_DEG_GO_Molecular_Function_2018$PRMT5i_6d),]
)

write.xlsx(
GO_Biological_Process_2018_DEG_compared,
"Analysis/Gene_sets_analysis/DEG/GO_Biological_Process_2018_gene_set_DEG_compared.xlsx"
)
write.xlsx(
GO_Molecular_Function_2018_DEG_compared,
"Analysis/Gene_sets_analysis/DEG/GO_Molecular_Function_2018_gene_set_DEG_compared.xlsx"
)
```

## 3.4. Differentially spliced gene sets vs differentially expressed gene sets
Compare the gene sets detected with diffentially spliced genes and differentially expressed genes. Below shows the comparison on the GO_Biological_Process_2018 gene sets collection.
```{r GO_Biological_Process_2018_DGE_DSPG_gene_sets_compare,results = 'asis'}
PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared <-
merge(
sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Biological_Process_2018[, c('Term', 'Adjusted.P.value')],
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Biological_Process_2018[, c('Term', 'Adjusted.P.value')],
by = 'Term',
all = TRUE
)
colnames(PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared) <-
c('GO_terms', 'DSPG', 'DGE')

CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared <-
merge(
sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Biological_Process_2018[, c('Term', 'Adjusted.P.value')],
sig_CDK4i_6d_DMSO_72hr_DEG_GO_Biological_Process_2018[, c('Term', 'Adjusted.P.value')],
by = 'Term',
all = TRUE
)
colnames(CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared) <-
c('GO_terms', 'DSPG', 'DGE')
```
```{r GO_Molecular_Function_2018_DGE_DSPG_gene_sets_compare, echo=FALSE}
PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared <-
merge(
sig_PRMT5i_72hr_DMSO_72hr_DSPG_GO_Molecular_Function_2018[, c('Term', 'Adjusted.P.value')],
sig_PRMT5i_72hr_DMSO_72hr_DEG_GO_Molecular_Function_2018[, c('Term', 'Adjusted.P.value')],
by = 'Term',
all = TRUE
)
colnames(PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared) <-
c('GO_terms', 'DSPG', 'DGE')

CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared <-
merge(
sig_CDK4i_6d_DMSO_72hr_DSPG_GO_Molecular_Function_2018[, c('Term', 'Adjusted.P.value')],
sig_CDK4i_6d_DMSO_72hr_DEG_GO_Molecular_Function_2018[, c('Term', 'Adjusted.P.value')],
by = 'Term',
all = TRUE
)
colnames(CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared) <-
c('GO_terms', 'DSPG', 'DGE')
```

Visualise the gene sets detected between differentially spliced genes and differentially expressed genes with venn diagrams.
```{r plot_gene_sets_between_DSPG_&_DEG,results = 'asis'}
par(mfrow = c(2, 2))
vennDiagram(ifelse(
is.na(
PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared[, 2:3]
),
0,
1
), cex = 1.0, main = 'GO_Biological_Process_2018_PRMT5i')
vennDiagram(ifelse(
is.na(
CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared[, 2:3]
),
0,
1
), cex = 1.0, main = 'GO_Biological_Process_2018_CDK4i_6d')
vennDiagram(ifelse(
is.na(
PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared[, 2:3]
),
0,
1
), cex = 1.0, main = 'GO_Molecular_Function_2018_PRMT5i')
vennDiagram(ifelse(
is.na(
CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared[, 2:3]
),
0,
1
), cex = 1.0, main = 'GO_Molecular_Function_2018_CDK4i_6d')
```

Four files are generated. The GO_Biological_Process_2018_PRMT5i_DSPG_DEG_compared.xlsx, GO_Biological_Process_2018_CDK4i_DSPG_DEG_compared.xlsx,
GO_Molecular_Function_2018_PRMT5i_DSPG_DEG_compared.xlsx and GO_Molecular_Function_2018_CDK4i_DSPG_DEG_compared.xlsx contains gene sets that are commonly and distinctly detected with both differentially spliced genes and differentially expressed genes.
```{r write_distinct_&_common_splicing_&_DE_gene_sets,echo=FALSE}
GO_Biological_Process_2018_PRMT5i_DSPG_DEG_compared <-
list(
'common' = na.omit(
PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared
),
'Distinct under DSPG' = PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared[is.na(PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared$DGE), ],
'Distinct under DEG' = PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared[is.na(PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared$DSPG), ]
)

GO_Biological_Process_2018_CDK4i_DSPG_DEG_compared <-
list(
'common' = na.omit(
CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared
),
'Distinct under DSPG' = CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared[is.na(CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared$DGE), ],
'Distinct under DEG' = CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared[is.na(CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Biological_Process_2018_compared$DSPG), ]
)

GO_Molecular_Function_2018_PRMT5i_DSPG_DEG_compared <-
list(
'common' = na.omit(
PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared
),
'Distinct under DSPG' = PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared[is.na(PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared$DGE), ],
'Distinct under DEG' = PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared[is.na(PRMT5i_72hr_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared$DSPG), ]
)

GO_Molecular_Function_2018_CDK4i_DSPG_DEG_compared <-
list(
'common' = na.omit(
CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared
),
'Distinct under DSPG' = CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared[is.na(CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared$DGE), ],
'Distinct under DEG' = CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared[is.na(CDK4i_6d_DMSO_72hr_DSPG_DGE_GO_Molecular_Function_2018_compared$DSPG), ]
)

write.xlsx(
GO_Biological_Process_2018_PRMT5i_DSPG_DEG_compared,
"Analysis/Gene_sets_analysis/DSPG vs DEG/GO_Biological_Process_2018_PRMT5i_DSPG_DEG_compared.xlsx"
)
write.xlsx(
GO_Biological_Process_2018_CDK4i_DSPG_DEG_compared,
"Analysis/Gene_sets_analysis/DSPG vs DEG/GO_Biological_Process_2018_CDK4i_DSPG_DEG_compared.xlsx"
)
write.xlsx(
GO_Molecular_Function_2018_PRMT5i_DSPG_DEG_compared,
"Analysis/Gene_sets_analysis/DSPG vs DEG/GO_Molecular_Function_2018_PRMT5i_DSPG_DEG_compared.xlsx"
)
write.xlsx(
GO_Molecular_Function_2018_CDK4i_DSPG_DEG_compared,
"Analysis/Gene_sets_analysis/DSPG vs DEG/GO_Molecular_Function_2018_CDK4i_DSPG_DEG_compared.xlsx"
)
```

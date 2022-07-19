# library -----------------------------------------------------------------
library('UpSetR')
library('ComplexHeatmap')

# Set working directory ---------------------------------------------------
setwd("/Users/lokhangc/Desktop/Master Research/Data/rMATS/RL149_AL08_hg19_VaryRL_novelSS/")

# read rMATS results ------------------------------------------------------
CDK4i_72hr_A3SS <- read.delim("CDK4i_72hr_DMSO_72hr/A3SS.MATS.JC.txt")
CDK4i_72hr_A5SS <- read.delim("CDK4i_72hr_DMSO_72hr/A5SS.MATS.JC.txt")
CDK4i_72hr_MXE  <- read.delim("CDK4i_72hr_DMSO_72hr/MXE.MATS.JC.txt")
CDK4i_72hr_RI <- read.delim("CDK4i_72hr_DMSO_72hr/RI.MATS.JC.txt")
CDK4i_72hr_SE <- read.delim("CDK4i_72hr_DMSO_72hr/SE.MATS.JC.txt")

CDK4i_6day_A3SS <- read.delim("CDK4i_6d_DMSO_72hr/A3SS.MATS.JC.txt")
CDK4i_6day_A5SS <- read.delim("CDK4i_6d_DMSO_72hr/A5SS.MATS.JC.txt")
CDK4i_6day_MXE  <- read.delim("CDK4i_6d_DMSO_72hr/MXE.MATS.JC.txt")
CDK4i_6day_RI <- read.delim("CDK4i_6d_DMSO_72hr/RI.MATS.JC.txt")
CDK4i_6day_SE <- read.delim("CDK4i_6d_DMSO_72hr/SE.MATS.JC.txt")

PRMT5i_72hr_A3SS <- read.delim("PRMT5i_72hr_DMSO_72hr/A3SS.MATS.JC.txt")
PRMT5i_72hr_A5SS <- read.delim("PRMT5i_72hr_DMSO_72hr/A5SS.MATS.JC.txt")
PRMT5i_72hr_MXE  <- read.delim("PRMT5i_72hr_DMSO_72hr/MXE.MATS.JC.txt")
PRMT5i_72hr_RI <- read.delim("PRMT5i_72hr_DMSO_72hr/RI.MATS.JC.txt")
PRMT5i_72hr_SE <- read.delim("PRMT5i_72hr_DMSO_72hr/SE.MATS.JC.txt")


# filter by FDR 0.05 ------------------------------------------------------
CDK4i_72hr_A3SS_FDR <- CDK4i_72hr_A3SS[CDK4i_72hr_A3SS$FDR<0.05,]
CDK4i_72hr_A5SS_FDR <- CDK4i_72hr_A5SS[CDK4i_72hr_A5SS$FDR<0.05,]
CDK4i_72hr_MXE_FDR <- CDK4i_72hr_MXE[CDK4i_72hr_MXE$FDR<0.05,]
CDK4i_72hr_RI_FDR <- CDK4i_72hr_RI[CDK4i_72hr_RI$FDR<0.05,]
CDK4i_72hr_SE_FDR <- CDK4i_72hr_SE[CDK4i_72hr_SE$FDR<0.05,]

CDK4i_6day_A3SS_FDR <- CDK4i_6day_A3SS[CDK4i_6day_A3SS$FDR<0.05,]
CDK4i_6day_A5SS_FDR <- CDK4i_6day_A5SS[CDK4i_6day_A5SS$FDR<0.05,]
CDK4i_6day_MXE_FDR <- CDK4i_6day_MXE[CDK4i_6day_MXE$FDR<0.05,]
CDK4i_6day_RI_FDR <- CDK4i_6day_RI[CDK4i_6day_RI$FDR<0.05,]
CDK4i_6day_SE_FDR <- CDK4i_6day_SE[CDK4i_6day_SE$FDR<0.05,]

PRMT5i_72hr_A3SS_FDR <- PRMT5i_72hr_A3SS[PRMT5i_72hr_A3SS$FDR<0.05,]
PRMT5i_72hr_A5SS_FDR <- PRMT5i_72hr_A5SS[PRMT5i_72hr_A5SS$FDR<0.05,]
PRMT5i_72hr_MXE_FDR <- PRMT5i_72hr_MXE[PRMT5i_72hr_MXE$FDR<0.05,]
PRMT5i_72hr_RI_FDR <- PRMT5i_72hr_RI[PRMT5i_72hr_RI$FDR<0.05,]
PRMT5i_72hr_SE_FDR <- PRMT5i_72hr_SE[PRMT5i_72hr_SE$FDR<0.05,]


# up-regulated splice events ----------------------------------------------
CDK4i_72hr_A3SS_FDR_up <- CDK4i_72hr_A3SS_FDR[CDK4i_72hr_A3SS_FDR$IncLevelDifference>0,]
CDK4i_72hr_A5SS_FDR_up <- CDK4i_72hr_A5SS_FDR[CDK4i_72hr_A5SS_FDR$IncLevelDifference>0,]
CDK4i_72hr_MXE_FDR_up <- CDK4i_72hr_MXE_FDR[CDK4i_72hr_MXE_FDR$IncLevelDifference>0,]
CDK4i_72hr_RI_FDR_up <- CDK4i_72hr_RI_FDR[CDK4i_72hr_RI_FDR$IncLevelDifference>0,]
CDK4i_72hr_SE_FDR_up <- CDK4i_72hr_SE_FDR[CDK4i_72hr_SE_FDR$IncLevelDifference>0,]

CDK4i_6day_A3SS_FDR_up <- CDK4i_6day_A3SS_FDR[CDK4i_6day_A3SS_FDR$IncLevelDifference>0,]
CDK4i_6day_A5SS_FDR_up <- CDK4i_6day_A5SS_FDR[CDK4i_6day_A5SS_FDR$IncLevelDifference>0,]
CDK4i_6day_MXE_FDR_up <- CDK4i_6day_MXE_FDR[CDK4i_6day_MXE_FDR$IncLevelDifference>0,]
CDK4i_6day_RI_FDR_up <- CDK4i_6day_RI_FDR[CDK4i_6day_RI_FDR$IncLevelDifference>0,]
CDK4i_6day_SE_FDR_up <- CDK4i_6day_SE_FDR[CDK4i_6day_SE_FDR$IncLevelDifference>0,]

PRMT5i_72hr_A3SS_FDR_up <- PRMT5i_72hr_A3SS_FDR[PRMT5i_72hr_A3SS_FDR$IncLevelDifference>0,]
PRMT5i_72hr_A5SS_FDR_up <- PRMT5i_72hr_A5SS_FDR[PRMT5i_72hr_A5SS_FDR$IncLevelDifference>0,]
PRMT5i_72hr_MXE_FDR_up <- PRMT5i_72hr_MXE_FDR[PRMT5i_72hr_MXE_FDR$IncLevelDifference>0,]
PRMT5i_72hr_RI_FDR_up <- PRMT5i_72hr_RI_FDR[PRMT5i_72hr_RI_FDR$IncLevelDifference>0,]
PRMT5i_72hr_SE_FDR_up <- PRMT5i_72hr_SE_FDR[PRMT5i_72hr_SE_FDR$IncLevelDifference>0,]


# Down-regulated splice events --------------------------------------------
CDK4i_72hr_A3SS_FDR_dn <- CDK4i_72hr_A3SS_FDR[CDK4i_72hr_A3SS_FDR$IncLevelDifference<0,]
CDK4i_72hr_A5SS_FDR_dn <- CDK4i_72hr_A5SS_FDR[CDK4i_72hr_A5SS_FDR$IncLevelDifference<0,]
CDK4i_72hr_MXE_FDR_dn <- CDK4i_72hr_MXE_FDR[CDK4i_72hr_MXE_FDR$IncLevelDifference<0,]
CDK4i_72hr_RI_FDR_dn <- CDK4i_72hr_RI_FDR[CDK4i_72hr_RI_FDR$IncLevelDifference<0,]
CDK4i_72hr_SE_FDR_dn <- CDK4i_72hr_SE_FDR[CDK4i_72hr_SE_FDR$IncLevelDifference<0,]

CDK4i_6day_A3SS_FDR_dn <- CDK4i_6day_A3SS_FDR[CDK4i_6day_A3SS_FDR$IncLevelDifference<0,]
CDK4i_6day_A5SS_FDR_dn <- CDK4i_6day_A5SS_FDR[CDK4i_6day_A5SS_FDR$IncLevelDifference<0,]
CDK4i_6day_MXE_FDR_dn <- CDK4i_6day_MXE_FDR[CDK4i_6day_MXE_FDR$IncLevelDifference<0,]
CDK4i_6day_RI_FDR_dn <- CDK4i_6day_RI_FDR[CDK4i_6day_RI_FDR$IncLevelDifference<0,]
CDK4i_6day_SE_FDR_dn <- CDK4i_6day_SE_FDR[CDK4i_6day_SE_FDR$IncLevelDifference<0,]

PRMT5i_72hr_A3SS_FDR_dn <- PRMT5i_72hr_A3SS_FDR[PRMT5i_72hr_A3SS_FDR$IncLevelDifference<0,]
PRMT5i_72hr_A5SS_FDR_dn <- PRMT5i_72hr_A5SS_FDR[PRMT5i_72hr_A5SS_FDR$IncLevelDifference<0,]
PRMT5i_72hr_MXE_FDR_dn <- PRMT5i_72hr_MXE_FDR[PRMT5i_72hr_MXE_FDR$IncLevelDifference<0,]
PRMT5i_72hr_RI_FDR_dn <- PRMT5i_72hr_RI_FDR[PRMT5i_72hr_RI_FDR$IncLevelDifference<0,]
PRMT5i_72hr_SE_FDR_dn <- PRMT5i_72hr_SE_FDR[PRMT5i_72hr_SE_FDR$IncLevelDifference<0,]


# Up-regulated events pre-process -----------------------------------------
## A3SS --------------------------------------------------------------------
CDK4i_72hr_A3SS_FDR_up_combine <- CDK4i_72hr_A3SS_FDR_up
CDK4i_72hr_A3SS_FDR_up_combine$combine <- paste(CDK4i_72hr_A3SS_FDR_up$geneSymbol,
                                                CDK4i_72hr_A3SS_FDR_up$longExonStart_0base,
                                                CDK4i_72hr_A3SS_FDR_up$longExonEnd,
                                                CDK4i_72hr_A3SS_FDR_up$shortES,
                                                CDK4i_72hr_A3SS_FDR_up$shortEE,
                                                CDK4i_72hr_A3SS_FDR_up$flankingES,
                                                CDK4i_72hr_A3SS_FDR_up$flankingEE,
                                                sep = '_')

CDK4i_6day_A3SS_FDR_up_combine <- CDK4i_6day_A3SS_FDR_up
CDK4i_6day_A3SS_FDR_up_combine$combine <- paste(CDK4i_6day_A3SS_FDR_up$geneSymbol,
                                                CDK4i_6day_A3SS_FDR_up$longExonStart_0base,
                                                CDK4i_6day_A3SS_FDR_up$longExonEnd,
                                                CDK4i_6day_A3SS_FDR_up$shortES,
                                                CDK4i_6day_A3SS_FDR_up$shortEE,
                                                CDK4i_6day_A3SS_FDR_up$flankingES,
                                                CDK4i_6day_A3SS_FDR_up$flankingEE,
                                                sep = '_')

PRMT5i_72hr_A3SS_FDR_up_combine <- PRMT5i_72hr_A3SS_FDR_up
PRMT5i_72hr_A3SS_FDR_up_combine$combine <- paste(PRMT5i_72hr_A3SS_FDR_up$geneSymbol,
                                                PRMT5i_72hr_A3SS_FDR_up$longExonStart_0base,
                                                PRMT5i_72hr_A3SS_FDR_up$longExonEnd,
                                                PRMT5i_72hr_A3SS_FDR_up$shortES,
                                                PRMT5i_72hr_A3SS_FDR_up$shortEE,
                                                PRMT5i_72hr_A3SS_FDR_up$flankingES,
                                                PRMT5i_72hr_A3SS_FDR_up$flankingEE,
                                                sep = '_')


## A5SS --------------------------------------------------------------------
CDK4i_72hr_A5SS_FDR_up_combine <- CDK4i_72hr_A5SS_FDR_up
CDK4i_72hr_A5SS_FDR_up_combine$combine <- paste(CDK4i_72hr_A5SS_FDR_up$geneSymbol,
                                                CDK4i_72hr_A5SS_FDR_up$longExonStart_0base,
                                                CDK4i_72hr_A5SS_FDR_up$longExonEnd,
                                                CDK4i_72hr_A5SS_FDR_up$shortES,
                                                CDK4i_72hr_A5SS_FDR_up$shortEE,
                                                CDK4i_72hr_A5SS_FDR_up$flankingES,
                                                CDK4i_72hr_A5SS_FDR_up$flankingEE,
                                                sep = '_')

CDK4i_6day_A5SS_FDR_up_combine <- CDK4i_6day_A5SS_FDR_up
CDK4i_6day_A5SS_FDR_up_combine$combine <- paste(CDK4i_6day_A5SS_FDR_up$geneSymbol,
                                                CDK4i_6day_A5SS_FDR_up$longExonStart_0base,
                                                CDK4i_6day_A5SS_FDR_up$longExonEnd,
                                                CDK4i_6day_A5SS_FDR_up$shortES,
                                                CDK4i_6day_A5SS_FDR_up$shortEE,
                                                CDK4i_6day_A5SS_FDR_up$flankingES,
                                                CDK4i_6day_A5SS_FDR_up$flankingEE,
                                                sep = '_')

PRMT5i_72hr_A5SS_FDR_up_combine <- PRMT5i_72hr_A5SS_FDR_up
PRMT5i_72hr_A5SS_FDR_up_combine$combine <- paste(PRMT5i_72hr_A5SS_FDR_up$geneSymbol,
                                                 PRMT5i_72hr_A5SS_FDR_up$longExonStart_0base,
                                                 PRMT5i_72hr_A5SS_FDR_up$longExonEnd,
                                                 PRMT5i_72hr_A5SS_FDR_up$shortES,
                                                 PRMT5i_72hr_A5SS_FDR_up$shortEE,
                                                 PRMT5i_72hr_A5SS_FDR_up$flankingES,
                                                 PRMT5i_72hr_A5SS_FDR_up$flankingEE,
                                                 sep = '_')

## MXE --------------------------------------------------------------------
CDK4i_72hr_MXE_FDR_up_combine <- CDK4i_72hr_MXE_FDR_up
CDK4i_72hr_MXE_FDR_up_combine$combine <- paste(CDK4i_72hr_MXE_FDR_up$geneSymbol,
                                                CDK4i_72hr_MXE_FDR_up$X1stExonStart_0base,
                                                CDK4i_72hr_MXE_FDR_up$X1stExonEnd,
                                                CDK4i_72hr_MXE_FDR_up$X2ndExonStart_0base,
                                                CDK4i_72hr_MXE_FDR_up$X2ndExonEnd,
                                                CDK4i_72hr_MXE_FDR_up$upstreamES,
                                                CDK4i_72hr_MXE_FDR_up$upstreamEE,
                                               CDK4i_72hr_MXE_FDR_up$downstreamES,
                                               CDK4i_72hr_MXE_FDR_up$downstreamEE,
                                                sep = '_')

CDK4i_6day_MXE_FDR_up_combine <- CDK4i_6day_MXE_FDR_up
CDK4i_6day_MXE_FDR_up_combine$combine <- paste(CDK4i_6day_MXE_FDR_up$geneSymbol,
                                               CDK4i_6day_MXE_FDR_up$X1stExonStart_0base,
                                               CDK4i_6day_MXE_FDR_up$X1stExonEnd,
                                               CDK4i_6day_MXE_FDR_up$X2ndExonStart_0base,
                                               CDK4i_6day_MXE_FDR_up$X2ndExonEnd,
                                               CDK4i_6day_MXE_FDR_up$upstreamES,
                                               CDK4i_6day_MXE_FDR_up$upstreamEE,
                                               CDK4i_6day_MXE_FDR_up$downstreamES,
                                               CDK4i_6day_MXE_FDR_up$downstreamEE,
                                               sep = '_')

PRMT5i_72hr_MXE_FDR_up_combine <- PRMT5i_72hr_MXE_FDR_up
PRMT5i_72hr_MXE_FDR_up_combine$combine <- paste(PRMT5i_72hr_MXE_FDR_up$geneSymbol,
                                               PRMT5i_72hr_MXE_FDR_up$X1stExonStart_0base,
                                               PRMT5i_72hr_MXE_FDR_up$X1stExonEnd,
                                               PRMT5i_72hr_MXE_FDR_up$X2ndExonStart_0base,
                                               PRMT5i_72hr_MXE_FDR_up$X2ndExonEnd,
                                               PRMT5i_72hr_MXE_FDR_up$upstreamES,
                                               PRMT5i_72hr_MXE_FDR_up$upstreamEE,
                                               PRMT5i_72hr_MXE_FDR_up$downstreamES,
                                               PRMT5i_72hr_MXE_FDR_up$downstreamEE,
                                               sep = '_')


## RI --------------------------------------------------------------------
CDK4i_72hr_RI_FDR_up_combine <- CDK4i_72hr_RI_FDR_up
CDK4i_72hr_RI_FDR_up_combine$combine <- paste(CDK4i_72hr_RI_FDR_up$geneSymbol,
                                                CDK4i_72hr_RI_FDR_up$longExonStart_0base,
                                                CDK4i_72hr_RI_FDR_up$longExonEnd,
                                                CDK4i_72hr_RI_FDR_up$shortES,
                                                CDK4i_72hr_RI_FDR_up$shortEE,
                                                CDK4i_72hr_RI_FDR_up$flankingES,
                                                CDK4i_72hr_RI_FDR_up$flankingEE,
                                                sep = '_')

CDK4i_6day_RI_FDR_up_combine <- CDK4i_6day_RI_FDR_up
CDK4i_6day_RI_FDR_up_combine$combine <- paste(CDK4i_6day_RI_FDR_up$geneSymbol,
                                                CDK4i_6day_RI_FDR_up$longExonStart_0base,
                                                CDK4i_6day_RI_FDR_up$longExonEnd,
                                                CDK4i_6day_RI_FDR_up$shortES,
                                                CDK4i_6day_RI_FDR_up$shortEE,
                                                CDK4i_6day_RI_FDR_up$flankingES,
                                                CDK4i_6day_RI_FDR_up$flankingEE,
                                                sep = '_')

PRMT5i_72hr_RI_FDR_up_combine <- PRMT5i_72hr_RI_FDR_up
PRMT5i_72hr_RI_FDR_up_combine$combine <- paste(PRMT5i_72hr_RI_FDR_up$geneSymbol,
                                                 PRMT5i_72hr_RI_FDR_up$longExonStart_0base,
                                                 PRMT5i_72hr_RI_FDR_up$longExonEnd,
                                                 PRMT5i_72hr_RI_FDR_up$shortES,
                                                 PRMT5i_72hr_RI_FDR_up$shortEE,
                                                 PRMT5i_72hr_RI_FDR_up$flankingES,
                                                 PRMT5i_72hr_RI_FDR_up$flankingEE,
                                                 sep = '_')

## SE --------------------------------------------------------------------
CDK4i_72hr_SE_FDR_up_combine <- CDK4i_72hr_SE_FDR_up
CDK4i_72hr_SE_FDR_up_combine$combine <- paste(CDK4i_72hr_SE_FDR_up$geneSymbol,
                                              CDK4i_72hr_SE_FDR_up$longExonStart_0base,
                                              CDK4i_72hr_SE_FDR_up$longExonEnd,
                                              CDK4i_72hr_SE_FDR_up$shortES,
                                              CDK4i_72hr_SE_FDR_up$shortEE,
                                              CDK4i_72hr_SE_FDR_up$flankingES,
                                              CDK4i_72hr_SE_FDR_up$flankingEE,
                                              sep = '_')

CDK4i_6day_SE_FDR_up_combine <- CDK4i_6day_SE_FDR_up
CDK4i_6day_SE_FDR_up_combine$combine <- paste(CDK4i_6day_SE_FDR_up$geneSymbol,
                                              CDK4i_6day_SE_FDR_up$longExonStart_0base,
                                              CDK4i_6day_SE_FDR_up$longExonEnd,
                                              CDK4i_6day_SE_FDR_up$shortES,
                                              CDK4i_6day_SE_FDR_up$shortEE,
                                              CDK4i_6day_SE_FDR_up$flankingES,
                                              CDK4i_6day_SE_FDR_up$flankingEE,
                                              sep = '_')

PRMT5i_72hr_SE_FDR_up_combine <- PRMT5i_72hr_SE_FDR_up
PRMT5i_72hr_SE_FDR_up_combine$combine <- paste(PRMT5i_72hr_SE_FDR_up$geneSymbol,
                                               PRMT5i_72hr_SE_FDR_up$longExonStart_0base,
                                               PRMT5i_72hr_SE_FDR_up$longExonEnd,
                                               PRMT5i_72hr_SE_FDR_up$shortES,
                                               PRMT5i_72hr_SE_FDR_up$shortEE,
                                               PRMT5i_72hr_SE_FDR_up$flankingES,
                                               PRMT5i_72hr_SE_FDR_up$flankingEE,
                                               sep = '_')




# Down-regulated events pre-process -----------------------------------------
## A3SS --------------------------------------------------------------------
CDK4i_72hr_A3SS_FDR_dn_combine <- CDK4i_72hr_A3SS_FDR_dn
CDK4i_72hr_A3SS_FDR_dn_combine$combine <- paste(CDK4i_72hr_A3SS_FDR_dn$geneSymbol,
                                                CDK4i_72hr_A3SS_FDR_dn$longExonStart_0base,
                                                CDK4i_72hr_A3SS_FDR_dn$longExonEnd,
                                                CDK4i_72hr_A3SS_FDR_dn$shortES,
                                                CDK4i_72hr_A3SS_FDR_dn$shortEE,
                                                CDK4i_72hr_A3SS_FDR_dn$flankingES,
                                                CDK4i_72hr_A3SS_FDR_dn$flankingEE,
                                                sep = '_')

CDK4i_6day_A3SS_FDR_dn_combine <- CDK4i_6day_A3SS_FDR_dn
CDK4i_6day_A3SS_FDR_dn_combine$combine <- paste(CDK4i_6day_A3SS_FDR_dn$geneSymbol,
                                                CDK4i_6day_A3SS_FDR_dn$longExonStart_0base,
                                                CDK4i_6day_A3SS_FDR_dn$longExonEnd,
                                                CDK4i_6day_A3SS_FDR_dn$shortES,
                                                CDK4i_6day_A3SS_FDR_dn$shortEE,
                                                CDK4i_6day_A3SS_FDR_dn$flankingES,
                                                CDK4i_6day_A3SS_FDR_dn$flankingEE,
                                                sep = '_')

PRMT5i_72hr_A3SS_FDR_dn_combine <- PRMT5i_72hr_A3SS_FDR_dn
PRMT5i_72hr_A3SS_FDR_dn_combine$combine <- paste(PRMT5i_72hr_A3SS_FDR_dn$geneSymbol,
                                                 PRMT5i_72hr_A3SS_FDR_dn$longExonStart_0base,
                                                 PRMT5i_72hr_A3SS_FDR_dn$longExonEnd,
                                                 PRMT5i_72hr_A3SS_FDR_dn$shortES,
                                                 PRMT5i_72hr_A3SS_FDR_dn$shortEE,
                                                 PRMT5i_72hr_A3SS_FDR_dn$flankingES,
                                                 PRMT5i_72hr_A3SS_FDR_dn$flankingEE,
                                                 sep = '_')


## A5SS --------------------------------------------------------------------
CDK4i_72hr_A5SS_FDR_dn_combine <- CDK4i_72hr_A5SS_FDR_dn
CDK4i_72hr_A5SS_FDR_dn_combine$combine <- paste(CDK4i_72hr_A5SS_FDR_dn$geneSymbol,
                                                CDK4i_72hr_A5SS_FDR_dn$longExonStart_0base,
                                                CDK4i_72hr_A5SS_FDR_dn$longExonEnd,
                                                CDK4i_72hr_A5SS_FDR_dn$shortES,
                                                CDK4i_72hr_A5SS_FDR_dn$shortEE,
                                                CDK4i_72hr_A5SS_FDR_dn$flankingES,
                                                CDK4i_72hr_A5SS_FDR_dn$flankingEE,
                                                sep = '_')

CDK4i_6day_A5SS_FDR_dn_combine <- CDK4i_6day_A5SS_FDR_dn
CDK4i_6day_A5SS_FDR_dn_combine$combine <- paste(CDK4i_6day_A5SS_FDR_dn$geneSymbol,
                                                CDK4i_6day_A5SS_FDR_dn$longExonStart_0base,
                                                CDK4i_6day_A5SS_FDR_dn$longExonEnd,
                                                CDK4i_6day_A5SS_FDR_dn$shortES,
                                                CDK4i_6day_A5SS_FDR_dn$shortEE,
                                                CDK4i_6day_A5SS_FDR_dn$flankingES,
                                                CDK4i_6day_A5SS_FDR_dn$flankingEE,
                                                sep = '_')

PRMT5i_72hr_A5SS_FDR_dn_combine <- PRMT5i_72hr_A5SS_FDR_dn
PRMT5i_72hr_A5SS_FDR_dn_combine$combine <- paste(PRMT5i_72hr_A5SS_FDR_dn$geneSymbol,
                                                 PRMT5i_72hr_A5SS_FDR_dn$longExonStart_0base,
                                                 PRMT5i_72hr_A5SS_FDR_dn$longExonEnd,
                                                 PRMT5i_72hr_A5SS_FDR_dn$shortES,
                                                 PRMT5i_72hr_A5SS_FDR_dn$shortEE,
                                                 PRMT5i_72hr_A5SS_FDR_dn$flankingES,
                                                 PRMT5i_72hr_A5SS_FDR_dn$flankingEE,
                                                 sep = '_')

## MXE --------------------------------------------------------------------
CDK4i_72hr_MXE_FDR_dn_combine <- CDK4i_72hr_MXE_FDR_dn
CDK4i_72hr_MXE_FDR_dn_combine$combine <- paste(CDK4i_72hr_MXE_FDR_dn$geneSymbol,
                                               CDK4i_72hr_MXE_FDR_dn$X1stExonStart_0base,
                                               CDK4i_72hr_MXE_FDR_dn$X1stExonEnd,
                                               CDK4i_72hr_MXE_FDR_dn$X2ndExonStart_0base,
                                               CDK4i_72hr_MXE_FDR_dn$X2ndExonEnd,
                                               CDK4i_72hr_MXE_FDR_dn$dnstreamES,
                                               CDK4i_72hr_MXE_FDR_dn$dnstreamEE,
                                               CDK4i_72hr_MXE_FDR_dn$downstreamES,
                                               CDK4i_72hr_MXE_FDR_dn$downstreamEE,
                                               sep = '_')

CDK4i_6day_MXE_FDR_dn_combine <- CDK4i_6day_MXE_FDR_dn
CDK4i_6day_MXE_FDR_dn_combine$combine <- paste(CDK4i_6day_MXE_FDR_dn$geneSymbol,
                                               CDK4i_6day_MXE_FDR_dn$X1stExonStart_0base,
                                               CDK4i_6day_MXE_FDR_dn$X1stExonEnd,
                                               CDK4i_6day_MXE_FDR_dn$X2ndExonStart_0base,
                                               CDK4i_6day_MXE_FDR_dn$X2ndExonEnd,
                                               CDK4i_6day_MXE_FDR_dn$dnstreamES,
                                               CDK4i_6day_MXE_FDR_dn$dnstreamEE,
                                               CDK4i_6day_MXE_FDR_dn$downstreamES,
                                               CDK4i_6day_MXE_FDR_dn$downstreamEE,
                                               sep = '_')

PRMT5i_72hr_MXE_FDR_dn_combine <- PRMT5i_72hr_MXE_FDR_dn
PRMT5i_72hr_MXE_FDR_dn_combine$combine <- paste(PRMT5i_72hr_MXE_FDR_dn$geneSymbol,
                                                PRMT5i_72hr_MXE_FDR_dn$X1stExonStart_0base,
                                                PRMT5i_72hr_MXE_FDR_dn$X1stExonEnd,
                                                PRMT5i_72hr_MXE_FDR_dn$X2ndExonStart_0base,
                                                PRMT5i_72hr_MXE_FDR_dn$X2ndExonEnd,
                                                PRMT5i_72hr_MXE_FDR_dn$dnstreamES,
                                                PRMT5i_72hr_MXE_FDR_dn$dnstreamEE,
                                                PRMT5i_72hr_MXE_FDR_dn$downstreamES,
                                                PRMT5i_72hr_MXE_FDR_dn$downstreamEE,
                                                sep = '_')


## RI --------------------------------------------------------------------
CDK4i_72hr_RI_FDR_dn_combine <- CDK4i_72hr_RI_FDR_dn
CDK4i_72hr_RI_FDR_dn_combine$combine <- paste(CDK4i_72hr_RI_FDR_dn$geneSymbol,
                                              CDK4i_72hr_RI_FDR_dn$longExonStart_0base,
                                              CDK4i_72hr_RI_FDR_dn$longExonEnd,
                                              CDK4i_72hr_RI_FDR_dn$shortES,
                                              CDK4i_72hr_RI_FDR_dn$shortEE,
                                              CDK4i_72hr_RI_FDR_dn$flankingES,
                                              CDK4i_72hr_RI_FDR_dn$flankingEE,
                                              sep = '_')

CDK4i_6day_RI_FDR_dn_combine <- CDK4i_6day_RI_FDR_dn
CDK4i_6day_RI_FDR_dn_combine$combine <- paste(CDK4i_6day_RI_FDR_dn$geneSymbol,
                                              CDK4i_6day_RI_FDR_dn$longExonStart_0base,
                                              CDK4i_6day_RI_FDR_dn$longExonEnd,
                                              CDK4i_6day_RI_FDR_dn$shortES,
                                              CDK4i_6day_RI_FDR_dn$shortEE,
                                              CDK4i_6day_RI_FDR_dn$flankingES,
                                              CDK4i_6day_RI_FDR_dn$flankingEE,
                                              sep = '_')

PRMT5i_72hr_RI_FDR_dn_combine <- PRMT5i_72hr_RI_FDR_dn
PRMT5i_72hr_RI_FDR_dn_combine$combine <- paste(PRMT5i_72hr_RI_FDR_dn$geneSymbol,
                                               PRMT5i_72hr_RI_FDR_dn$longExonStart_0base,
                                               PRMT5i_72hr_RI_FDR_dn$longExonEnd,
                                               PRMT5i_72hr_RI_FDR_dn$shortES,
                                               PRMT5i_72hr_RI_FDR_dn$shortEE,
                                               PRMT5i_72hr_RI_FDR_dn$flankingES,
                                               PRMT5i_72hr_RI_FDR_dn$flankingEE,
                                               sep = '_')

## SE --------------------------------------------------------------------
CDK4i_72hr_SE_FDR_dn_combine <- CDK4i_72hr_SE_FDR_dn
CDK4i_72hr_SE_FDR_dn_combine$combine <- paste(CDK4i_72hr_SE_FDR_dn$geneSymbol,
                                              CDK4i_72hr_SE_FDR_dn$longExonStart_0base,
                                              CDK4i_72hr_SE_FDR_dn$longExonEnd,
                                              CDK4i_72hr_SE_FDR_dn$shortES,
                                              CDK4i_72hr_SE_FDR_dn$shortEE,
                                              CDK4i_72hr_SE_FDR_dn$flankingES,
                                              CDK4i_72hr_SE_FDR_dn$flankingEE,
                                              sep = '_')

CDK4i_6day_SE_FDR_dn_combine <- CDK4i_6day_SE_FDR_dn
CDK4i_6day_SE_FDR_dn_combine$combine <- paste(CDK4i_6day_SE_FDR_dn$geneSymbol,
                                              CDK4i_6day_SE_FDR_dn$longExonStart_0base,
                                              CDK4i_6day_SE_FDR_dn$longExonEnd,
                                              CDK4i_6day_SE_FDR_dn$shortES,
                                              CDK4i_6day_SE_FDR_dn$shortEE,
                                              CDK4i_6day_SE_FDR_dn$flankingES,
                                              CDK4i_6day_SE_FDR_dn$flankingEE,
                                              sep = '_')

PRMT5i_72hr_SE_FDR_dn_combine <- PRMT5i_72hr_SE_FDR_dn
PRMT5i_72hr_SE_FDR_dn_combine$combine <- paste(PRMT5i_72hr_SE_FDR_dn$geneSymbol,
                                               PRMT5i_72hr_SE_FDR_dn$longExonStart_0base,
                                               PRMT5i_72hr_SE_FDR_dn$longExonEnd,
                                               PRMT5i_72hr_SE_FDR_dn$shortES,
                                               PRMT5i_72hr_SE_FDR_dn$shortEE,
                                               PRMT5i_72hr_SE_FDR_dn$flankingES,
                                               PRMT5i_72hr_SE_FDR_dn$flankingEE,
                                               sep = '_')




# list of up-regulated splice events for UpSet plot ----------------------------
A3SS_up_lt = list(CDK4i_72hr = CDK4i_72hr_A3SS_FDR_up_combine$combine,
          CDK4i_6day = CDK4i_6day_A3SS_FDR_up_combine$combine,
          PRMT5i_72hr = PRMT5i_72hr_A3SS_FDR_up_combine$combine)

A5SS_up_lt = list(CDK4i_72hr = CDK4i_72hr_A5SS_FDR_up_combine$combine,
                  CDK4i_6day = CDK4i_6day_A5SS_FDR_up_combine$combine,
                  PRMT5i_72hr = PRMT5i_72hr_A5SS_FDR_up_combine$combine)

MXE_up_lt = list(CDK4i_72hr = CDK4i_72hr_MXE_FDR_up_combine$combine,
                  CDK4i_6day = CDK4i_6day_MXE_FDR_up_combine$combine,
                  PRMT5i_72hr = PRMT5i_72hr_MXE_FDR_up_combine$combine)

RI_up_lt = list(CDK4i_72hr = CDK4i_72hr_RI_FDR_up_combine$combine,
                  CDK4i_6day = CDK4i_6day_RI_FDR_up_combine$combine,
                  PRMT5i_72hr = PRMT5i_72hr_RI_FDR_up_combine$combine)

SE_up_lt = list(CDK4i_72hr = CDK4i_72hr_SE_FDR_up_combine$combine,
                  CDK4i_6day = CDK4i_6day_SE_FDR_up_combine$combine,
                  PRMT5i_72hr = PRMT5i_72hr_SE_FDR_up_combine$combine)

# list of dn-regulated splice events for dnSet plot ----------------------------
A3SS_dn_lt = list(CDK4i_72hr = CDK4i_72hr_A3SS_FDR_dn_combine$combine,
                  CDK4i_6day = CDK4i_6day_A3SS_FDR_dn_combine$combine,
                  PRMT5i_72hr = PRMT5i_72hr_A3SS_FDR_dn_combine$combine)

A5SS_dn_lt = list(CDK4i_72hr = CDK4i_72hr_A5SS_FDR_dn_combine$combine,
                  CDK4i_6day = CDK4i_6day_A5SS_FDR_dn_combine$combine,
                  PRMT5i_72hr = PRMT5i_72hr_A5SS_FDR_dn_combine$combine)

MXE_dn_lt = list(CDK4i_72hr = CDK4i_72hr_MXE_FDR_dn_combine$combine,
                 CDK4i_6day = CDK4i_6day_MXE_FDR_dn_combine$combine,
                 PRMT5i_72hr = PRMT5i_72hr_MXE_FDR_dn_combine$combine)

RI_dn_lt = list(CDK4i_72hr = CDK4i_72hr_RI_FDR_dn_combine$combine,
                CDK4i_6day = CDK4i_6day_RI_FDR_dn_combine$combine,
                PRMT5i_72hr = PRMT5i_72hr_RI_FDR_dn_combine$combine)

SE_dn_lt = list(CDK4i_72hr = CDK4i_72hr_SE_FDR_dn_combine$combine,
                CDK4i_6day = CDK4i_6day_SE_FDR_dn_combine$combine,
                PRMT5i_72hr = PRMT5i_72hr_SE_FDR_dn_combine$combine)


# Up-regulated splice events UpSet plots ----------------------------------
A3SS_up_m= make_comb_mat(A3SS_up_lt,mode = 'intersect')
A5SS_up_m= make_comb_mat(A5SS_up_lt,mode = 'intersect')
MXE_up_m= make_comb_mat(MXE_up_lt,mode = 'intersect')
RI_up_m= make_comb_mat(RI_up_lt,mode = 'intersect')
SE_up_m= make_comb_mat(SE_up_lt,mode = 'intersect')

# Down-regulated splice events dnSet plots ----------------------------------
A3SS_dn_m= make_comb_mat(A3SS_dn_lt,mode = 'intersect')
A5SS_dn_m= make_comb_mat(A5SS_dn_lt,mode = 'intersect')
MXE_dn_m= make_comb_mat(MXE_dn_lt,mode = 'intersect')
RI_dn_m= make_comb_mat(RI_dn_lt,mode = 'intersect')
SE_dn_m= make_comb_mat(SE_dn_lt,mode = 'intersect')


# Create UpSet Plots ------------------------------------------------------
## UpSet plot up-regulated -------------------------------------------------
UpSet_plot_upReg = rowAnnotation(
  "A3SS" = anno_barplot(comb_size(A3SS_up_m), 
                           gp = gpar(fill = "black"), width = unit(3, "cm")), 
  "A5SS" = anno_barplot(comb_size(A5SS_up_m), 
                             gp = gpar(fill = "black"), width = unit(3, "cm")), 
  "MXE" = anno_barplot(comb_size(MXE_up_m), 
                         gp = gpar(fill = "black"), width = unit(3, "cm")),
  "RI" = anno_barplot(comb_size(MXE_up_m), 
                       gp = gpar(fill = "black"), width = unit(3, "cm")),
  "SE" = anno_barplot(comb_size(MXE_up_m), 
                       gp = gpar(fill = "black"), width = unit(3, "cm")),
  gap = unit(2, "mm"), annotation_name_side = "bottom")

UpSet(t(A3SS_up_m), right_annotation = UpSet_plot_upReg)


## UpSet plot down-regulated -----------------------------------------------
UpSet_plot_dnReg = rowAnnotation(
  "A3SS" = anno_barplot(comb_size(A3SS_dn_m), 
                        gp = gpar(fill = "black"), width = unit(3, "cm")), 
  "A5SS" = anno_barplot(comb_size(A5SS_dn_m), 
                        gp = gpar(fill = "black"), width = unit(3, "cm")), 
  "MXE" = anno_barplot(comb_size(MXE_dn_m), 
                       gp = gpar(fill = "black"), width = unit(3, "cm")),
  "RI" = anno_barplot(comb_size(MXE_dn_m), 
                      gp = gpar(fill = "black"), width = unit(3, "cm")),
  "SE" = anno_barplot(comb_size(MXE_dn_m), 
                      gp = gpar(fill = "black"), width = unit(3, "cm")),
  gap = unit(2, "mm"), annotation_name_side = "bottom")

UpSet(t(A3SS_dn_m), right_annotation = UpSet_plot_dnReg)


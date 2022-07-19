library(msigdbr)
library(limma)
library(edgeR)
library(openxlsx)
library(EGSEA)

C5 <- msigdbr(species = "Homo sapiens", category = "C5")
C5_BP <- C5[C5$gs_subcat=="GO:BP",]
msigdb <- unique(C5_BP[,7:8])
head(msigdb)

#CDK4i DGE
CDK4i_DEGs <- read.xlsx("/Differential_gene_expression/edgeR/edgeR_sigchange_dge.xlsx",sheet=4)

unique_CDK4i_DEGs <- unique(CDK4i_DEGs[,1])
unique_CDK4i_DEGs_df <- data.frame("Gene_symbol" = unique_CDK4i_DEGs,"test" = 1)

unique_CDK4i_DEGs_ID <- merge(unique_CDK4i_DEGs_df,msigdb,by.x="Gene_symbol",by.y="human_gene_symbol")
unique_CDK4i_DEGs_ID <- unique_CDK4i_DEGs_ID[,c(1,3)]

gs.annots = buildIdx(entrezIDs=unique_CDK4i_DEGs_ID$human_entrez_gene, species="human", msigdb.gsets=c("c5"))
gs.annots

gsa.ora_CDK4i <- egsea.ora(unique_CDK4i_DEGs_ID$human_entrez_gene,gs.annots = gs.annots,report.dir = NULL,report = FALSE)
EGSEA_CDK4i_ora <- topSets(gsa.ora_CDK4i, gs.label = 'c5', sort.by = "p.adj", number = Inf, names.only = FALSE)

EGSEA_CDK4i_ora_adj <- EGSEA_CDK4i_ora[EGSEA_CDK4i_ora$p.adj<0.05,]
write.csv(EGSEA_CDK4i_ora_adj,"/EGSEA_CDK4i_ora_adj_BP_DEG.csv")


#PRMT5i DGE
PRMT5i_DEGs <- read.xlsx("/Differential_gene_expression/edgeR/edgeR_sigchange_dge.xlsx",sheet=2)

unique_PRMT5i_DEGs <- unique(PRMT5i_DEGs[,1])
unique_PRMT5i_DEGs_df <- data.frame("Gene_symbol" = unique_PRMT5i_DEGs,"test" = 1)

unique_PRMT5i_DEGs_ID <- merge(unique_PRMT5i_DEGs_df,msigdb,by.x="Gene_symbol",by.y="human_gene_symbol")
unique_PRMT5i_DEGs_ID <- unique_PRMT5i_DEGs_ID[,c(1,3)]

gs.annots = buildIdx(entrezIDs=unique_PRMT5i_DEGs_ID$human_entrez_gene, species="human", msigdb.gsets=c("c5"))
gs.annots

gsa.ora_PRMT5i <- egsea.ora(unique_PRMT5i_DEGs_ID$human_entrez_gene,gs.annots = gs.annots,report.dir = NULL,report = FALSE)
EGSEA_PRMT5i_ora <- topSets(gsa.ora_PRMT5i, gs.label = 'c5', sort.by = "p.adj", number = Inf, names.only = FALSE)

EGSEA_PRMT5i_ora_adj <- EGSEA_PRMT5i_ora[EGSEA_PRMT5i_ora$p.adj<0.05,]
write.csv(EGSEA_PRMT5i_ora_adj,"/EGSEA_PRMT5i_ora_adj_BP_DEG.csv")


#CDK4i DSG
CDK4i_6d_A3SS <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/CDK4i_6d_DMSO_72hr/A3SS.MATS.JC.txt")
CDK4i_6d_A5SS <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/CDK4i_6d_DMSO_72hr/A5SS.MATS.JC.txt")
CDK4i_6d_MXE <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/CDK4i_6d_DMSO_72hr/MXE.MATS.JC.txt")
CDK4i_6d_RI <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/CDK4i_6d_DMSO_72hr/RI.MATS.JC.txt")
CDK4i_6d_SE <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/CDK4i_6d_DMSO_72hr/SE.MATS.JC.txt")

CDK4i_6d_A3SS <- CDK4i_6d_A3SS[CDK4i_6d_A3SS$FDR<0.05,]
CDK4i_6d_A5SS <- CDK4i_6d_A5SS[CDK4i_6d_A5SS$FDR<0.05,]
CDK4i_6d_MXE <- CDK4i_6d_MXE[CDK4i_6d_MXE$FDR<0.05,]
CDK4i_6d_RI <- CDK4i_6d_RI[CDK4i_6d_RI$FDR<0.05,]
CDK4i_6d_SE <- CDK4i_6d_SE[CDK4i_6d_SE$FDR<0.05,]

CDK4i_6d_DSGs <- rbind(CDK4i_6d_A3SS[,c(3,23)],
                       CDK4i_6d_A5SS[,c(3,23)],
                       CDK4i_6d_MXE[,c(3,25)],
                       CDK4i_6d_RI[,c(3,23)],
                       CDK4i_6d_SE[,c(3,23)])
unique_CDK4i_6d_DSGs <- unique(CDK4i_6d_DSGs$geneSymbol)
unique_CDK4i_6d_DSGs_df <- data.frame("Gene_symbol" = unique_CDK4i_6d_DSGs,"test" = 1)

unique_CDK4i_6d_DSGs_ID <- merge(unique_CDK4i_6d_DSGs_df,msigdb,by.x="Gene_symbol",by.y="human_gene_symbol")
unique_CDK4i_6d_DSGs_ID <- unique_CDK4i_6d_DSGs_ID[,c(1,3)]

gs.annots = buildIdx(entrezIDs=unique_CDK4i_6d_DSGs_ID$human_entrez_gene, species="human", msigdb.gsets=c("C5"))
gs.annots

gsa.ora <- egsea.ora(unique_CDK4i_6d_DSGs_ID$human_entrez_gene,gs.annots = gs.annots,report.dir = NULL,report = FALSE)
EGSEA_CDK4i_ora <- topSets(gsa.ora, gs.label = 'c5', sort.by = "p.adj", number = Inf, names.only = FALSE)

EGSEA_CDK4i_ora_adj <- EGSEA_CDK4i_ora[EGSEA_CDK4i_ora$p.adj<0.05,]
write.csv(EGSEA_CDK4i_ora_adj,"/EGSEA_CDK4i_ora_adj_BP.csv")


#PRMT5i DSG
PRMT5i_A3SS <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/PRMT5i_72hr_DMSO_72hr/A3SS.MATS.JC.txt")
PRMT5i_A5SS <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/PRMT5i_72hr_DMSO_72hr/A5SS.MATS.JC.txt")
PRMT5i_MXE <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/PRMT5i_72hr_DMSO_72hr/MXE.MATS.JC.txt")
PRMT5i_RI <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/PRMT5i_72hr_DMSO_72hr/RI.MATS.JC.txt")
PRMT5i_SE <- read.delim("/rMATS/RL149_AL08_hg19_VaryRL_novelSS/PRMT5i_72hr_DMSO_72hr/SE.MATS.JC.txt")

PRMT5i_A3SS <- PRMT5i_A3SS[PRMT5i_A3SS$FDR<0.05,]
PRMT5i_A5SS <- PRMT5i_A5SS[PRMT5i_A5SS$FDR<0.05,]
PRMT5i_MXE <- PRMT5i_MXE[PRMT5i_MXE$FDR<0.05,]
PRMT5i_RI <- PRMT5i_RI[PRMT5i_RI$FDR<0.05,]
PRMT5i_SE <- PRMT5i_SE[PRMT5i_SE$FDR<0.05,]

PRMT5i_DSGs <- rbind(PRMT5i_A3SS[,c(3,23)],
                     PRMT5i_A5SS[,c(3,23)],
                     PRMT5i_MXE[,c(3,25)],
                     PRMT5i_RI[,c(3,23)],
                     PRMT5i_SE[,c(3,23)])

unique_PRMT5i_DSGs <- unique(PRMT5i_DSGs$geneSymbol)
unique_PRMT5i_6d_DSGs_df <- data.frame("Gene_symbol" = unique_PRMT5i_DSGs,"test" = 1)

unique_PRMT5i_6d_DSGs_ID <- merge(unique_PRMT5i_6d_DSGs_df,msigdb,by.x="Gene_symbol",by.y="human_gene_symbol")
unique_PRMT5i_6d_DSGs_ID <- unique_PRMT5i_6d_DSGs_ID[,c(1,3)]

gs.annots = buildIdx(entrezIDs=unique_PRMT5i_6d_DSGs_ID$human_entrez_gene, species="human", msigdb.gsets=c("c5"))
gs.annots

gsa.ora_PRMT5i <- egsea.ora(unique_PRMT5i_6d_DSGs_ID$human_entrez_gene,gs.annots = gs.annots,report.dir = NULL,report = FALSE)
EGSEA_PRMT5i_ora <- topSets(gsa.ora_PRMT5i, gs.label = 'c5', sort.by = "p.adj", number = Inf, names.only = FALSE)

EGSEA_PRMT5i_ora_adj <- EGSEA_PRMT5i_ora[EGSEA_PRMT5i_ora$p.adj<0.05,]
write.csv(EGSEA_PRMT5i_ora_adj,"/EGSEA_PRMT5i_ora_adj_BP.csv")

# 2. Differential gene expression analysis

## 2.1.Set Up
#The limma and edgeR packages are required for differential gene expression analysis.The pheatmap is used to create heatmap to visualize the differential gene expression results. The limma and edgeR uses different statistical model to determine differentially expressed genes. The following document shows steps required for edgeR and limma differential gene expression analysis.
library(limma)
library(edgeR)
library(pheatmap)
library(openxlsx)
library(knitr)


#Insert read count matrix resulted from HTSeq
setwd('../Read_quantification/')
count <-
  read.csv("HTSeq_matrix.csv",
           row.names = 1)


## 2.2. Pre-processing
#reate a DGEList object from the read count matrix using the DGEList function.The group option is used to group replicates with same conditions.The DGEList object contains the read count matrix and the information of each samples, including the grouping conditions, the libary size and the normalisation factors.
group <-
  factor(c(
    'DMSO',
    'DMSO',
    'PD991',
    'PD991',
    'CTX',
    'CTX',
    'PD991_6d',
    'PD991_6d'
  ))
dge_list <- DGEList(count, group = group)


#Lowly expressed genes in DGEList object are then filtered.Grouping information is included and the DGEList can be processed directly by the filterByExpr function.If grouping information is not included, it can be inserted using the group option or the design option. For the used data, the library size have dropped slightly after filtering.
keep <- filterByExpr(dge_list)
dge_list <- dge_list[keep, , keep.lib.sizes = FALSE]


#The DGEList object is then normalised with the default TMM method.Other method such as RLE and upperquartile are also availible.
dge_list <- calcNormFactors(dge_list)


#Create model matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


#Create contrast matrix
cont.matrix <- makeContrasts(
  CDK4i_72hr_DMSO_72hr = PD991 - DMSO,
  PRMT5i_72hr_DMSO_72hr = CTX - DMSO,
  CDK4i_72hr_PRMT5i_72hr = (PD991 - DMSO) - (CTX - DMSO),
  CDK4i_6d_DMSO_72hr = PD991_6d - DMSO,
  CDK4i_6d_CDK4i_72hr = (PD991_6d - DMSO) - (PD991 - DMSO),
  CDK4i_6d_PRMT5i_72hr = (PD991_6d - DMSO) - (CTX - DMSO),
  levels = design
)
kable(cont.matrix)


## 2.3. edgeR differential gene expression analysis
### Estimate dispersion
#Estimate the dispersion between sample genes, which is required for fitting a generalised linear model. The estimateDisp() function estimates the common, trended and tagwise dispersion. Alternatively, it can be specified using the estimateGLMCommonDisp, estimateGLMTrendedDisp and estimateGLMTagwiseDisp.
edgeR_dge_list <- estimateDisp(dge_list, design, robust = TRUE)


### Fitting model
#glmQLFit() fits the read count of the genes in to a negative binomial model and produce object of the class DGEGLM.Details of the DGEGLM class can be found under the help section by searching DGEGLM -class.
edgeR_fit <- glmQLFit(edgeR_dge_list, design, robust = TRUE)


#The DGEGLM object can be passed to the glmQLFTest() function to test for differentially expressed genes.The test results for each comparisons are stored individually.
edgeR_deg <- glmQLFTest(edgeR_fit, contrast = cont.matrix)
edgeR_deg_result <-
  topTags(
    edgeR_deg,
    n = nrow(edgeR_deg$table),
    sort.by = 'none',
    p.value = 0.05
  )$table

edgeR_dge_CDK4i_72hr_DMSO_72hr <-
  glmQLFTest(edgeR_fit, contrast = cont.matrix[, 1])
edgeR_dge_PRMT5i_72hr_DMSO_72hr <-
  glmQLFTest(edgeR_fit, contrast = cont.matrix[, 2])
edgeR_dge_CDK4i_72hr_PRMT5i_72hr <-
  glmQLFTest(edgeR_fit, contrast = cont.matrix[, 3])
edgeR_dge_CDK4i_6d_DMSO_72hr <-
  glmQLFTest(edgeR_fit, contrast = cont.matrix[, 4])
edgeR_dge_CDK4i_6d_CDK4i_72hr <-
  glmQLFTest(edgeR_fit, contrast = cont.matrix[, 5])
edgeR_dge_CDK4i_6d_PRMT5i_72hr <-
  glmQLFTest(edgeR_fit, contrast = cont.matrix[, 6])


### Extract results
#The number of differentially expressed genes is summarized, including the number of non -significant, significantly up regulated and significantly down regulated differentially expressed genes.
edgeR_dge_summary <-
  cbind(
    summary(decideTests(edgeR_dge_CDK4i_72hr_DMSO_72hr)),
    summary(decideTests(edgeR_dge_PRMT5i_72hr_DMSO_72hr)),
    summary(decideTests(edgeR_dge_CDK4i_72hr_PRMT5i_72hr)),
    summary(decideTests(edgeR_dge_CDK4i_6d_DMSO_72hr)),
    summary(decideTests(edgeR_dge_CDK4i_6d_CDK4i_72hr)),
    summary(decideTests(edgeR_dge_CDK4i_6d_PRMT5i_72hr))
  )
kable(edgeR_dge_summary)


#The entire list of differentially expressed genes can be extracted using the topTags() function. FDR is calculated for each gene. Significantly changed genes are selected with FDR < 0.05.
edgeR_res_CDK4i_72hr_DMSO_72hr <-
  topTags(edgeR_dge_CDK4i_72hr_DMSO_72hr,
          n = nrow(edgeR_dge_CDK4i_72hr_DMSO_72hr))

edgeR_sigchange_CDK4i_72hr_DMSO_72hr <-
  edgeR_res_CDK4i_72hr_DMSO_72hr$table[edgeR_res_CDK4i_72hr_DMSO_72hr$table$FDR < 0.05,]

edgeR_res_PRMT5i_72hr_DMSO_72hr <-
  topTags(edgeR_dge_PRMT5i_72hr_DMSO_72hr,
          n = nrow(edgeR_dge_PRMT5i_72hr_DMSO_72hr))
edgeR_sigchange_PRMT5i_72hr_DMSO_72hr <-
  edgeR_res_PRMT5i_72hr_DMSO_72hr$table[edgeR_res_PRMT5i_72hr_DMSO_72hr$table$FDR < 0.05,]

edgeR_res_CDK4i_72hr_PRMT5i_72hr <-
  topTags(edgeR_dge_CDK4i_72hr_PRMT5i_72hr,
          n = nrow(edgeR_dge_CDK4i_72hr_PRMT5i_72hr))
edgeR_sigchange_CDK4i_72hr_PRMT5i_72hr <-
  edgeR_res_CDK4i_72hr_PRMT5i_72hr$table[edgeR_res_CDK4i_72hr_PRMT5i_72hr$table$FDR < 0.05,]

edgeR_res_CDK4i_6d_DMSO_72hr <-
  topTags(edgeR_dge_CDK4i_6d_DMSO_72hr,
          n = nrow(edgeR_dge_CDK4i_6d_DMSO_72hr))
edgeR_sigchange_CDK4i_6d_DMSO_72hr <-
  edgeR_res_CDK4i_6d_DMSO_72hr$table[edgeR_res_CDK4i_6d_DMSO_72hr$table$FDR < 0.05,]

edgeR_res_CDK4i_6d_CDK4i_72hr <-
  topTags(edgeR_dge_CDK4i_6d_CDK4i_72hr,
          n = nrow(edgeR_dge_CDK4i_6d_CDK4i_72hr))
edgeR_sigchange_CDK4i_6d_CDK4i_72hr <-
  edgeR_res_CDK4i_6d_CDK4i_72hr$table[edgeR_res_CDK4i_6d_CDK4i_72hr$table$FDR < 0.05,]

edgeR_res_CDK4i_6d_PRMT5i_72hr <-
  topTags(edgeR_dge_CDK4i_6d_PRMT5i_72hr,
          n = nrow(edgeR_dge_CDK4i_6d_PRMT5i_72hr))
edgeR_sigchange_CDK4i_6d_PRMT5i_72hr <-
  edgeR_res_CDK4i_6d_PRMT5i_72hr$table[edgeR_res_CDK4i_6d_PRMT5i_72hr$table$FDR < 0.05,]


#The significantly differentially expressed genes are further splitted by the regulation directions, which is determined by the fold change direction.
edgeR_sigDownReg_CDK4i_72hr_DMSO_72hr <-
  edgeR_sigchange_CDK4i_72hr_DMSO_72hr[edgeR_sigchange_CDK4i_72hr_DMSO_72hr$logFC < 0,]
edgeR_sigUpReg_CDK4i_72hr_DMSO_72hr <-
  edgeR_sigchange_CDK4i_72hr_DMSO_72hr[edgeR_sigchange_CDK4i_72hr_DMSO_72hr$logFC > 0,]

edgeR_sigDownReg_PRMT5i_72hr_DMSO_72hr <-
  edgeR_sigchange_PRMT5i_72hr_DMSO_72hr[edgeR_sigchange_PRMT5i_72hr_DMSO_72hr$logFC < 0,]
edgeR_sigUpReg_PRMT5i_72hr_DMSO_72hr <-
  edgeR_sigchange_PRMT5i_72hr_DMSO_72hr[edgeR_sigchange_PRMT5i_72hr_DMSO_72hr$logFC > 0,]

edgeR_sigDownReg_CDK4i_72hr_PRMT5i_72hr <-
  edgeR_sigchange_CDK4i_72hr_PRMT5i_72hr[edgeR_sigchange_CDK4i_72hr_PRMT5i_72hr$logFC < 0,]
edgeR_sigUpReg_CDK4i_72hr_PRMT5i_72hr <-
  edgeR_sigchange_CDK4i_72hr_PRMT5i_72hr[edgeR_sigchange_CDK4i_72hr_PRMT5i_72hr$logFC > 0,]

edgeR_sigDownReg_CDK4i_6d_DMSO_72hr <-
  edgeR_sigchange_CDK4i_6d_DMSO_72hr[edgeR_sigchange_CDK4i_6d_DMSO_72hr$logFC < 0,]
edgeR_sigUpReg_CDK4i_6d_DMSO_72hr <-
  edgeR_sigchange_CDK4i_6d_DMSO_72hr[edgeR_sigchange_CDK4i_6d_DMSO_72hr$logFC > 0,]

edgeR_sigDownReg_CDK4i_6d_CDK4i_72hr <-
  edgeR_sigchange_CDK4i_6d_CDK4i_72hr[edgeR_sigchange_CDK4i_6d_CDK4i_72hr$logFC < 0,]
edgeR_sigUpReg_CDK4i_6d_CDK4i_72hr <-
  edgeR_sigchange_CDK4i_6d_CDK4i_72hr[edgeR_sigchange_CDK4i_6d_CDK4i_72hr$logFC > 0,]

edgeR_sigDownReg_CDK4i_6d_PRMT5i_72hr <-
  edgeR_sigchange_CDK4i_6d_PRMT5i_72hr[edgeR_sigchange_CDK4i_6d_PRMT5i_72hr$logFC < 0,]
edgeR_sigUpReg_CDK4i_6d_PRMT5i_72hr <-
  edgeR_sigchange_CDK4i_6d_PRMT5i_72hr[edgeR_sigchange_CDK4i_6d_PRMT5i_72hr$logFC > 0,]


#Three excel files containing list of significantly differentially expressed genes under different contrast are generated.The edgeR_sigchange_dge.xlsx contains all significantly differentially expressed genes, the edgeR_sigDownReg_dge.xlsx contains down regulated significantly differentially expressed genes and the edgeR_sigUpReg_dge.xlsx contains up regulated significantly differentially expressed genes.The file contains log fold change, log count per million across all samples, F statistics, p - values and the false discovery rate.
edgeR_sigchange_dge <-
  list(
    'CDK4i_72hr_DMSO_72hr' = edgeR_sigchange_CDK4i_72hr_DMSO_72hr,
    'PRMT5i_72hr_DMSO_72hr' = edgeR_sigchange_PRMT5i_72hr_DMSO_72hr,
    'CDK4i_72hr_PRMT5i_72hr' = edgeR_sigchange_CDK4i_72hr_PRMT5i_72hr,
    'CDK4i_6d_DMSO_72hr' = edgeR_sigchange_CDK4i_6d_DMSO_72hr,
    'CDK4i_6d_CDK4i_72hr' = edgeR_sigchange_CDK4i_6d_CDK4i_72hr,
    'CDK4i_6d_PRMT5i_72hr' = edgeR_sigchange_CDK4i_6d_PRMT5i_72hr
  )

edgeR_sigDownReg_dge <-
  list(
    'CDK4i_72hr_DMSO_72hr' = edgeR_sigDownReg_CDK4i_72hr_DMSO_72hr,
    'PRMT5i_72hr_DMSO_72hr' = edgeR_sigDownReg_PRMT5i_72hr_DMSO_72hr,
    'CDK4i_72hr_PRMT5i_72hr' = edgeR_sigDownReg_CDK4i_72hr_PRMT5i_72hr,
    'CDK4i_6d_DMSO_72hr' = edgeR_sigDownReg_CDK4i_6d_DMSO_72hr,
    'CDK4i_6d_CDK4i_72hr' = edgeR_sigDownReg_CDK4i_6d_CDK4i_72hr,
    'CDK4i_6d_PRMT5i_72hr' = edgeR_sigDownReg_CDK4i_6d_PRMT5i_72hr
  )

edgeR_sigUpReg_dge <-
  list(
    'CDK4i_72hr_DMSO_72hr' = edgeR_sigUpReg_CDK4i_72hr_DMSO_72hr,
    'PRMT5i_72hr_DMSO_72hr' = edgeR_sigUpReg_PRMT5i_72hr_DMSO_72hr,
    'CDK4i_72hr_PRMT5i_72hr' = edgeR_sigUpReg_CDK4i_72hr_PRMT5i_72hr,
    'CDK4i_6d_DMSO_72hr' = edgeR_sigUpReg_CDK4i_6d_DMSO_72hr,
    'CDK4i_6d_CDK4i_72hr' = edgeR_sigUpReg_CDK4i_6d_CDK4i_72hr,
    'CDK4i_6d_PRMT5i_72hr' = edgeR_sigUpReg_CDK4i_6d_PRMT5i_72hr
  )

write.xlsx(edgeR_sigchange_dge,
           'edgeR/edgeR_sigchange_dge.xlsx',
           row.names = TRUE)
write.xlsx(edgeR_sigDownReg_dge,
           'edgeR/edgeR_sigDownReg_dge.xlsx',
           row.names = TRUE)
write.xlsx(edgeR_sigUpReg_dge,
           'edgeR/edgeR_sigUpReg_dge.xlsx',
           row.names = TRUE)


## 2.4. limma differential gene expression analysis
### voom transformation
#Transform to log (count per million).
v = voom(dge_list, design = design)


### Fitting model
#lmFit fits linear model for each gene using weighted least squares for each gene
limma_fit = lmFit(v, design)
head(coef(limma_fit))


#Estimate contrast for each gene and each comparison
tmp <- contrasts.fit(limma_fit, cont.matrix)


#Empirical Bayes smoothed the standard errors
tmp <- eBayes(tmp)


### Extract results
#Extract significantly differentially expressed genes
limma_deg <-
  topTable(tmp,
           number = Inf,
           sort.by = 'none',
           p.value = 0.05)


#The number of differentially expressed genes is summarized, including the number of non significant, significantly up regulated and significantly down regulated differentially expressed genes.
limma_summary <- decideTests(tmp)
kable(summary(limma_summary))


#The entire list of differentially expressed genes can be extracted using the topTags() function. FDR is not calculated, significantly changed genes are selected with p value < 0.05.
limma_res_CDK4i_72hr_DMSO_72hr <-
  topTable(
    tmp,
    coef = 1,
    number = Inf,
    sort.by = 'none',
    p.value = 0.05
  )

limma_res_PRMT5i_72hr_DMSO_72hr <-
  topTable(
    tmp,
    coef = 2,
    number = Inf,
    sort.by = 'none',
    p.value = 0.05
  )

limma_res_CDK4i_72hr_PRMT5i_72hr <-
  topTable(
    tmp,
    coef = 3,
    number = Inf,
    sort.by = 'none',
    p.value = 0.05
  )

limma_res_CDK4i_6d_DMSO_72hr <-
  topTable(
    tmp,
    coef = 4,
    number = Inf,
    sort.by = 'none',
    p.value = 0.05
  )

limma_res_CDK4i_6d_CDK4i_72hr <-
  topTable(
    tmp,
    coef = 5,
    number = Inf,
    sort.by = 'none',
    p.value = 0.05
  )

limma_res_CDK4i_6d_PRMT5i_72hr <-
  topTable(
    tmp,
    coef = 6,
    number = Inf,
    sort.by = 'none',
    p.value = 0.05
  )


#The significantly differentially expressed genes are further splitted by the regulation directions, which is determined by the fold change direction.
limma_sigDownReg_CDK4i_72hr_DMSO_72hr <-
  limma_res_CDK4i_72hr_DMSO_72hr[limma_res_CDK4i_72hr_DMSO_72hr$logFC < 0,]
limma_sigUpReg_CDK4i_72hr_DMSO_72hr <-
  limma_res_CDK4i_72hr_DMSO_72hr[limma_res_CDK4i_72hr_DMSO_72hr$logFC > 0,]

limma_sigDownReg_PRMT5i_72hr_DMSO_72hr <-
  limma_res_PRMT5i_72hr_DMSO_72hr[limma_res_PRMT5i_72hr_DMSO_72hr$logFC < 0,]
limma_sigUpReg_PRMT5i_72hr_DMSO_72hr <-
  limma_res_PRMT5i_72hr_DMSO_72hr[limma_res_PRMT5i_72hr_DMSO_72hr$logFC > 0,]

limma_sigDownReg_CDK4i_72hr_PRMT5i_72hr <-
  limma_res_CDK4i_72hr_PRMT5i_72hr[limma_res_CDK4i_72hr_PRMT5i_72hr$logFC < 0,]
limma_sigUpReg_CDK4i_72hr_PRMT5i_72hr <-
  limma_res_CDK4i_72hr_PRMT5i_72hr[limma_res_CDK4i_72hr_PRMT5i_72hr$logFC > 0,]

limma_sigDownReg_CDK4i_6d_DMSO_72hr <-
  limma_res_CDK4i_6d_DMSO_72hr[limma_res_CDK4i_6d_DMSO_72hr$logFC < 0,]
limma_sigUpReg_CDK4i_6d_DMSO_72hr <-
  limma_res_CDK4i_6d_DMSO_72hr[limma_res_CDK4i_6d_DMSO_72hr$logFC > 0,]

limma_sigDownReg_CDK4i_6d_CDK4i_72hr <-
  limma_res_CDK4i_6d_CDK4i_72hr[limma_res_CDK4i_6d_CDK4i_72hr$logFC < 0,]
limma_sigUpReg_CDK4i_6d_CDK4i_72hr <-
  limma_res_CDK4i_6d_CDK4i_72hr[limma_res_CDK4i_6d_CDK4i_72hr$logFC > 0,]

limma_sigDownReg_CDK4i_6d_PRMT5i_72hr <-
  limma_res_CDK4i_6d_PRMT5i_72hr[limma_res_CDK4i_6d_PRMT5i_72hr$logFC < 0,]
limma_sigUpReg_CDK4i_6d_PRMT5i_72hr <-
  limma_res_CDK4i_6d_PRMT5i_72hr[limma_res_CDK4i_6d_PRMT5i_72hr$logFC > 0,]


#Three excel files containing list of significantly differentially expressed genes under different contrast are generated.The limma_sigchange_dge.xlsx contains all significantly differentially expressed genes, the limma_sigDownReg_dge.xlsx contains down regulated significantly differentially expressed genes and the limma_sigUpReg_dge.xlsx contains up regulated significantly differentially expressed genes.The file contains log fold change, average expression across samples, t statistics, p - values and the B - statistics.
limma_sigchange_dge <-
  list(
    'CDK4i_72hr_DMSO_72hr' = limma_res_CDK4i_72hr_DMSO_72hr,
    'PRMT5i_72hr_DMSO_72hr' = limma_res_PRMT5i_72hr_DMSO_72hr,
    'CDK4i_72hr_PRMT5i_72hr' = limma_res_CDK4i_72hr_PRMT5i_72hr,
    'CDK4i_6d_DMSO_72hr' = limma_res_CDK4i_6d_DMSO_72hr,
    'CDK4i_6d_CDK4i_72hr' = limma_res_CDK4i_6d_CDK4i_72hr,
    'CDK4i_6d_PRMT5i_72hr' = limma_res_CDK4i_6d_PRMT5i_72hr
  )

limma_sigDownReg_dge <-
  list(
    'CDK4i_72hr_DMSO_72hr' = limma_sigDownReg_CDK4i_72hr_DMSO_72hr,
    'PRMT5i_72hr_DMSO_72hr' = limma_sigDownReg_PRMT5i_72hr_DMSO_72hr,
    'CDK4i_72hr_PRMT5i_72hr' = limma_sigDownReg_CDK4i_72hr_PRMT5i_72hr,
    'CDK4i_6d_DMSO_72hr' = limma_sigDownReg_CDK4i_6d_DMSO_72hr,
    'CDK4i_6d_CDK4i_72hr' = limma_sigDownReg_CDK4i_6d_CDK4i_72hr,
    'CDK4i_6d_PRMT5i_72hr' = limma_sigDownReg_CDK4i_6d_PRMT5i_72hr
  )

limma_sigUpReg_dge <-
  list(
    'CDK4i_72hr_DMSO_72hr' = limma_sigUpReg_CDK4i_72hr_DMSO_72hr,
    'PRMT5i_72hr_DMSO_72hr' = limma_sigUpReg_PRMT5i_72hr_DMSO_72hr,
    'CDK4i_72hr_PRMT5i_72hr' = limma_sigUpReg_CDK4i_72hr_PRMT5i_72hr,
    'CDK4i_6d_DMSO_72hr' = limma_sigUpReg_CDK4i_6d_DMSO_72hr,
    'CDK4i_6d_CDK4i_72hr' = limma_sigUpReg_CDK4i_6d_CDK4i_72hr,
    'CDK4i_6d_PRMT5i_72hr' = limma_sigUpReg_CDK4i_6d_PRMT5i_72hr
  )

write.xlsx(limma_sigchange_dge,
           'limma/limma_sigchange_dge.xlsx',
           row.names = TRUE)
write.xlsx(limma_sigDownReg_dge,
           'limma/limma_sigDownReg_dge.xlsx',
           row.names = TRUE)
write.xlsx(limma_sigUpReg_dge,
           'limma/limma_sigUpReg_dge.xlsx',
           row.names = TRUE)


## 2.5. Result visualisation
#A heatmap showing all differentially expressed genes  and their fold change detected by edgeR.
pheatmap(
  edgeR_deg_result[, c(1, 2, 4)],
  cellwidth = 50,
  cluster_rows = 0,
  annotation_names_row = 1,
  labels_col = c(
    'CDK4i_72hr_DMSO_72hr',
    'PRMT5i_72hr_DMSO_72hr',
    'CDK4i_6d_DMSO_72hr'
  ),
  main = 'edgeR DEGs logFC'
)


#A heatmap showing all differentially expressed genes  and their fold change detected by limma.
pheatmap(
  limma_deg[, c(1, 2, 4)],
  cellwidth = 50,
  cluster_rows = 0,
  annotation_names_row = 1,
  main = 'limma DEGs logFC'
)


## 2.6. limma result comparison
#The number of differentially expressed genes detected by edgeR and limma are very similar. Gene set enrichment analysis is  performed next using the limma based tool, EGSEA. Thus, the differential gene expression result from limma is selected and compared.Venn diagrams are presented to compare the differentially expressed genes under different inhibitions.
colnames(limma_summary)[c(1, 2, 4)] <-
  c('CDK4i_72hr', 'PRMT5i_72hr', 'CDK4i_6d')

par(mfrow = c(2, 2))
vennDiagram(
  limma_summary[, 1:2],
  include = c('up', 'down'),
  cex = 0.8,
  main = 'CDK4i_72hr vs PRMT5i_72hr'
)

vennDiagram(
  limma_summary[, c(1, 4)],
  include = c('up', 'down'),
  cex = 0.8,
  main = 'CDK4i_72hr vs CDK4i_6d'
)

vennDiagram(
  limma_summary[, c(4, 2)],
  include = c('up', 'down'),
  cex = 0.8,
  main = 'CDK4i_6d vs PRMT5i_72hr'
)

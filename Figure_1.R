library(openxlsx)
library(ggplot2)

#Figure 1b
CDK4i_72hr_down <- read.xlsx("../edgeR/edgeR_sigDownReg_dge.xlsx",sheet=1)
PRMT5i_72hr_down <- read.xlsx("../edgeR_sigDownReg_dge.xlsx",sheet=2)
CDK4i_6d_down <- read.xlsx("../edgeR_sigDownReg_dge.xlsx",sheet=4)

CDK4i_72hr_up <- read.xlsx("../edgeR_sigUpReg_dge.xlsx",sheet=1)
PRMT5i_72hr_up <- read.xlsx("../edgeR_sigUpReg_dge.xlsx",sheet=2)
CDK4i_6d_up <- read.xlsx("../edgeR_sigUpReg_dge.xlsx",sheet=4)

DGE_count <- c(nrow(CDK4i_72hr_down),
               nrow(CDK4i_72hr_up),
               nrow(CDK4i_6d_down),
               nrow(CDK4i_6d_up),
               nrow(PRMT5i_72hr_down),
               nrow(PRMT5i_72hr_up))
Treatment <- c(rep("CDK4/6i (72hr)",2),rep("CDK4/6i (6d)",2),rep("PRMT5i (72hr)",2))
direction <- rep(c("down","up"),3)
DGE_count_df <- data.frame("Treatment"=Treatment,
                           "Direction"=direction,
                           "Number of DEGs"=DGE_count)
DGE_count_df$Treatment <- factor(DGE_count_df$Treatment, unique(DGE_count_df$Treatment),ordered = TRUE)

cbp3 <- c("orange1","deepskyblue2")
ggplot(DGE_count_df,aes(fill=direction,y=DGE_count,x=Treatment))+
  geom_bar(position="dodge", stat="identity",width = 0.8)+
  scale_fill_manual(values = cbp3)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  xlab("")+
  ylab("Number of DEGs")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title.y=element_text(size=15,face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        panel.background = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.1,0.85),
        legend.text = element_text(size = 15),
        legend.title = element_blank())

#Figure 1c
CDk4i_72hr_CDK4i_6d_overlap_up <- length(c(CDK4i_72hr_up[,1])[c(CDK4i_72hr_up[,1] %in% CDK4i_6d_up[,1])])
CDk4i_72hr_CDK4i_6d_overlap_down <- length(c(CDK4i_72hr_down[,1])[c(CDK4i_72hr_down[,1] %in% CDK4i_6d_down[,1])])

CDk4i_72hr_only_not_CDK4i_6d_up <- length(c(CDK4i_72hr_up[,1])[!c(CDK4i_72hr_up[,1] %in% CDK4i_6d_up[,1])])
CDk4i_72hr_only_not_CDK4i_6d_down <- length(c(CDK4i_72hr_down[,1])[!c(CDK4i_72hr_down[,1] %in% CDK4i_6d_down[,1])])

CDk4i_6d_only_not_CDK4i_72hr_up <- length(c(CDK4i_6d_up[,1])[!c(CDK4i_6d_up[,1] %in% CDK4i_72hr_up[,1])])
CDk4i_6d_only_not_CDK4i_72hr_down <- length(c(CDK4i_6d_down[,1])[!c(CDK4i_6d_down[,1] %in% CDK4i_72hr_down[,1])])


#Figure 1c
PRMT5i_72hr_CDK4i_6d_overlap_up <- length(c(PRMT5i_72hr_up[,1])[c(PRMT5i_72hr_up[,1] %in% CDK4i_6d_up[,1])])
PRMT5i_72hr_CDK4i_6d_overlap_down <- length(c(PRMT5i_72hr_down[,1])[c(PRMT5i_72hr_down[,1] %in% CDK4i_6d_down[,1])])

PRMT5i_72hr_only_not_CDK4i_6d_up <- length(c(PRMT5i_72hr_up[,1])[!c(PRMT5i_72hr_up[,1] %in% CDK4i_6d_up[,1])])
PRMT5i_72hr_only_not_CDK4i_6d_down <- length(c(PRMT5i_72hr_down[,1])[!c(PRMT5i_72hr_down[,1] %in% CDK4i_6d_down[,1])])

CDk4i_6d_only_not_PRMT5i_72hr_up <- length(c(CDK4i_6d_up[,1])[!c(CDK4i_6d_up[,1] %in% PRMT5i_72hr_up[,1])])
CDk4i_6d_only_not_PRMT5i_72hr_down <- length(c(CDK4i_6d_down[,1])[!c(CDK4i_6d_down[,1] %in% PRMT5i_72hr_down[,1])])


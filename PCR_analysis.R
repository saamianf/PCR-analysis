library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
library(pheatmap)

## Edit wd and file name as needed 
setwd('/Users/sf/Library/CloudStorage/Box-Box/Spring2026/McConkey Lab Rotation/PCRs/03102026/')
pcr_data <- read.csv('SF_03102026_Quantification_Cq_Results.csv') 

pcr_data_analysis <- pcr_data[,c(4,6:8)] %>% group_by(Target, Sample, Biological.Set.Name) %>% 
  mutate(id = row_number()) %>% 
  pivot_wider(names_from = id, values_from = Cq, names_prefix = "Cq_") 
colnames(pcr_data_analysis) <- c('Target', 'Sample', 'Condition', 'Cq_1', 'Cq_2', 'Cq_3')

pcr_data_analysis$Cq_avg <- rowMeans(pcr_data_analysis[,c(4:6)], na.rm = TRUE)
pcr_data_analysis <- pcr_data_analysis %>% group_by(Sample, Condition) %>%
  mutate(Cq_ref = Cq_avg[Target == "GAPDH"], dT = Cq_avg - Cq_ref) %>% ungroup()
pcr_data_analysis <- pcr_data_analysis %>% group_by(Target, Sample) %>% ## using the 2D sample's mean as ctrl
  mutate(control_dT = if (any(Condition == "2D")) {
    mean(dT[Condition == '2D'], na.rm = TRUE) } else { NA_real_},
    ddT = dT - control_dT) %>% ungroup()
pcr_data_analysis$FC <- 2^(-pcr_data_analysis$ddT)

## update the file save name as needed 
write_xlsx(pcr_data_analysis, 'PCR_data_R_SF_03102026.xlsx')

pcr_data_analysis$Name <- paste(pcr_data_analysis$Sample,pcr_data_analysis$Condition)
graph_data <- pcr_data_analysis[c(2:7, 9:14, 16:21, 23:28),c(1:3,13,9,11)] %>% group_by(Target) 

graph_data_dT <- graph_data[,c(1,4:5)] %>% pivot_wider(names_from = Name, values_from = dT) %>% as.data.frame(graph_data_dT)
rownames(graph_data_dT) <- graph_data_dT$Target
graph_data_dT <- graph_data_dT[,c(2:5)]

graph_data_ddT <- graph_data[,c(1,4,6)] %>% pivot_wider(names_from = Name, values_from = ddT) %>% as.data.frame(graph_data_ddT)
rownames(graph_data_ddT) <- graph_data_ddT$Target
graph_data_ddT <- graph_data_ddT[,c(2:5)]

pheatmap(graph_data_dT, cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 20,
  filename = 'Graphs/dT_heatmap.png', width=8, height=8)
pheatmap(graph_data_ddT, cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 20,
  filename = 'Graphs/ddT_heatmap.png', width=8, height=8)

primers <- c('ZEB1', 'p63 #38', 'p63 #39', 'p63 #40', 'KRT20')
cols <- c( "#56B4E9", "#009E73",)
for (i in primers){
  dT <- ggplot(subset(graph_data, Target == i), aes(x=Sample, y=dT, fill = Condition)) + 
    geom_col(position = "dodge") + ylab('dT Value') + xlab('Samples') + scale_fill_manual(values = cols) +
    ggtitle(paste('PCR Expression Levels for', i)) + theme(text = element_text(size = 20))
  ggsave(path='Graphs', paste(i, '_bargraph_dT.png', sep=''), width=8, height=8)
}
for (i in primers){
  dT <- ggplot(subset(graph_data, Target == i), aes(x=Sample, y=abs(ddT), fill = Condition)) + 
    geom_col(position = "dodge") + ylab('ddT Value') + xlab('Samples') + scale_fill_manual(values = cols) +
    ggtitle(paste('PCR Expression Levels for', i)) + theme(text = element_text(size = 20))
  ggsave(path='Graphs', paste(i, '_bargraph_ddT.png', sep=''), width=8, height=8)
}


# https://toptipbio.com/delta-delta-ct-pcr/

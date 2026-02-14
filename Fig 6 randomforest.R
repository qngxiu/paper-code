##五样本细菌pcoa绘制

rm(list = ls())
getwd()
# setwd("H:\\东北地理与农业生态研究所\\重金属投稿\\原始数据--图表--对应\\原始数据--图表--对应\\Figure5 网络图和zipi\\网络图/真菌文件")
xiufu <- read.csv("total_fu_拓扑.csv",header = T,row.names = 1)
 xiufu <- xiufu[-14,-c(1,2,11,15)]
df1 <- xiufu
merged_df <- df1
print(merged_df)
# set.seed(1122)
library(vegan)
 data <- t(merged_df)
 # data <- decostand(data, method = "hellinger")
# group <- read.csv("group.csv",header = T,row.names = 1)
# colnames(group) <- "Group"
bray <- vegdist(data, method = 'bray')
bray <- as.matrix(bray)
pcoa <- cmdscale(bray, k = 3, eig = T)
pcoa_data <- data.frame({pcoa$point})
pcoa_data$Sample_ID <- rownames(pcoa_data)
names(pcoa_data)[1:3] <- paste0("PCoA", 1:3)




 rm(list = ls())
set.seed(123)
setwd("H:\\东北地理与农业生态研究所\\重金属投稿\\原始数据--图表--对应\\原始数据--图表--对应\\Figure6 SEM\\随机森林")
# # plant <- read.csv("mantel_pcoa.csv",header = T,row.names = 1)
# plant <- as.data.frame(pcoa_data[,1])
# colnames(plant) <- "plant"
 env <- read.csv("env_fu.csv",header = T,row.names = 1)
# # env <- env[-c(1,2,11,15),]
# 
otu <- env


library(rfPermute)
library(ggplot2)
library(tidyverse)
# 
# #install.packages("conflicted")#由于filter() 和 lag()有冲突，所以需要这四行代码
# library(conflicted)
# conflict_prefer("filter", "dplyr")
# conflict_prefer("lag", "dplyr")

#计算显著性,"rfPermute("后面为因变量的名称，如PCoA1gf
otu_rfP <- rfPermute(plant_age~., data = otu, importance = TRUE, ntree = 500, num.rep = 999, num.cores = 1)
otu_rfP
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)
importance_otu.scale

library(tidyverse)
# 添加一个新的变量，这个变量根据 `%IncMSE.pval` 的值来决定显示的星号数量
importance_otu.scale <- importance_otu.scale %>%
  mutate(sig = case_when(
    `%IncMSE.pval` <= 0.001 ~ "***",
    `%IncMSE.pval` < 0.01 ~ "**",
    `%IncMSE.pval` < 0.05 ~ "*",
    TRUE ~ ""
  ))

pq <- ggplot(importance_otu.scale, aes(reorder(rownames(importance_otu.scale),`%IncMSE`),`%IncMSE`)) +   
  geom_col(aes(fill = `%IncMSE` > 0), width = 0.5, color = NA) + 
  geom_text(aes(label=sig), hjust=-0.3, size=5) +
  scale_fill_manual(values = c("#6da9cc", "#d79f93")) + 
  labs(title = NULL, x = NULL, y = 'Increase in MSE (%)', fill = NULL,size = 15) +  
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, color = "black", size = 18),
        axis.text.y = element_text(angle = 0, hjust = 1, color = "black", size = 18),
        axis.line.x = element_line(color = "black", size = 0.5),
        text = element_text(size = 15),
        legend.position = "none") +  
  scale_y_continuous(expand = c(0, 0), limit = c(-3, 20)) +
  geom_hline(yintercept = 0, linetype="solid", color = "black")+
  coord_flip()  
pq
# ggsave("randomforest-alph-bacteria.pdf",device = "pdf",dpi = 600,width = 10,height = 8)





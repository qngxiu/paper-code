
setwd("H:\\东北地理与农业生态研究所\\重金属实验\\修正后数据\\mantel")
rm(list = ls())
plant <- read.csv("mantel_pcoa.csv",header = T,row.names = 1)
plant <- as.data.frame(plant[,1])
colnames(plant) <- "plant_age"
env <- read.csv("env.csv",header = T,row.names = 1)
env <- env[-c(1,2,11,15),]
set.seed(1122)
otu <- cbind(plant,env)


library(rfPermute)
library(ggplot2)
library(tidyverse)

#install.packages("conflicted")#由于filter() 和 lag()有冲突，所以需要这四行代码
# library(conflicted)
# conflict_prefer("filter", "dplyr")
# conflict_prefer("lag", "dplyr")

#计算显著性,"rfPermute("后面为因变量的名称，如PCoA1gf
otu_rfP <- rfPermute(plant_age~., data = otu, importance = TRUE, ntree = 500, num.rep = 9999, num.cores = 1)
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
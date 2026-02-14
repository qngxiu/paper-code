##mantel分析

rm(list=ls())
setwd("C:/Users\\19721\\Desktop\\重金属实验\\修正后数据\\mantel")
data_bio <- read.csv("mantel_bio使用属水平数据.csv",row.names = 1)
data_env <- read.csv("env.csv",row.names = 1)
data_bio <- as.data.frame(t(data_bio))
# 1. 删除 cyanide 列中值为 "ND" 的行
data_env  <- data_env[-c(1:2,11,15),]
data_bio  <- data_bio[-c(1:2,11,15),]
# 2. 将 cyanide 列转换为数值型
data_env $cyanide <- as.numeric(data_env $cyanide)
class(data_bio)

library(ggcor)
library(ggcorrplot)
library(dplyr)
library(vegan)
library(ggplot2)
library(ade4)

# Mantel.test 检验计算矩阵相关性
mantel <- mantel_test(data_bio, data_env, mantel.fun = 'mantel.randtest',spec.dist.method = 'bray', env.dist.method = 'euclidean', 
                      spec.select = list(Nematode = 1:44,
                                         Fugui = 45:984,
                                         Bacteria = 985:2020
                      )) %>% 
  mutate(r_value = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf), 
                       labels = c('<0.25', '0.25-0.5', '>=0.5'), right = FALSE),
         p_value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                       labels = c('<0.001', '0.001-0.01', '0.01-0.05', '>=0.05'), right = FALSE))


quickcor(data_env, type = "upper") +
  geom_square() +
  anno_link(aes(colour = p_value, size = r_value), data = mantel,curvature = -0.2) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
ggsave('ggcor_network_test.pdf',width = 10,height = 7)



quickcor(data_env, type = "upper", show.diag = F) +
  geom_square() +
  anno_link(aes(colour = p_value, size = r_value), data = mantel, curvature = -0.2) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  scale_fill_gradientn(
    # 调整颜色顺序：蓝色(负相关) -> 白色(零相关) -> 红色(正相关)
    colors = c("#4575B4", "#E0F3F8", "#FFFFBF", "#FC8D59", "#D73027"),
    limits = c(-1, 1),
    na.value = "white"
  ) +
  guides(
    size = guide_legend(title = "Mantel's r",
                        override.aes = list(colour = "grey35"), 
                        order = 2),
    colour = guide_legend(title = "Mantel's p", 
                          override.aes = list(size = 3), 
                          order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  )




write.csv(mantel,"物种_mantel.csv")







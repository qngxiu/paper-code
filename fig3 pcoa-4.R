#细菌绘图代码
rm(list = ls())
getwd()
setwd("C:/Users/19721/Desktop/重金属实验/修正后数据/PCoA")
xiufu <- read.csv("ba.csv",header = T,row.names = 1)
merged_df <- xiufu
print(merged_df)
set.seed(123)
library(vegan)
data <- t(merged_df)
data <- decostand(data, method = "hellinger")
group <- read.csv("group.csv",header = T,row.names = 1)
colnames(group) <- "Group"
bray <- vegdist(data, method = 'bray')
bray <- as.matrix(bray)
pcoa <- cmdscale(bray, k = 3, eig = T)
pcoa_data <- data.frame({pcoa$point})
pcoa_data$Sample_ID <- rownames(pcoa_data)
names(pcoa_data)[1:3] <- paste0("PCoA", 1:3)
eig = pcoa$eig
sum_eig <- sum(pcoa$eig)
eig_percent <- round(pcoa$eig/sum_eig*100,1)
poi = pcoa$points
poi = as.data.frame(poi)
pcoa_result <- cbind(pcoa_data,group )
head(pcoa_result)
dune.div <- adonis2(data ~ Group, data = group, permutations = 999, method="bray")
dune.div
dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
dune_adonis
p <- ggplot(pcoa_result, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(aes(color = Group,fill= Group), size = 5,shape = 21) +
  stat_ellipse(data = pcoa_result, 
               geom = "polygon", 
               level = 0.9, 
               linetype = 2, 
               linewidth = 0.5, 
               aes(fill = Group), 
               alpha = 0.3, 
               show.legend = TRUE) +
  labs(
    x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
    y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""),
    title = dune_adonis
  ) +
  scale_colour_manual(values = c("#c0c9e7", "#eedb8c", "#8a69a8", "#cdd6b6")) +
  scale_fill_manual(values = c("#c0c9e7", "#eedb8c", "#8a69a8", "#cdd6b6")) +
  theme(
    legend.position = c(0.4, 0.2),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 18, face = "bold", colour = "black", family = "Times New Roman"),
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 18, face = "bold", family = "Times New Roman"),
    axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", family = "Times New Roman"),
    axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", family = "Times New Roman")
  ) +
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "dashed")
p

#真菌绘图代码

rm(list = ls())
getwd()
setwd("C:/Users/19721/Desktop/重金属实验/修正后数据/PCoA")
xiufu <- read.csv("fu.csv",header = T,row.names = 1)
merged_df <- xiufu
print(merged_df)
set.seed(123)
library(vegan)
data <- t(merged_df)
data <- decostand(data, method = "hellinger")
group <- read.csv("group.csv",header = T,row.names = 1)
colnames(group) <- "Group"
bray <- vegdist(data, method = 'bray')
bray <- as.matrix(bray)
pcoa <- cmdscale(bray, k = 3, eig = T)
pcoa_data <- data.frame({pcoa$point})
pcoa_data$Sample_ID <- rownames(pcoa_data)
names(pcoa_data)[1:3] <- paste0("PCoA", 1:3)
eig = pcoa$eig
sum_eig <- sum(pcoa$eig)
eig_percent <- round(pcoa$eig/sum_eig*100,1)
poi = pcoa$points
poi = as.data.frame(poi)
pcoa_result <- cbind(pcoa_data,group )
head(pcoa_result)
dune.div <- adonis2(data ~ Group, data = group, permutations = 999, method="bray")
dune.div
dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
dune_adonis
p <- ggplot(pcoa_result, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(aes(color = Group,fill= Group), size = 5,shape = 21) +
  stat_ellipse(data = pcoa_result, 
               geom = "polygon", 
               level = 0.9, 
               linetype = 2, 
               linewidth = 0.5, 
               aes(fill = Group), 
               alpha = 0.3, 
               show.legend = TRUE) +
  labs(
    x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
    y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""),
    title = dune_adonis
  ) +
  scale_colour_manual(values = c("#c0c9e7", "#eedb8c", "#8a69a8", "#cdd6b6")) +
  scale_fill_manual(values = c("#c0c9e7", "#eedb8c", "#8a69a8", "#cdd6b6")) +
  theme(
    legend.position = c(0.4, 0.2),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 18, face = "bold", colour = "black", family = "Times New Roman"),
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 18, face = "bold", family = "Times New Roman"),
    axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", family = "Times New Roman"),
    axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", family = "Times New Roman")
  ) +
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "dashed")
p

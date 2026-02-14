#R语言进行多组回归分析拟合和回归显著性计算
#测试数据
#链接：https://pan.baidu.com/s/177hQzQXo8B3kr4JYEzbQfA 提取码：qhs0
#链接：https://pan.baidu.com/s/1NnBquWKGp1cMH2HZXGw-jA 提取码：hayq
#数据按照对应的OTU表,物种文件和分组文件
rm(list = ls())
otu<-read.csv("total_pcoa.csv",row.names = 1)
# design<-read.csv("test_design.csv",header = T,row.names = 1)

# 加载包
library(ggplot2)
library(ggpubr)
library(car)
#构建数据框1，用于回归分析
df1<-otu[1:32,c(1,2,6)]
colnames(df1)<-c("shannon","Env","Group")
# df1$Group[df1$Group>=0]<-"Group_1"
#构建数据框2
df2<-otu[33:64,c(1,2,6)]
colnames(df2)<-c("shannon","Env","Group")
# df2$Group[df2$Group>=0]<-"Group_2"
#构建数据框3
df3<-otu[65:96,c(1,2,6)]
colnames(df3)<-c("shannon","Env","Group")
# df3$Group[df3$Group>=0]<-"Group_3"
#合并三个数据框
df<-rbind(df1,df2,df3)
# df1 <- df1[df1$Env != "ND", ] 
# df1$Env <- as.numeric(df1$Env)
# 
# df2 <- df2[df2$Env != "ND", ] 
# df2$Env <- as.numeric(df2$Env)
# df3 <- df3[df3$Env != "ND", ] 
# df3$Env <- as.numeric(df3$Env)
# df <- df[df$Env != "ND", ]
# df$Env <- as.numeric(df$Env)
fit_lm1 <- lm(shannon ~ Env, data = df1)
fit_lm2 <- lm(shannon ~ Env, data = df2)
fit_lm3 <- lm(shannon ~ Env, data = df3)


class(df$shannon)


# 提取R2和P值的函数（保持不变）
extract_stats <- function(fit) {
  summary_fit <- summary(fit)
  r_squared <- format(round(summary_fit$r.squared, 2), nsmall = 2)
  p_value <- format(round(summary_fit$coefficients[2, 4], 3), nsmall = 3)
  return(list(r_squared = r_squared, p_value = p_value))
}

stats1 <- extract_stats(fit_lm1)
stats2 <- extract_stats(fit_lm2)
stats3 <- extract_stats(fit_lm3)

# 画图（修改注释部分）
p1 <- ggplot(df, aes(x = Env, y = shannon)) +
  geom_point(size = 4, aes(fill = Group), shape = 21, color = "gray2", alpha = 1) +
  geom_smooth(aes(color = Group), method = 'lm', se = T, level = 0.95, size = 1.0) +
  scale_fill_manual(values = c("#F2B379", "#DD5F60", "#9BCD9B")) +
  scale_color_manual(values = c("#F2B379", "#DD5F60", "#9BCD9B")) +
  theme_bw() + labs(x = "elements", y = "PCoA") +
  annotate('text', 
           label = paste("R² =", stats1$r_squared, ", p =", stats1$p_value),
           x = min(df$Env)+3 , 
           y = max(df$shannon) * 0.8+0.15, 
           size = 3.9, 
           color = "#F2B379", 
           hjust = 1) +
  annotate('text', 
           label = paste("R² =", stats2$r_squared, ", p =", stats2$p_value),
           x = min(df$Env)+3, 
           y = max(df$shannon) * 0.8+0.1, 
           size = 3.9, 
           color = "#DD5F60", 
           hjust = 1) +
  annotate('text', 
           label = paste("R² =", stats3$r_squared, ", p =", stats3$p_value),
           x = min(df$Env)+3, 
           y = max(df$shannon) * 0.8+0.05, 
           size = 3.9, 
           color = "#9BCD9B", 
           hjust = 1) +
  theme(axis.text = element_text(colour = 'black', size = 20)) +
  theme(text = element_text(family = "Times New Roman", face = "bold", size = 20),
        axis.text = element_text(colour = 'black', size = 16),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
  # 使用字符串拼接代替bquote
  # annotate('text', 
  #          label = paste("R² =", stats1$r_squared, ", p =", stats1$p_value),
  #          x = max(df$Env) , 
  #          y = min(df$shannon) * 1.4, 
  #          size = 3.9, 
  #          color = "#F2B379", 
  #          hjust = 1) +
  # annotate('text', 
  #          label = paste("R² =", stats2$r_squared, ", p =", stats2$p_value),
  #          x = max(df$Env), 
  #          y = min(df$shannon) * 1.20, 
  #          size = 3.9, 
  #          color = "#DD5F60", 
  #          hjust = 1) +
  # annotate('text', 
  #          label = paste("R² =", stats3$r_squared, ", p =", stats3$p_value),
  #          x = max(df$Env), 
  #          y = min(df$shannon) * 1.0, 
  #          size = 3.9, 
  #          color = "#9BCD9B", 
  #          hjust = 1) +
  # theme(axis.text = element_text(colour = 'black', size = 20)) +
  # theme(text = element_text(family = "Times New Roman", face = "bold", size = 20),
  #       axis.text = element_text(colour = 'black', size = 16),
  #       axis.title = element_text(size = 20),
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank())


p1
















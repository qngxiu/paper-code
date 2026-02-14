# 安装并加载vegan包（首次使用时需要安装）
if (!require(vegan)) install.packages("vegan")
library(vegan)

# 设置工作目录（替换为你的数据所在的目录）
#setwd("H:\\东北地理与农业生态研究所\\重金属投稿\\原始数据--图表--对应\\原始数据--图表--对应\\Table S3")

# 读取数据：你的数据需要是【行为样本，列为物种】的格式，第一列是样本名
# 数据值为物种的丰度（比如OTU丰度、物种数量等，不能有负数）
data <- read.csv("ba.csv", row.names = 1)
data <- as.data.frame(t(data))
# 计算香农多样性指数：margin=1表示按行（每个样本）计算
# index="shannon"指定计算香农指数，默认使用自然对数
shannon_index <- diversity(data, index = "shannon")

# 将结果转换为数据框，方便查看和保存
shannon_result <- data.frame(
  sample_name = names(shannon_index),
  shannon_diversity = shannon_index,
  row.names = NULL
)

# 保存结果到CSV文件
write.csv(shannon_result, "shannon_diversity_result.csv", row.names = FALSE)

# 打印结果
print(shannon_result)
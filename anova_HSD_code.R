# 安装并加载所需R包
if (!require(agricolae)) install.packages("agricolae")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(openxlsx)) install.packages("openxlsx")
if (!require(tibble)) install.packages("tibble")

library(agricolae)
library(dplyr)
library(tidyr)
library(openxlsx)
library(tibble)

# 设置工作目录
setwd("H:\\东北地理与农业生态研究所\\重金属投稿\\原始数据--图表--对应\\原始数据--图表--对应\\Table S3")
rm(list = ls())

# 读取Excel文件：第一行作为列名，第一列作为行名
soil_data <- read.xlsx("1.xlsx", rowNames = TRUE)

# 去掉所有包含缺失值的行和列
data <- na.omit(soil_data)

# 将group列转换为因子类型，确保分组正确
data$group <- as.factor(data$group)

# 确保前19列是数值型，避免方差分析出错
data[, 1:19] <- lapply(data[, 1:19], function(x) as.numeric(as.character(x)))

# 获取前19列的变量名，严格按照数据里第一列到19列的顺序
var_cols <- colnames(data)[1:19]

# 定义你需要的分组顺序：F、G、P、R
group_order <- c("F", "G", "P", "R")

# 创建空列表存储每个变量的结果
results_list <- list()

# 循环处理每个变量，按照数据里的列顺序处理
for (col in var_cols) {
  # 构建方差分析公式：当前变量 ~ group
  formula <- as.formula(paste(col, "~ group"))
  
  # 进行单因素方差分析
  aov_model <- aov(formula, data = data)
  
  # 替换为Tukey HSD事后检验，p.adj="fdr"表示使用FDR调整p值
  tukey_result <- HSD.test(aov_model, "group")
  
  # 正确提取显著性字母：tukey_result$groups的groups列是显著性字母，rownames是分组名
  sig_vector <- tukey_result$groups$groups
  names(sig_vector) <- rownames(tukey_result$groups)
  
  # 将向量转换为数据框
  sig_letters <- sig_vector %>%
    tibble::enframe(name = "group", value = "sig_letter") %>%
    mutate(variable = col) %>%
    select(variable, group, sig_letter)
  
  # 按照指定的F、G、P、R顺序排列分组
  sig_letters <- sig_letters %>%
    mutate(group = factor(group, levels = group_order)) %>%
    arrange(group)
  
  # 将当前变量的结果加入列表
  results_list[[col]] <- sig_letters
}

# 合并所有变量的结果
final_results <- bind_rows(results_list)

# 转换为宽格式：变量为行，严格按照数据里第一列到19列的顺序排列
final_wide <- final_results %>%
  pivot_wider(names_from = group, values_from = sig_letter) %>%
  # 按照数据里的列顺序排列变量，不使用字母排序
  mutate(variable = factor(variable, levels = var_cols)) %>%
  arrange(variable)

# 保存为Excel表格
write.xlsx(final_wide, "anova_tukey_sig_final.xlsx", rowNames = FALSE)

# 打印提示信息
cat("已生成最终的Tukey显著性表格：anova_tukey_sig_final.xlsx\n")
cat("1. 变量顺序严格按照数据里第一列到19列的顺序排列\n")
cat("2. 分组列按照F、G、P、R的顺序排列\n")
cat("3. 每个单元格是对应变量在该分组的显著性字母，完全匹配你的需求\n")
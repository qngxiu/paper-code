rm(list = ls())
setwd("C:\\Users\\19721\\Desktop\\重金属实验\\修正后数据\\venn")
library(VennDiagram)
library(tidyverse)

library(ggvenn)

library(RColorBrewer)
data <- read.csv("nematode.csv",header = T,row.names = 1)
# 提取列名的第一个字母作为分组依据
col_groups <- substr(names(data), 1, 1)

# 获取唯一的分组
unique_groups <- unique(col_groups)
# 初始化一个空的数据框来存储结果，行名与原始数据相同
result_df <- data.frame(matrix(nrow = nrow(data), ncol = length(unique_groups)))
rownames(result_df) <- rownames(data)
colnames(result_df) <- unique_groups

# 对每个分组进行遍历并求和
for (group in unique_groups) {
  group_cols <- names(data)[col_groups == group]
  result_df[, group] <- rowSums(data[, group_cols, drop = FALSE])
}

# 输出结果
result_df

# 创建与result_df相同结构的空数据框
dada <- result_df

# 遍历每一行
for (i in 1:nrow(result_df)) {
  # 获取当前行名
  row_name <- rownames(result_df)[i]
  
  # 将非零值替换为行名，零值设为NA
  dada[i, result_df[i, ] != 0] <- row_name
  dada[i, result_df[i, ] == 0] <- NA
}

# 删除全为NA的行
dada <- dada[rowSums(!is.na(dada)) > 0, ]

# 查看结果
dada

# write.csv(dada,"nematode_venn.csv")


# 
# # 提取每个分组中的元素（非NA值）
# venn_data <- list()
# for (group in unique_groups) {
#   venn_data[[group]] <- unique(dada[!is.na(dada[, group]), group])
# }
# 
# # 绘制四元维恩图
# filename <- "four_set_venn.pdf"
# venn.plot <- venn.diagram(
#   x = venn_data,  # 传入四个分组的元素
#   filename = filename,
#   col = "transparent",  # 背景透明
#   fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),  # 四个集合的填充色
#   alpha = 0.5,  # 透明度
#   label.col = "black",  # 标签颜色
#   cex = 1.2,  # 文本大小
#   fontfamily = "serif",
#   cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),  # 分组标签颜色
#   cat.cex = 1.5,  # 分组标签大小
#   cat.pos = 0,  # 标签位置
#   cat.dist = 0.07,  # 标签与圆的距离
#   rotation.degree = 0,  # 旋转角度
#   margin = 0.2  # 边距
# )


dada_logical <- dada %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ !is.na(.)))
c <- ggvenn(dada_logical, # 数据列表
       
       columns = c("F","G","P","R","O"),
       
       show_percentage = T, # 显示每一组的百分比
       
       digits = 1, # 百分比的小数点位数
       
       fill_color = c("#c0c9e7", "#eedb8c", "#8a69a8", "#cdd6b6","grey"))

c
# write.csv(dada,"veen_fu.csv")


# 筛选所有样本中都存在的物种
universal_species <- dada_logical[rowSums(dada_logical) == ncol(dada_logical), ]

# 如果需要获取物种名称（行名）
universal_species_names <- rownames(universal_species)

# 从原始数据中提取这些物种的数据
universal_data <- data[universal_species_names, ]

# 查看结果
print(paste("共筛选出", nrow(universal_data), "个存在于所有样本中的物种"))
universal_data
 # write.csv(universal_data,"fu_共有_venn.csv")
# tax <- read.csv("tax_ba.csv",row.names = 1)





rm(list = ls())

 universal_data <- read.csv("C:\\Users\\19721\\Desktop\\重金属实验\\修正后数据\\venn/五种韦恩/ba_共有_venn.csv",row.names = 1)

zipi <- read.csv("C:\\Users\\19721\\Desktop\\重金属实验\\修正后数据\\网络图/细菌文件/total_ba_zipi.csv",row.names = 1)
#
zipi1 <- subset(zipi, type != "Peripherals")
# row.names(zipi1) <- zipi1$genus 
common_rownames <- intersect(rownames(universal_data), rownames(zipi1))

keystone <- universal_data[common_rownames, ]
keystone
# 查看结果
print(paste("共找到", nrow(keystone), "个共同的行名"))
head(keystone)


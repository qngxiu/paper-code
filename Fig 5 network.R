#############基于丰度相关性的微生物共发生网络
#准备文件相对丰度表格，分类信息表格
##计算微生物丰度间的相关系数
library(Hmisc)
##获得网络
setwd("C:\\Users\\19721\\Desktop\\重金属实验\\修正后数据\\网络图")

library(igraph)
rm(list = ls())
#以属水平丰度为例，“genus_table.txt” 是一个属水平的微生物丰度表
# genus <- read.delim('genus_table.txt', row.name = 1, check.names = FALSE,skipNul=TRUE)
genus <- read.csv("ba.csv",header = T,row.names = 1)
genus <- as.data.frame((genus))
# neam <- read.csv("nema.csv",header = T,row.names = 1)
# nema <- as.data.frame(((neam)))
# names(genus) <- names(nema)
# 
# genus <- rbind(genus,nema)
# genus <- as.data.frame((genus[,1:8]))
#可选事先过滤一些低丰度或低频的类群
genus <- genus[which(rowSums(genus) >= 0.001), ]    #例如只保留相对丰度总和高于 0.005 的属

##过滤低丰度 OTUs 类群，它们对分类贡献度低，且影响计算效率
# #genus <- genus[which(rowSums(genus) >= 30), ]
# 
# # 
# genus1 <- genus
# genus1[genus1>0] <- 1
# genus <- genus[which(rowSums(genus1) >= 5), ]    #例如只保留在 5 个及以上样本中出现的属
class(genus)
#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_corr <- rcorr(t(genus), type = 'spearman')

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.6
r <- genus_corr$r
r[abs(r) < 0.7] <- 0

#输出相关系数矩阵（相关分析）
# write.table(data.frame(r, check.names = FALSE), 'r.txt', col.names = NA, sep = '\t', quote = FALSE)

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_corr$P

#输出p值矩阵（相关分析）
# write.table(data.frame(p, check.names = FALSE), 'p.txt', col.names = NA, sep = '\t', quote = FALSE)

#p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
#p <- p.adjust(p, method = 'fdr')    #可选 p 值校正，这里使用 FDR 法校正 p 值
# 
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
# write.table(data.frame(z, check.names = FALSE), 'genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)



#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g

#自相关也可以通过该式去除
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
E(g)$color <- ifelse(E(g)$correlation>=0,"pink","grey")
E(g)$cor <- ifelse(E(g)$correlation>=0,1,-1)
#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“taxa.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
#tax <- read.delim('taxonomy.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- read.csv("tax_ba.csv",header = T,row.names = 1)
#将空格转换成unassigned
tax[tax==""] <- "unassigned"

tax <- tax[as.character(V(g)$name), ]

#V(g)$kingdom <- tax$kingdom
V(g)$phylum <- tax$phylum
# V(g)$class <- tax$class
# V(g)$order <- tax$order
# V(g)$family <- tax$family
V(g)$genus <- tax$genus
tax

#计算节点度
V(g)$degree <- degree(g)

#模块划分，详情 ?cluster_fast_greedy，有多种模型
#cfg <- cluster_fast_greedy(g)#划分模块
#plot(cfg,g)#整体网络展示模块
#cfg#显示每一个模块的节点
#nodes <- V(g)[cfg$membership == 1]#提取模块1节点
#g1 <- induced_subgraph(g, nodes)#提取模块1
#plot(g1, layout = layout_in_circle)#绘制模块1

set.seed(123)
V(g)$modularity <- membership(cluster_fast_greedy(g))
#查看网络图
g
plot(g)
write.graph(g, 'total_ba_network.graphml', format = 'graphml')######                                             ##输出网络图




##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')

#重新读取转换的邻接矩阵类型的网络文件
# adjacency_unweight <- read.delim('network_ba.adj_matrix.txt', row.names = 1, sep = '\t', check.names = FALSE)
adjacency_unweight <- as.matrix(as_adjacency_matrix(g, attr = 'correlation'))
#节点属性列表，包含节点所划分的模块
# nodes_list <- read.delim('nodes_ba_list.txt', row.names = 1, sep = '\t', check.names = FALSE)
nodes_list <- data.frame(
  nodes_id = V(g)$name, 
  degree = V(g)$degree, 
  modularity = as.vector(V(g)$modularity)
)
rownames(nodes_list) <- nodes_list$nodes_id
#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

# write.table(zi_pi, 'zi_pi_ba_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

# ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
#   geom_point(aes(color = type), alpha = 0.5, size = 2) +
#   scale_color_manual(values = c('gray','red','blue','purple'), 
#                      limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
#   theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
#         panel.background = element_blank(), legend.key = element_blank()) +
#   labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
#   geom_vline(xintercept = 0.62) +
#   geom_hline(yintercept = 2.5)
# 
# 
# 
# library(ggplot2)
# library(ggrepel)  # 用于避免标签重叠
# 
# ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
#   geom_point(aes(color = type), alpha = 0.5, size = 2) +
#   scale_color_manual(values = c('gray','red','blue','purple'), 
#                      limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs')) +
#   theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
#         panel.background = element_blank(), legend.key = element_blank()) +
#   labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
#   geom_vline(xintercept = 0.62) +
#   geom_hline(yintercept = 2.5) +
#   # 添加标签层，仅显示非灰色节点
#   geom_text_repel(
#     data = subset(zi_pi, type != "Peripherals"),  # 筛选非灰色节点
#     aes(label = nodes_id),  # 使用node_id列作为标签
#     max.overlaps = Inf,    # 显示所有标签
#     box.padding = 0.5,     # 标签与点的间距
#     point.padding = 0.5,   # 标签与其他元素的间距
#     segment.color = "gray50",  # 连接线颜色
#     size = 3               # 标签字体大小
#   )
# 
# 
# 
# 


library(ggplot2)
library(ggrepel)
library(extrafont)

# # 导入系统字体（首次运行需等待几分钟）
# font_import(pattern = "Times New Roman")  # 仅导入新罗马字体
# y
ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs')) +
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(colour = 'black'), 
    panel.background = element_blank(), 
    legend.key = element_blank(),
    # 设置全局字体
    text = element_text(family = "Times New Roman", size = 16, face = "bold", color = "black"),
    axis.title = element_text(family = "Times New Roman", size = 16, face = "bold", color = "black"),
    axis.text = element_text(family = "Times New Roman", size = 16, face = "bold", color = "black"),
    legend.text = element_text(family = "Times New Roman", size = 16, face = "bold", color = "black"),
    legend.title = element_text(family = "Times New Roman", size = 16, face = "bold", color = "black"),
    plot.title = element_text(family = "Times New Roman", size = 16, face = "bold", color = "black")
  ) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5) +
  # 修改标签字体
  geom_text_repel(
    data = subset(zi_pi, type != "Peripherals"),
    aes(label = nodes_id),
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "gray50",
    # 设置标签字体为斜体、16号、新罗马
    family = "Times New Roman",
    size = 4,  # ggplot2中size单位为mm，16pt约等于16/2.835mm
    fontface = "bold.italic"
  )


write.csv(zi_pi,"total_ba_zipi.csv")

####计算拓扑性质


#加载包
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(igraph)
library(RColorBrewer)
library(reshape2)

set.seed(123)

ZiPiPlot1  <-  function(igraph = igraph,method = "cluster_fast_greedy")  #z值和p值,展开之后修改,修改完成之后再点击最左边行数的折叠箭头选择折叠
{
  
  if (method == "cluster_walktrap" ) {
    fc <- igraph::cluster_walktrap(igraph,weights =  abs(igraph::E(igraph)$weight))
  }
  if (method == "cluster_edge_betweenness" ) {
    fc <- igraph::cluster_edge_betweenness(igraph,weights =  abs(igraph::E(igraph)$weight))
  }
  if (method == "cluster_fast_greedy" ) {
    fc <- igraph::cluster_fast_greedy(igraph,weights =  abs(igraph::E(igraph)$weight))
  }
  if (method == "cluster_spinglass" ) {
    fc <- igraph::cluster_spinglass(igraph,weights =  abs(igraph::E(igraph)$weight))
  }
  
  modularity <- igraph::modularity(igraph,igraph::membership(fc))
  # 模块化程度
  # 按照模块为节点配色，这里我们可以加入到nodes中
  comps <- igraph::membership(fc)
  
  # comps
  igraph::V(igraph)$module <- as.character(comps)
  
  taxa.roles <- module.roles(igraph)
  
  taxa.roles$label = row.names(taxa.roles)
  for (i in 1:nrow(taxa.roles))if(taxa.roles[i,3]> 0.62|taxa.roles[i,1]> 1.0) {
    taxa.roles[i,5]=taxa.roles[i,5]
  }else{
    taxa.roles[i,5]= ""
  }
  taxa.roles$role_7 = taxa.roles$roles
  
  taxa.roles <- na.omit(taxa.roles)   # remove NA values
  taxa.roles[which(taxa.roles$z < 2.5 & taxa.roles$p < 0.62),'roles'] <- 'Peripherals'
  taxa.roles[which(taxa.roles$z < 2.5 & taxa.roles$p >= 0.62),'roles'] <- 'Connectors'
  taxa.roles[which(taxa.roles$z >= 2.5 & taxa.roles$p < 0.62),'roles'] <- 'Module hubs'
  taxa.roles[which(taxa.roles$z >= 2.5 & taxa.roles$p >= 0.62),'roles'] <- 'Network hubs'
  
  
  
  p <- plot_roles2(taxa.roles) +
    ggrepel::geom_text_repel(data = taxa.roles,
                             aes(x = p, y = z, color = module,label=taxa.roles$label),size=4)#
  p
  #geom_text(data = taxa.roles, aes(x = p, y = z, color = module,label=taxa.roles$label),size=4)
  # print(p)
  
  return(list(p,taxa.roles))
}

net_properties.3 <- function (igraph, n.hub = FALSE)    #n.hub设置为TRUE或者FALSE,T的话则计算关键节点,F的话不计算关键节点.
{
  num.edges <- length(igraph::E(igraph))
  num.edges
  num.vertices <- length(igraph::V(igraph))
  num.vertices
  connectance <- igraph::edge_density(igraph, loops = FALSE)
  average.degree <- mean(igraph::degree(igraph))
  average.degree
  if (!is.null(igraph::E(igraph)$weight)) {
    igraph.weight <- igraph::E(igraph)$weight
    igraph::E(igraph)$weight = abs(igraph::E(igraph)$weight)
  }
  average.path.length <- igraph::average.path.length(igraph)
  average.path.length
  diameter <- igraph::diameter(igraph, directed = FALSE, unconnected = TRUE, 
                               weights = NULL)
  diameter
  if (!is.null(igraph::E(igraph)$weight)) {
    igraph::E(igraph)$weight = igraph.weight
  }
  edge.connectivity <- igraph::edge_connectivity(igraph)
  edge.connectivity
  clustering.coefficient <- igraph::transitivity(igraph, type = "average")
  clustering.coefficient
  no.clusters <- igraph::no.clusters(igraph)
  no.clusters
  centralization.degree <- igraph::centralization.degree(igraph)$centralization
  centralization.degree
  centralization.betweenness <- igraph::centralization.betweenness(igraph)$centralization
  centralization.betweenness
  centralization.closeness <- igraph::centralization.closeness(igraph)$centralization
  centralization.closeness
  if (!is.null(igraph::E(igraph)$weight)) {
    num.pos.edges <- sum(igraph.weight > 0)
    num.neg.edges <- sum(igraph.weight < 0)
  }
  else {
    num.pos.edges <- 0
    num.neg.edges <- 0
  }
  modularity_igraph = function(net, method = "cluster_walktrap") {
    if (method == "cluster_walktrap") {
      fc <- igraph::cluster_walktrap(net, weights = abs(igraph::E(igraph)$weight))
    }
    if (method == "cluster_edge_betweenness") {
      fc <- igraph::cluster_edge_betweenness(net, weights = abs(igraph::E(igraph)$weight))
    }
    if (method == "cluster_fast_greedy") {
      fc <- igraph::cluster_fast_greedy(net, weights = abs(igraph::E(igraph)$weight))
    }
    if (method == "cluster_spinglass") {
      fc <- igraph::cluster_spinglass(net, weights = abs(igraph::E(igraph)$weight))
    }
    modularity <- igraph::modularity(net, membership(fc))
    return(modularity)
  }
  mod1 = modularity_igraph(igraph, method = "cluster_walktrap")
  rand.g <- igraph::erdos.renyi.game(length(V(igraph)), length(E(igraph)), 
                                     type = "gnm")
  mod2 = modularity_igraph(rand.g, method = "cluster_walktrap")
  RM = (mod1 - mod2)/mod2
  if (n.hub) {
    res = ZiPiPlot1(igraph = tem.g, method = "cluster_walktrap")
    data = res[[2]]
    head(data)
    n.hub = data$roles[data$roles != "Peripherals"] %>% length()
  }
  else {
    n.hub = "Not.calculated"
  }
  igraph.network.pro <- rbind(num.edges, num.pos.edges, num.neg.edges, 
                              num.vertices, connectance, average.degree, average.path.length, 
                              diameter, edge.connectivity, clustering.coefficient, 
                              no.clusters, centralization.degree, centralization.betweenness, 
                              centralization.closeness, RM, n.hub)
  rownames(igraph.network.pro) <- c("num.edges(L)", "num.pos.edges", 
                                    "num.neg.edges", "num.vertices(n)", "Connectance(edge_density)", 
                                    "average.degree(Average K)", "average.path.length", "diameter", 
                                    "edge.connectivity", "mean.clustering.coefficient(Average.CC)", 
                                    "no.clusters", "centralization.degree", "centralization.betweenness", 
                                    "centralization.closeness", "RM(relative.modularity)", 
                                    "the.number.of.keystone.nodes")
  colnames(igraph.network.pro) <- "value"
  return(igraph.network.pro)
}


##读取数据

otutab <- as.data.frame((genus))
tax <- read.csv("tax_ba.csv",row.names = 1)
sam <- read.csv("sam.csv",row.names = 1)
#以上数据读取自己需要分析的数据，其中tax数据非必要
ps = phyloseq(sample_data(sam),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(tax))
)

igraph = make_igraph(z)

dat = igraph::V(igraph)
names(dat) %>% length()
#--弄清楚每个样本包含的OTU数量
# pst =  ps %>%
#   scale_micro("rela") %>%
#   phyloseq::subset_samples(Group %in% c("KO","WT","OE")) %>%
#   filter_OTU_ps(500) 
otu = ps %>% subset_taxa(row.names(tax_table(ps)) %in% names(dat)) %>% vegan_otu() %>% t()

# otu = ps %>% 
#   phyloseq::subset_samples(Group %in% c("KO","WT","OE")) %>%
#   # filter_OTU_ps(500) %>%
#   subset_taxa(row.names(tax_table(ps)) %in% names(dat)) %>%
#   vegan_otu() %>% 
#   t() 
dim(otu)

otu[otu > 1] = 1
dim(otu)
A = list()
i = 1
for (i in 1:length(colnames(otu))) {
  tem = otu[,colnames(otu)[i]][otu[,colnames(otu)[i]] > 0 ] %>% names()
  A[[colnames(otu)[i]]] = tem
  #-计算性质
  tem.2 = A[[colnames(otu)[i]]]
  tem.g = igraph::induced_subgraph(igraph,tem.2)
  dat = net_properties.3(tem.g,n.hub = T)   ##是否计算关键节点
  head(dat,n = 16)
  
  dat = as.data.frame(dat)
  dat$value = as.numeric(dat$value)
  colnames(dat) = colnames(otu)[i]
  if (i == 1) {
    dat.f = dat
  } else {
    dat.f = cbind(dat.f,dat)
  }
}

dat.f 
write.csv(dat.f,"total_ba_拓扑.csv")
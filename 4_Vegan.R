# install.packages("vegan")
library(stringr)
library(devtools)
# install_github("microbiota/amplicon")
# install_github("vqv/ggbiplot")
# library(ggbiplot) # M1/R 版本用不了这个包
library(ade4)

# library(picante)	#用于计算 PD_whole_tree，若不计算它就无需加载。
library(ggplot2)	#用于 ggplot2 作图
library(ggsignif)
library(doBy) 	#用于分组统计
library(ggalt)	#用于绘制拟合曲线
library(ggpubr)
library(dplyr)
library(scales)
library(ggsci)
library(RColorBrewer)
#计算多种 Alpha 多样性指数，结果返回至向量
library(data.table)
library(tibble)
library(vegan)	#用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样

outlier_process <- function(x){
  temp <- apply(mat, 1, function(x){
    q <- quantile(x)
    iqr <- q[4]-q[2]
    median <- median(x)
    x[x > as.numeric(q[4] + 3 * iqr) | x < as.numeric(q[2] - 3 * iqr)] <- median # 1.5 or 3
    return(x)
  })
  return(t(temp))
}

alpha_index <- function(x,  tree = NULL, base = exp(1)) {
  result <- data.frame(richness = rowSums(x > 0), chao1 = estimateR(x)[3, ], ace = estimateR(x)[5, ],
                       shannon = diversity(x, index = 'shannon'), simpson = diversity(x, index = 'simpson'),
                       pielou = diversity(x, index = 'shannon') / log(estimateR(x)[1, ], base), 
                       goods_coverage = 1 - rowSums(x == 1) / rowSums(x))
  if(!is.null(tree)) {#PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- cbind(result, pd[ ,1])
  }
  return(result)
}

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")
# 读取输入文件
mat <- read.table("species.txt")
groupinfo <- read.table("meta.txt")
groupinfo <- groupinfo[groupinfo$discover==1,]
mat <- mat[,rownames(groupinfo)]

x <- round(t(mat))
alpha.index <- data.frame(richness = rowSums(x > 0), chao1 = estimateR(x)[3, ], ace = estimateR(x)[5, ],
                          shannon = diversity(x, index = 'shannon'), simpson = diversity(x, index = 'simpson'),
                          pielou = diversity(x, index = 'shannon') / log(estimateR(x)[1, ], exp(1)), 
                          goods_coverage = 1 - rowSums(x == 1) / rowSums(x),
                          group = groupinfo$Group, cancer = groupinfo$Cancer,
                          multigroup = groupinfo$multigroup)

table(rownames(alpha.index) == colnames(mat))
# alpha.index <- outlier_process(alpha.index)
manual_color_vector = c("#009FFD", "#FFA400", "#D00000")

gginput <- reshape2::melt(alpha.index, id.vars = c("group", "cancer", "multigroup"), measure.vars = c("richness", "chao1", "ace", "shannon", "simpson", "pielou", "goods_coverage"))
gginput$group <- factor(gginput$group, levels = c("Normal", "Overweight", "Obesity"))
gginput$multigroup <- factor(groupinfo$multigroup, levels = c("HN", "CN", "HOv", "COv", "HOb", "COb"))

my_comparison <- combn(as.character(unique(groupinfo$multigroup)), 2, simplify=FALSE)
unique(gginput$variable)
# my_comparison <- my_comparison[c(1,2,3,6,8,12,13,14,15)]
plotlist <- list()
wiltestpvalue <- NA

for(i in 1:length(unique(gginput$variable))){
  # i = 4
  test <- gginput[gginput$variable==unique(gginput$variable)[i],]
  idj <- NA
  for(j in 1:length(my_comparison)) {
    wiltestpvalue[j] <- t.test(test[test$multigroup == my_comparison[[j]][1],]$value, test[test$multigroup==my_comparison[[j]][2],]$value)$p.value
  }

  # wiltestpvalue <- p.adjust(wiltestpvalue, method = "bonferroni")
  idj <- ifelse(wiltestpvalue < 0.05, 1, 0)
  
  plotlist[[i]] <- ggplot(test,aes(x=multigroup, y=value, fill=group))+
    geom_violin(trim=FALSE, aes(linetype=NA)) +
    geom_boxplot(width = 0.25, outlier.size = 0.25) +
    # geom_point(position = position_jitterdodge(),size=0.3)+
    stat_compare_means(comparisons =my_comparison[idj==1],
                       method = "t.test",
                       label = "p.signif",
                       hide.ns = TRUE)+
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5))+
    labs(title = unique(gginput$variable)[i]) +
    scale_fill_manual(values = manual_color_vector)
}

length(plotlist)
filename <- "Alpha_index_boxplot_unadj.pdf"
pdf(filename, width = 8, height = 8)
print(ggarrange(plotlist = plotlist[1:7], nrow = 2, ncol = 4))
dev.off()


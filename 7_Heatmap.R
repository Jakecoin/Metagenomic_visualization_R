library(lattice)
library(Formula)
library(readxl)
library(ggpubr)
library(grid)
library(vcd)
library(tibble)
library(RColorBrewer)

library(caret)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(devtools)

library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(convertid)
library(ComplexHeatmap)
library(circlize)
# library(pheatmap)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")
getwd()

groupinfo <- group_input()
mat <- metagenome_input()
mat <- metabolite_input()

mat <- GO_KEGG_input(i = 6)
mat <- mat[-1,]
mat <- as.data.frame(t(apply(mat, 1, function(x){scale(x)})))
# colnames(mat) <- colnames(mat)

heatmap_out <- data.frame(item = rownames(mat),
                          H_N = apply(mat[, groupinfo$multigroup=="H_N"], 1, function(x){mean(x)}),
                          C_N = apply(mat[, groupinfo$multigroup=="C_N"], 1, function(x){mean(x)}),
                          H_Ov = apply(mat[, groupinfo$multigroup=="H_Ov"], 1, function(x){mean(x)}),
                          C_Ov = apply(mat[, groupinfo$multigroup=="C_Ov"], 1, function(x){mean(x)}),
                          H_Ob = apply(mat[, groupinfo$multigroup=="H_Ob"], 1, function(x){mean(x)}),
                          C_Ob = apply(mat[, groupinfo$multigroup=="C_Ob"], 1, function(x){mean(x)})
                          )
heatmap_out[1:5,1:5]

rownames(heatmap_out) <- gsub(heatmap_out$item, pattern = "s__", replacement = "")

heatmap_out <- read.table("/Users/lijinming/Downloads/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt",
                          header = T)
# simple heatmap
{
  pdf("heatmap.pdf")

  Heatmap(heatmap_out[,-c(1:2)], name = "Z-score", show_column_names = TRUE, show_row_names = TRUE,
          column_names_side = "top", show_column_dend = TRUE, show_row_dend = TRUE,
          row_names_gp = gpar(fontsize = 10),
          column_names_rot = 0, column_names_centered = TRUE,
          column_order = levels(groupinfo$multigroup),
          # top_annotation = HeatmapAnnotation(BMI = anno_block(gp = gpar(fill = c("#0F8B8D", "#EC9A29", "#A8201A")),
          #                                                     labels = c("Normal", "Overweight", "Obesity"),
          #                                                     labels_gp = gpar(col = "white", fontsize = 10))),
          # column_split = c(1,1,2,2,3,3),
          column_title = ""
  )

  dev.off()
  #column_order = order(pro_taxgroup$group)
  #heatmap_out[,order(pro_taxgroup$Group)]
}

signdiff <- data.frame(item = rownames(mat),
                       H_N = rep(NA,dim(heatmap_out)[1]), C_N = rep(NA,dim(heatmap_out)[1]), 
                       H_Ov = rep(NA,dim(heatmap_out)[1]), C_Ov = rep(NA,dim(heatmap_out)[1]),
                       H_Ob = rep(NA,dim(heatmap_out)[1]), C_Ob = rep(NA,dim(heatmap_out)[1]))
head(signdiff)

diff <- read.csv("Species_diff_pvalue_padjust.csv")
# diff <- read.csv("met_diff_pvalue_padjust.csv")
diff <- read.csv("KOEntry_diff.csv") # "GO_Function_diff.csv", "GO_ID_diff.csv", "KOEntry_diff.csv"

# 确认行名称匹配
diff <- diff[match(heatmap_out$item, diff$item),]
head(diff[,1:5])

diff[1:5,27:29]
diff[1:5,42:44]

# temp <- signdiff
# signdiff <- temp

for(i in 1:dim(diff)[1]) {
  for(j in 27:29){ # log2fc column
    k <- j+15 # fdr column 42:44
    # k <- j+6 # p column 33:35
    if(!is.na(diff[i,j]) & !is.na(diff[i,k])){ # +15 = fdr column
      psign <- diff[i,j]
      psign_num <- NA
      
      if(diff[i,k] < 0.0001) {
        psign_num <- 4
      }else if(diff[i,k] < 0.001) {
        psign_num <- 3
      }else if(diff[i,k] < 0.01) {
        psign_num <- 2
      }else if(diff[i,k] < 0.05) {
        psign_num <- 1
      }
      
      if(!is.na(psign) & !is.na(psign_num)){
        signdiff[i, ((j-26)*2+1)] <- paste(rep(psign, psign_num), collapse ="") # 3, 5, 7
      }
    }
  }
}

#总体有差异
idx <- (diff$NC_H_fdr < 0.05 | diff$OvC_H_fdr < 0.05 | diff$ObC_H_fdr < 0.05)
table(idx)

# idx <- !is.na(diff)
# table(idx)
# diff <- na.omit(diff)
# signdiff <- signdiff[diff$item,]
# heatmap_out <- heatmap_out[idx,]

# species #######################################################
# signature 1
# 与健康人相比，仅肥胖肠癌人群有差异，而正常体重和超重无差异的菌群
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr > 0.05 & diff$ObC_H_fdr < 0.05 & abs(diff$ObC_H_log2fc) > 1)
idx[is.na(idx)] <- FALSE
table(idx)

# signature 2
# 与健康人相比，仅正常体重肠癌人群无差异，而肥胖和超重有差异的菌群
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05 & 
        abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1)
idx[is.na(idx)] <- FALSE
table(idx)

# signature 3
# 与健康人相比，仅肥胖肠癌人群无差异，而正常体重和超重有差异的菌群
idx <- (diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr > 0.05 &
          abs(diff$NC_H_log2fc) > 2 & abs(diff$OvC_H_log2fc) > 2)
idx[is.na(idx)] <- FALSE
table(idx)

# metabolite #######################################################
# signature 1 不存在，改变策略，只看超重，但纳入肥胖有差异的代谢物
# 与健康人相比，仅肥胖肠癌人群有差异，而正常体重和超重无差异的
idx <- (diff$NC_H_p > 0.05 & diff$OvC_H_p > 0.05 & diff$ObC_H_p < 0.05 & abs(diff$ObC_H_log2fc) > 1)
table(idx)

# signature 2
# 与健康人相比，仅正常体重肠癌人群无差异，而超重有差异的菌群
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr < 0.05 & abs(diff$OvC_H_log2fc) > 1) | diff$ObC_H_fdr < 0.05

# pvalue
idx <- (diff$NC_H_p > 0.05 & diff$OvC_H_p < 0.05 & diff$ObC_H_p < 0.05) &
           (abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1)
table(idx)

# signature 3
# 与健康人相比，仅肥胖肠癌人群无差异，而正常体重和超重有差异的菌群
idx <- (diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr > 0.05) &
        (abs(diff$NC_H_log2fc) > 1 | abs(diff$OvC_H_log2fc) > 1)
table(idx)
# pvalue
idx <- idx <- (diff$NC_H_p < 0.05 & diff$OvC_H_p < 0.05 & diff$ObC_H_p > 0.05 &
                 abs(diff$NC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1)
table(idx)

# signature 4
idx <- (diff$COb_Ov_p < 0.05 | diff$COb_N_p < 0.05 | diff$COv_N_p < 0.05) & 
  (diff$HOb_Ov_p > 0.05 | diff$HOb_N_p > 0.05 | diff$HOv_N_p > 0.05)
table(idx)

#GO_ID######################
# signature 1
# 与健康人相比，仅肥胖肠癌人群有差异，而正常体重和超重无差异的菌群
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr > 0.05 & diff$ObC_H_fdr < 0.05)
table(idx)

# signature 2
# 与健康人相比，仅正常体重肠癌人群无差异，而肥胖和超重有差异的菌群
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05 & 
          abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1)
table(idx)

# signature 3
# 与健康人相比，仅肥胖肠癌人群无差异，而正常体重和超重有差异的菌群
idx <- (diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr > 0.05 &
          abs(diff$NC_H_log2fc) > 2 & abs(diff$OvC_H_log2fc) > 2)
table(idx)

# KO_ENTRY
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr > 0.05 & diff$ObC_H_fdr < 0.05 & abs(diff$ObC_H_log2fc) > 1) | 
       (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05 & 
        abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1) | 
       (diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr > 0.05 &
       abs(diff$NC_H_log2fc) > 2 & abs(diff$OvC_H_log2fc) > 2) |
  (diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05)
  
idx <- (diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05) & (abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1 & abs(diff$NC_H_log2fc) > 1)

idx[is.na(idx)] <- FALSE
table(idx)

#
col_fun = colorRamp2(c(-1, 0, 1), c("#0F8B8D", "white", "#A8201A"))

# 开始
hout <- heatmap_out[idx,-1]
hsign <- signdiff[idx,-1]

# Relative abundance
pdf("KO_ENTRY.pdf", height = 1.5, width = 3)

print(
  Heatmap(hout, name = "Z-score", col = col_fun,
          show_column_names = TRUE,# show_row_names = FALSE,
          column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
          row_names_gp = gpar(fontsize = 7),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 0, column_names_centered = TRUE,
          column_order = levels(groupinfo$multigroup), 
          top_annotation = HeatmapAnnotation(BMI = anno_block(gp = gpar(fill = c("#0F8B8D", "#EC9A29", "#A8201A")),
                                                              labels = c("Normal", "Overweight", "Obesity"),
                                                              labels_gp = gpar(col = "white", fontsize = 10))),
          column_split = c(1,1,2,2,3,3),
          column_title = "",
          # width = unit(10, "cm"), height = unit(40, "cm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(!is.na(hsign[i,j]) == TRUE) {
              grid.text(hsign[i,j], x, y, gp = gpar(fontsize = 5))
            }
          }
  )
)

dev.off()

# log2fc
diff <- read.csv("Species_diff_pvalue_padjust.csv")
# diff <- read.csv("met_diff_pvalue_padjust.csv")
# diff <- read.csv("KOEntry_diff.csv")

rownames(diff) <- diff[,1]
# diff <- diff[,c(28:30, 40:42)]

diff <- diff[is.finite(diff$NC_H_log2fc) & is.finite(diff$OvC_H_log2fc) & is.finite(diff$ObC_H_log2fc),]

# species #######################################################
# signature 1
# 与健康人相比，仅肥胖肠癌人群有差异，而正常体重和超重无差异的菌群
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr > 0.05 & diff$ObC_H_fdr < 0.05 & abs(diff$ObC_H_log2fc) > 1)
idx[is.na(idx)] <- FALSE
table(idx)

# signature 2
# 与健康人相比，仅正常体重肠癌人群无差异，而肥胖和超重有差异的菌群
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05 & 
          abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1)
idx[is.na(idx)] <- FALSE
table(idx)

# signature 3
# 与健康人相比，仅肥胖肠癌人群无差异，而正常体重和超重有差异的菌群
idx <- (diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr > 0.05 &
          abs(diff$NC_H_log2fc) > 2 & abs(diff$OvC_H_log2fc) > 2)
idx[is.na(idx)] <- FALSE
table(idx)

diff <- diff[idx,]
dim(diff)
signdiff <- data.frame(item = rownames(diff),
                       N = NA,
                       Ov = NA,
                       Ob = NA)
rownames(signdiff) <- rownames(diff)
signdiff <- signdiff[,-1]
head(signdiff)

colnames(diff)

summary(diff$NC_H_log2fc)
summary(diff$OvC_H_log2fc)
summary(diff$ObC_H_log2fc)

col2 <- colorRamp2(c(-floor(max(diff)+1), 0, floor(max(diff)+1)), c("#0F8B8D", "white", "#A8201A"))

for(i in 1:dim(diff)[1]) {
  for(j in 30:32){ # log2fc column
    k <- j+12 # fdr column 42:44
    # k <- j+3 # p column 33:35
    if(!is.na(diff[i,j]) & !is.na(diff[i,k])){
      psign <- diff[i, j-3]
      psign_num <- NA
      
      if(diff[i,k] < 0.0001) {
        psign_num <- 4
      }else if(diff[i,k] < 0.001) {
        psign_num <- 3
      }else if(diff[i,k] < 0.01) {
        psign_num <- 2
      }else if(diff[i,k] < 0.05) {
        psign_num <- 1
      }
      
      if(!is.na(psign) & !is.na(psign_num)){
        signdiff[i, (j-28)] <- paste(rep(psign, psign_num), collapse ="")
      }
    }
  }
}

hout <- diff[,c(30:32)]
rownames(hout) <- diff$item
colnames(hout) <- c("N", "Ov", "Ob")
hsign <- signdiff[,-1]

# match(rownames(groupinfo), diff$item)

print(
  Heatmap(hout, name = "Z-score", col = col2,
          show_column_names = TRUE,# show_row_names = FALSE,
          column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
          row_names_gp = gpar(fontsize = 7),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 0, column_names_centered = TRUE,
          column_order = c("N", "Ov", "Ob"), 
          # width = unit(10, "cm"), height = unit(40, "cm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(!is.na(hsign[i,j]) == TRUE) {
              grid.text(hsign[i,j], x, y, gp = gpar(fontsize = 5))
            }
          }
  )
)

#########circle################

diff <- read.csv("Species_diff_pvalue_padjust.csv")
diff <- read.csv("met_diff_pvalue_padjust.csv")
diff <- read.csv("KOEntry_diff.csv")

rownames(diff) <- diff[,1]
diff <- as.data.frame(na.omit(diff[,-c(1,2)]))
diff <- diff[,c(28:30, 40:42)]

diff <- diff[is.finite(diff[,1]) & is.finite(diff[,2]) & is.finite(diff[,3]),]

# for(i in 1:nrow(diff)){
#   for(j in 1:3){
#     if(diff[i,j+3]>0.05){
#       diff[i,j] <- 0
#     }
#   }
# }

# species #######################################################
idx <- (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr > 0.05 & diff$ObC_H_fdr < 0.05 & abs(diff$ObC_H_log2fc) > 1) | 
       (diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05 & 
          abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1) |
       (diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr > 0.05 & 
          abs(diff$NC_H_log2fc) > 2 & abs(diff$OvC_H_log2fc) > 2)
idx[is.na(idx)] <- FALSE
table(idx)

diff <- diff[idx,]

# column_od <- hclust(dist(t(diff)))$order

col2 <- colorRamp2(c(-floor(max(diff)+1), 0, floor(max(diff)+1)), c("#0F8B8D", "white", "#A8201A"))

circos.par(gap.after = c(10))

split = NA
split[(diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr > 0.05 & diff$ObC_H_fdr < 0.05 & abs(diff$ObC_H_log2fc) > 1)] = 1
split[(diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05 & 
         abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1)] = 2
split[(diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr > 0.05 & 
         abs(diff$NC_H_log2fc) > 2 & abs(diff$OvC_H_log2fc) > 2)] = 3
split = factor(split, levels = 1:3)

circos.heatmap(diff[,1:3], col = col2, split = split,
               cluster = TRUE, # dend.side = "inside",
               rownames.side = "outside", rownames.cex = 0.2,
               track.height = 0.1)

circos.track(
  track.index = get.current.track.index(), 
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 1) {
      # 在最后一个扇形中
      cn = c("N", "Ov", "Ob")
      n = length(cn)
      circos.text(
        rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
        1:n - 0.5, cn, cex = 0.5, adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }, 
  bg.border = NA
)

lgd <- Legend(title = "log2FC", col_fun = col2)
grid.draw(lgd)

circos.clear()

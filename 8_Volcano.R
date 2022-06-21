library(lattice);library(Formula)
library(readxl)
library(ggpubr)
library(grid)
library(vcd)
library(data.table)
library(RColorBrewer)

library(caret)
require(lattice)
require(Formula)
require(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
# library(devtools)
# devtools::install_github("jaspershen/MetNormalizer")
# library(MetNormalizer)

library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(convertid)


meta <- read.table("meta.txt")
mat <- read.table("met.txt")
mat <- mat[,meta$discover==1]
metname <- read.csv("processed_data/metabolite_info.csv")
# mat$name <- rownames(mat)
mat[1:5,1:5]

# s.mat <- t(apply(mat, 1, function(x){scale(x)}))
# colnames(s.mat) <- colnames(mat)
s.mat <- mat

index = 2
target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")
group.test <- c("Normal", "Overweight", "Obesity")

i = 3

diff <- read.table(paste("Result/", save_path[index], "/", group.test[i], "_maaslin2.txt", sep = ""), 
                    header = T, row.names = 1)
VIP <- read.table(paste("VIP_", group.test[i], ".txt", sep = ""))

diff_output <- data.frame(
  name = metname$MS2Metabolite[match(rownames(s.mat), metname$MS2kegg)],
  concentration = apply(s.mat, 1, median),
  pvalue = apply(s.mat,1,function(x){wilcox.test(unlist(x[meta$Cancer == "CRC" & meta$Group == group.test[i]]),unlist(x[meta$Cancer == "Health" & meta$Group == group.test[i]]))$p.value}),
  log2fc = apply(s.mat,1,function(x){log2(median(na.omit(x[meta$Cancer == "CRC" & meta$Group == group.test[i]]))/(median(na.omit(x[meta$Cancer == "Health" & meta$Group == group.test[i]]))))}),
  group = rep(group.test[i], nrow(s.mat))
  )  %>% arrange(pvalue) %>% mutate(fdr=p.adjust(pvalue,method = "BH")) %>%
  mutate(label=ifelse(log2fc>threshold,"UP",ifelse(log2fc<(-threshold),"DOWN","Nosig"))) %>%
  mutate(TPplotlabel = ifelse((fdr<0.05 & abs(log2fc)>threshold), name, NA)) %>%
  mutate(TPplotcolor = as.factor(paste(label,ifelse(is.na(TPplotlabel),0,1),sep = "")))

# diff_output$pval.maaslin <- obdiff$pval[match(rownames(diff_output), diff$feature)]
# diff_output$qval.maaslin <- obdiff$qval[match(rownames(diff_output), diff$feature)]
diff_output$VIP <- VIP[,2][match(rownames(diff_output), VIP[,1])]

diff_output <- diff_output %>%
  mutate(TPplotlabel = ifelse((pvalue < 0.05 & abs(log2fc)>threshold & VIP > 1.5), name, NA)) %>%
  mutate(TPplotcolor = as.factor(paste(label,ifelse(is.na(TPplotlabel),0,1),sep = "")))

diff_output <- diff_output[!duplicated(diff_output$TPplotlabel) | is.na(diff_output$TPplotlabel),]

# write.csv(diff_output,"Ob_met_diff_CRCvsH.csv")
x.axis.lim <- round(max(abs(diff_output$log2fc))+1)

# pdf(paste("vocano_met_", group.test[i], ".pdf", sep = ""))
# ggplot(data=as.data.frame(diff_output), aes(x=log2fc, y=-log10(pvalue))) +
#   geom_point(aes(color = TPplotcolor, size = log10(concentration)), alpha = 0.8) +
#   scale_color_manual(labels = c("Nosig0"="No Significance","UP1"="Significantly Up","UP0"="Up","DOWN1"="Significantly Down","DOWN0"="Down"),
#                      values = c("UP0"='#FAC0AE',"DOWN0"='#9BCFF0',"UP1"='#FA2311',"DOWN1"='#6175DB',"Nosig0" ='gray')) +
#   labs(x="log2 Fold Change",  y="-log10 pvalue") +
#   theme(panel.grid = element_blank(),
#         panel.background = element_rect(color = 'black', fill = 'transparent'), 
#         legend.position = c(0.12, 0.75)) +
#   theme(legend.title = element_blank(), 
#         legend.key = element_rect(fill = 'transparent'),
#         legend.background = element_rect(fill = 'transparent')) +
#   geom_hline(yintercept = -log10(0.05),linetype=4, color = 'black', size = 0.5) +
#   geom_vline(xintercept = c(-threshold, threshold), linetype = 4, color = 'black', size = 0.5) +
#   geom_text_repel(data=diff_output, aes(x = log2fc, y = -log10(pvalue), label = TPplotlabel), size=5, colour = "black") +
#   xlim(-x.axis.lim, x.axis.lim)+ ##控制横坐标长度
#   theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
# dev.off()

# ndiff <- read.table(paste("Result/", save_path[index], "/", group.test[1], "_maaslin2.txt", sep = ""), 
#                     header = T, row.names = 1)
# ovdiff <- read.table(paste("Result/", save_path[index], "/", group.test[2], "_maaslin2.txt", sep = ""), 
#                      header = T, row.names = 1)
# obdiff <- read.table(paste("Result/", save_path[index], "/", group.test[3], "_maaslin2.txt", sep = ""), 
#                      header = T, row.names = 1)
# 
# num <- 6 # 6 pvalue / 8 fdr
# threshold = 0.05
# table(ndiff[,num] < threshold)
# table(ovdiff[,num] < threshold)
# table(obdiff[,num] < threshold)
# 
# ndiff <- ndiff[ndiff[,num] < threshold,]$feature
# ovdiff <- ovdiff[ovdiff[,num] < threshold,]$feature
# obdiff <- obdiff[obdiff[,num] < threshold,]$feature

####################################################################################################################
meta <- read.table("meta.txt")
mat <- read.table("met.txt")
mat <- mat[,meta$discover==1]
metname <- read.csv("processed_data/metabolite_info.csv")
# mat$name <- rownames(mat)
mat[1:5,1:5]

# s.mat <- t(apply(mat, 1, function(x){scale(x)}))
# colnames(s.mat) <- colnames(mat)
s.mat <- mat

index = 2
target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")
group.test <- c("Normal", "Overweight", "Obesity")

log2fc.threshold = 0.5
VIP.threshold = 1
pvalue.threshold = 0.05

make_diff <- function(s.mat, i, pvalue.threshold, log2fc.threshold, VIP.threshold) {
  threshold = log2fc.threshold
  VIP.threshold = VIP.threshold
  diff <- read.table(paste("Result/", save_path[index], "/", group.test[i], "_maaslin2.txt", sep = ""), 
                     header = T, row.names = 1)
  VIP <- read.table(paste("VIP_", group.test[i], ".txt", sep = ""))
  
  diff_output <- data.frame(
    feature <- rownames(s.mat),
    name = metname$MS2Metabolite[match(rownames(s.mat), metname$MS2kegg)],
    concentration = apply(s.mat, 1, median),
    pvalue = diff$pval[match(rownames(s.mat), diff$feature)],
    log2fc = apply(s.mat,1,function(x){log2(mean(na.omit(x[meta$Cancer == "CRC" & meta$Group == group.test[i]]))/(mean(na.omit(x[meta$Cancer == "Health" & meta$Group == group.test[i]]))))}),
    group = rep(group.test[i], nrow(s.mat))
  )  %>% arrange(pvalue) %>% mutate(fdr=p.adjust(pvalue,method = "BH")) %>%
    mutate(label=ifelse(log2fc>threshold,"UP",ifelse(log2fc<(-threshold),"DOWN","Nosig"))) %>%
    mutate(TPplotlabel = ifelse((fdr<0.05 & abs(log2fc)>threshold), name, NA)) %>%
    mutate(TPplotcolor = as.factor(paste(label,ifelse(is.na(TPplotlabel),0,1),sep = "")))
  
  diff_output$VIP <- VIP[,2][match(rownames(diff_output), VIP[,1])]
  
  diff_output <- diff_output %>%
    mutate(TPplotlabel = ifelse((pvalue < pvalue.threshold & abs(log2fc)>threshold & VIP > VIP.threshold), name, NA)) %>%
    mutate(TPplotcolor = as.factor(paste(label,ifelse(is.na(TPplotlabel),0,1),sep = "")))
  
  diff_output <- diff_output[!duplicated(diff_output$TPplotlabel) | is.na(diff_output$TPplotlabel),]
  return(diff_output)
}

ndiff <- make_diff(s.mat, 1, pvalue.threshold,log2fc.threshold, VIP.threshold)
ovdiff <- make_diff(s.mat, 2, pvalue.threshold, log2fc.threshold, VIP.threshold)
obdiff <- make_diff(s.mat, 3, pvalue.threshold, log2fc.threshold, VIP.threshold)

# temp <- rbind(obdiff[!is.na(obdiff$TPplotlabel),], ovdiff[!is.na(ovdiff$TPplotlabel),])
temp <- rbind(obdiff, ovdiff)
temp <- rbind(temp, ndiff)

temp$group[is.na(temp$TPplotlabel)] <- "Nosig0"
table(temp$group)

p <- ggplot(data=as.data.frame(temp), aes(x=log2fc, y=-log10(pvalue))) +
  geom_point(aes(color = group), alpha = 0.7) + # size = log10(concentration), 
  # scale_color_manual(labels = c("Nosig0"="No Significance","UP1"="Significantly Up","UP0"="Up","DOWN1"="Significantly Down","DOWN0"="Down"),
  #                    values = c("UP0"='#FAC0AE',"DOWN0"='#9BCFF0',"UP1"='#FA2311',"DOWN1"='#6175DB',"Nosig0" ='gray')) +
  scale_color_manual(values = c("Normal" = "#009FFD", "Overweight" = "#FFA400", "Obesity" = "#D00000"))+
  labs(x="log2 Fold Change",  y="-log10 pvalue") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.position = c(0.12, 0.75)) +
  theme(legend.title = element_blank(), 
        legend.key = element_rect(fill = 'transparent'),
        legend.background = element_rect(fill = 'transparent')) +
  geom_hline(yintercept = -log10(0.05),linetype=4, color = 'black', size = 0.5) +
  geom_vline(xintercept = c(-log2fc.threshold, log2fc.threshold), linetype = 4, color = 'black', size = 0.5) +
  geom_text_repel(data = temp, aes(x = log2fc, y = -log10(pvalue), label = TPplotlabel), size = 3, colour = "black",
                  max.overlaps = 30) +
  xlim(-x.axis.lim, x.axis.lim)##控制横坐标长度
  # theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
p

pdf(paste("vocano_met_total_mean_pvalue_", pvalue.threshold, "_VIP_", VIP.threshold, "_log2fc_", log2fc.threshold, ".pdf", sep = ""), width = 8, height = 8)
p
dev.off()

write.table(temp, paste("vocano_met_total_mean_pvalue_", pvalue.threshold, "_VIP_", VIP.threshold, "_log2fc_", log2fc.threshold, ".txt", sep = ""))

write.table(unique(temp$feature....rownames.s.mat.[temp$group != "Nosig0"]),
            "met_sig_volcano.txt")

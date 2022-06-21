library(corrplot)
library(dplyr)
require(readxl)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(data.table)
library(tibble)

library(igraph)
library(psych)
library(impute)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")

meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1,]

otu <- read.table("species.txt")
ev <- read.table("met.txt")
ko <- read.table("ko.txt")

otusig <- read.table("selected_species_274.txt") # species_heatmap_sig_0518.txt # final_selected_species.txt
evsig <- read.table("met_sig_volcano.txt")
kosig <- read.table("ko_sig_p0.05_full.txt")

table(rownames(otu) %in% otusig[,1])
rownames(otu) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(otu))

otu <- otu[otusig[,1],rownames(meta)]
ev <- ev[evsig[,1],rownames(meta)]
ko <- ko[kosig[,1],rownames(meta)]

taxinfo <- read.table("species_sig_taxonomy.txt")
taxinfo <- taxinfo[rownames(otu),]
taxorder <- order(taxinfo$domain)

# OTU vs OTU

cormat <- corr.test(t(otu[, meta$multigroup == group.test[1]]), method="spearman")
order.AOE = corrMatOrder(cormat$r, order = 'AOE')

lapply(1:6, function(num){
  # gotu <- otu[, meta$Group == group[num]]
  # gev <- ev[, meta$Group == group[num]]
  # gko <- ko[, meta$Group == group[num]]
  gotu <- otu[, meta$multigroup == group.test[num]]
  gev <- ev[, meta$multigroup == group.test[num]]
  gko <- ko[, meta$multigroup == group.test[num]]
  
  # mat <- as.matrix(rbind(gotu, gev))
  cormat <- corr.test(t(gotu), method="spearman")
  write.csv(cormat$r, file=paste("OTU_spearman_correlation_rmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  write.csv(cormat$p, file=paste("OTU_spearman_correlation_pmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  
  r <- -cormat$r[order.AOE, order.AOE] #taxorder
  p <- cormat$p[order.AOE, order.AOE]
  
  pdf(paste("OTU_spearman_correlation_", group.test[num],".pdf", sep=""), width = 40, height = 40)
  corrplot(r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = p, sig.level = 0.05, insig='label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu')
           #addCoef.col ='black', number.cex = 0.3
  ) # %>%
  #corrRect(name = c('.Eubacterium._rectale', 'C00074', 'K01190', 'K20257'))
  dev.off()
})

###############################
# otu <- read.table("species.txt")
# ev <- read.table("met.txt")
# ko <- read.table("ko.txt")
# 
# otu <- otu[,rownames(meta)]
# ev <- ev[,rownames(meta)]
# ko <- ko[,rownames(meta)]
# # otu <- otu["Anaerostipes_hadrus",]
# 
# # otu <- otu[, meta$Group == "Obesity"]
# # ev <- ev[, meta$Group == "Obesity"]
# # ko <- ko[, meta$Group == "Obesity"]
# 
# singlecorr <- corr.test(t(otu), use="pairwise",
#                         method="pearson", # 可选pearson/kendall
#                         adjust="fdr",
#                         alpha=0.05)
# 
# ob.r.single <- data.frame(spe.r = t(singlecorr$r), spe.p = t(singlecorr$p), spe.fdr = t(singlecorr$p.adj))
# colnames(ob.r.single) = c("ob.r", "ob.p", "ob.padj")
# comp <- round(cbind(r.single, ob.r.single),4)
# 
# comp$name <- metname$MS2Metabolite[match(rownames(r.single), metname$MS2kegg)]

##########################################################################33333
# Cross correlation/interaction
library(ggraph)
library(tidygraph)
library(igraph)

meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1,]

otu <- read.table("species.txt")
ev <- read.table("met.txt")
ko <- read.table("ko.txt")

otusig <- read.table("selected_species_274.txt") # species_heatmap_sig_0518.txt # final_selected_species.txt
evsig <- read.table("met_sig_volcano.txt")
kosig <- read.table("ko_sig_p0.05_full.txt")

table(rownames(otu) %in% otusig[,1])
rownames(otu) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(otu))

otu <- otu[otusig[,1],rownames(meta)]
ev <- ev[evsig[,1],rownames(meta)]
ko <- ko[kosig[,1],rownames(meta)]
# otu <- otu["Anaerostipes_hadrus",]
group <- c("Normal", "Overweight", "Obesity")
group.test <- c("COb", "HOb", "COv", "HOv", "CN", "HN")

rownames(ev) <- metname$MS2Metabolite[match(rownames(ev), metname$MS2kegg)]

taxinfo <- read.table("species_sig_taxonomy.txt")
taxinfo <- taxinfo[rownames(otu),]
taxorder <- order(taxinfo$domain)

otu <- otu[taxorder,]

lapply(1:6, function(num){
  # gotu <- otu[, meta$Group == group[num]]
  # gev <- ev[, meta$Group == group[num]]
  # gko <- ko[, meta$Group == group[num]]
  gotu <- otu[, meta$multigroup == group.test[num]]
  gev <- ev[, meta$multigroup == group.test[num]]
  gko <- ko[, meta$multigroup == group.test[num]]
  
  # mat <- as.matrix(rbind(gotu, gev))
  cormat <- corr.test(t(gotu), t(gev), method="spearman")
  write.csv(cormat$r, file=paste("OTUvsEV_spearman_correlation_rmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  write.csv(cormat$p, file=paste("OTUvsEV_spearman_correlation_pmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  
  r <- -cormat$r
  pdf(paste("OTUvsEV_spearman_correlation_", group.test[num],".pdf", sep=""), width = 30, height = 40)
  corrplot(-cormat$r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = cormat$p, sig.level = 0.05, insig='label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu')
           #addCoef.col ='black', number.cex = 0.3
           ) # %>%
    #corrRect(name = c('.Eubacterium._rectale', 'C00074', 'K01190', 'K20257'))
  dev.off()
  
  # mat <- as.matrix(rbind(gotu, gko))
  cormat <- corr.test(t(gotu), t(gko), method="spearman")
  write.csv(cormat$r, file=paste("OTUvsKO_spearman_correlation_rmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  write.csv(cormat$p, file=paste("OTUvsKO_spearman_correlation_pmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  
  r <- -cormat$r
  pdf(paste("OTUvsKO_spearman_correlation_", group.test[num],".pdf", sep=""), width = 40, height = 40)
  corrplot(r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = cormat$p, sig.level = 0.05, insig = 'label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu')
           #addCoef.col ='black', number.cex = 0.3
  ) # %>%
  #corrRect(name = c('.Eubacterium._rectale', 'C00074', 'K01190', 'K20257'))
  dev.off()
  
  cormat <- corr.test(t(gev), t(gko), method="spearman")
  write.csv(cormat$r, file=paste("EVvsKO_spearman_correlation_rmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  write.csv(cormat$p, file=paste("EVvsKO_spearman_correlation_pmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  
  r <- -cormat$r
  pdf(paste("EVvsKO_spearman_correlation_", group.test[num],".pdf", sep=""), width = 40, height = 30)
  corrplot(r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = cormat$p, sig.level = 0.05, insig = 'label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu')
           # addCoef.col ='black', number.cex = 0.3
  ) # %>%
  #corrRect(name = c('.Eubacterium._rectale', 'C00074', 'K01190', 'K20257'))
  dev.off()
  
})
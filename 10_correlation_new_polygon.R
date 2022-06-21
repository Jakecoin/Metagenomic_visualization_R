library(psych)
library(reshape2)
library(WGCNA)
library(ggrepel)

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}

target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")
group.test <- c("COb", "HOb", "COv", "HOv", "CN", "HN")

metname <- read.csv("processed_data/metabolite_info.csv")

meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1,]
otu <- read.table("species.txt")
ev <- read.table("met.txt")
ko <- read.table("ko.txt")

otusig <- read.table("selected_species_274.txt") # selected_species.txt 48 # selected_species_274.txt # spe_sig.txt # species_heatmap_sig_0518.txt # species_heatmap_sig_0524_simple.txt
evsig <- read.table("met_sig_volcano.txt") # met_sig_p0.05_0514.txt
kosig <- read.table("ko_sig_p0.05_0424.txt") #ko_sig_p0.05_full.txt

taxinfo <- read.table("species_sig_taxonomy.txt")

table(rownames(otu) %in% otusig[,1])
rownames(otu) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(otu))
table(rownames(otu) %in% otusig[,1])

otu <- otu[otusig[,1],rownames(meta)]
ev <- ev[evsig[,1],rownames(meta)]
ko <- ko[kosig[,1],rownames(meta)]

metname <- metname[rownames(ev) %in% metname$MS2kegg,]
t <- metname$MS2Metabolite[match(rownames(ev), metname$MS2kegg)]
# t[37] <- "Procaine"
# t[39] <- "Estradiol valerate"
# t[42] <- "Mesterolone"
rownames(ev) <- t

temp <- rbind(otu, ev)
temp <- rbind(temp, ko)



# start
rvalue <- 0.3
pvalue <- 0.005
point.del <- 1
withname <- FALSE
pdf(paste("network/", "ggB_F_48_r_", rvalue, "_p_", pvalue, "_",
          "_network_mean_0619.pdf", sep = ""),
    width = 10, height = 9)
# par(mfrow=c(3, 2)) #, mar=c(1,1,1,1)

lapply(1:6, function(i){
  idy <- meta$multigroup %in% group.test[i]
  # idy <- meta$Group %in% group.test[i]
  table(idy)
  mat <- t(temp[,idy])
  mat <- na.omit(mat)
  
  occor <- corr.test(mat,
                     use="pairwise",
                     method="spearman", # 可选pearson/kendall
                     adjust="fdr",
                     alpha=0.05)
  
  r_matrix <- occor$r
  p_matrix <- occor$p
  
  table(p_matrix < pvalue & abs(r_matrix) > rvalue)
  idx <- p_matrix < pvalue & abs(r_matrix) > rvalue
  r_matrix[!idx]=0
  idx1 <- 1:nrow(otu)
  idx2 <- (nrow(otu)+1):(nrow(otu)+nrow(ev)) #ev
  idx3 <- (nrow(otu)+nrow(ev)+1):nrow(temp) #ko
  r_matrix[idx1, idx1] = 0
  r_matrix[idx2, idx2] = 0
  r_matrix[idx3, idx3] = 0
  
  netClu = data.frame(ID = colnames(mat),
                      group = c(rep("Species", nrow(otu)), 
                                rep("Metabolite", nrow(ev)),
                                rep("KO", nrow(ko))))
  group1 <- taxinfo$domain[match(netClu$ID, taxinfo$species)]
  group1 <- ifelse(is.na(group1), netClu$group, group1)
  netClu$group = as.factor(group1)
  
  set.seed(12)
  
  result2 = PolygonClusterG(cor = r_matrix, nodeGroup = netClu, zoom = 0.8, zoom2 = 0.8) #PolygonRrClusterG
  node = result2[[1]]
  head(node)
  # ---node节点注释
  # nodes = nodeadd(plotcord = node, otu_table = otu_table, tax_table = tax_table)
  nodes <- node
  nodes$group <- netClu$group[match(node$elements, netClu$ID)]
  if(i <= 2){
    foldchange = apply(temp,1,function(x){log2(mean(na.omit(x[meta$multigroup == "COb"]))/(mean(na.omit(x[meta$multigroup == "HOb"]))))})
  } else if (i <= 4) {
    foldchange = apply(temp,1,function(x){log2(mean(na.omit(x[meta$multigroup == "COv"]))/(mean(na.omit(x[meta$multigroup == "HOv"]))))})
  } else if (i <= 6) {
    foldchange = apply(temp,1,function(x){log2(mean(na.omit(x[meta$multigroup == "CN"]))/(mean(na.omit(x[meta$multigroup == "HN"]))))})
  }
  
  nodes$updown <- ifelse(foldchange > 0, "Increased", "Decreased")
  
  # colnames(nodes)
  #-----计算边
  edge = edgeBuild(cor = r_matrix, node = node)
  head(edge)
  # edge$weight <- abseedge
  # colnames(edge)[8] = "cor"
  p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = cor),
                                data = edge, size = 0.2, alpha = 0.6) +
    geom_point(aes(X1, X2, fill = updown), pch = 21, data = nodes) + #, size = mean
    geom_text_repel(aes(X1, X2, label = elements), size = 2, nudge_y = -0.2, data = nodes,
                    max.overlaps = 20) +
    # geom_text(aes(X1, X2, label = elements), size = 2, nudge_y = -0.2, data = nodes) +
    scale_colour_manual(values = c("+" = "#D00000", "-" = "#009FFD")) +
    scale_fill_manual(values = c("Increased" = "#D00000", "Decreased" = "#009FFD")) +
    scale_size(range = c(4, 14)) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    ggtitle(group.test[i])
  p1
})

dev.off()


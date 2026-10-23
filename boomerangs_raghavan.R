pacman::p_load(data.table,dplyr,
               Seurat,
               ggplot2,
               ggrepel)

# csy_sigs <- fread("~/Projects/Johnsen_group/gene_signatures/basal_classical_sigs_CSY2020.tsv")
# 
# csy_indv <- fread("~/Projects/Johnsen_group/gene_signatures/csy_cla_bsl_sigs_1_2_6_10.tsv")
raghavan_sigs <- fread("~/Projects/Johnsen_group/gene_signatures/raghavan2021_scSubtype_scores.tsv")
#Overlapping_visualization
Chan_classical_all <- raghavan_sigs$scClassical#csy_sigs$classical
Chan_basal_all <- raghavan_sigs$scBasal#csy_indv$sig2bsl
IC_all <- raghavan_sigs$IC

Idents(seu) <- "scvi_samp_snn_res.0.2"

Ductal2 <- subset(seu, idents = c(1,10,17))#c(1,4,17,22,29,41,43))
Ductal2[["RNA"]] <- JoinLayers(Ductal2[["RNA"]])


# visualisation
DimPlot(Ductal2, raster=FALSE, reduction = "umap.scvi_samp", label = T) + 
  NoLegend()

## Module scoring with CSY signatures (same as tried above)
Ductal2 <- AddModuleScore(Ductal2, features = list(Chan_classical_all), 
                          name = "Chan_classical_all_score_", search = T)
Ductal2 <- AddModuleScore(Ductal2, features = list(Chan_basal_all), 
                          name = "Chan_basal_all_score_", search = T)
Ductal2 <- AddModuleScore(Ductal2, features = list(IC_all), 
                          name = "IC_all_score_", search = T)

FeaturePlot(Ductal2, features = c("Chan_classical_all_score_1", 
                                  "Chan_basal_all_score_1"), 
            order = T, label = T, reduction = "umap.scvi_samp") &
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n=11, 
                                                                name="RdBu")))

## Get expression matrix for signature genes
subtype_signature <- c(Chan_classical_all, Chan_basal_all, IC_all)
Ductal2_subtypegene <- Ductal2[subtype_signature,]
DefaultAssay(Ductal2_subtypegene) <- "RNA"
Ductal2_subtypegene


Ductal2_metadata <- Ductal2_subtypegene@meta.data
Ductal2_metadata <- Ductal2_metadata[ , c("sample_id", "Chan_classical_all_score_1", "Chan_basal_all_score_1", "Project")]
Ductal2_ExpMat <- GetAssayData(Ductal2_subtypegene, assay = "RNA", layer = "data")
Ductal2_ExpMat <- data.frame(Ductal2_ExpMat)
dim(Ductal2_ExpMat) 


## Unexpressed gene in scRNAseq
setdiff(subtype_signature, rownames(Ductal2_ExpMat))
# classical signature
intersect(Chan_classical_all, setdiff(subtype_signature, rownames(Ductal2_ExpMat)))
# basal-like signature
intersect(Chan_basal_all, setdiff(subtype_signature, rownames(Ductal2_ExpMat)))



### Correlation between subtype score and gene expression
Ductal2_ExpMat <- data.frame(t(Ductal2_ExpMat),
                               row.names = colnames(Ductal2_ExpMat))
colnames(Ductal2_ExpMat) <- gsub(pattern = "\\.", replacement = "-", colnames(Ductal2_ExpMat))
Ductal2_ExpMat <- cbind(Ductal2_metadata, Ductal2_ExpMat)
subtype_cor_df <- data.frame(matrix(ncol = 2, nrow = ncol(Ductal2_ExpMat)))  #depend on dim(Ductal2_ExpMat) # [1] 963 16715
colnames(subtype_cor_df) <- c("Chan_classical_all_score_1", "Chan_basal_all_score_1")


#need to change starting from gene names (IMPORTANT to change column number here instead of 5)
rownames(subtype_cor_df) <- colnames(Ductal2_ExpMat)#[1:ncol(Ductal2_ExpMat)]
for(i in colnames(Ductal2_ExpMat)[5:ncol(Ductal2_ExpMat)]){
  subtype_cor_df[i,"Chan_classical_all_score_1"] <- cor(Ductal2_ExpMat$Chan_classical_all_score_1, Ductal2_ExpMat[,i])
  subtype_cor_df[i,"Chan_basal_all_score_1"] <- cor(Ductal2_ExpMat$Chan_basal_all_score_1, Ductal2_ExpMat[,i])
}
subtype_cor_df <- na.omit(subtype_cor_df)



## Get average expression
Ductal2_subtypegene <- SetIdent(Ductal2_subtypegene, value = "Cell_type")
subtype_aveexp <- AverageExpression(Ductal2_subtypegene, assays = "RNA")
subtype_aveexp <- subtype_aveexp$RNA
subtype_aveexp <- data.frame(subtype_aveexp)
subtype_aveexp$Gene.name <- rownames(subtype_aveexp)
subtype_cor_df$Gene.name <- rownames(subtype_cor_df)
subtype_aveexp_cor <- inner_join(subtype_cor_df, subtype_aveexp, by="Gene.name")
rownames(subtype_aveexp_cor) <- subtype_aveexp_cor$Gene.name
subtype_aveexp_cor$Gene.name <- NULL



subtype_aveexp_cor_classical <- subtype_aveexp_cor[Chan_classical_all,]
subtype_aveexp_cor_classical <- na.omit(subtype_aveexp_cor_classical)
dim(subtype_aveexp_cor_classical) # [1] 478   3
subtype_aveexp_cor_basal <- subtype_aveexp_cor[Chan_basal_all,]
subtype_aveexp_cor_basal <- na.omit(subtype_aveexp_cor_basal)
dim(subtype_aveexp_cor_basal) # [1] 399  3

names(subtype_aveexp_cor_classical)


#low expression genes are not correlate subtype score
ggplot(subtype_aveexp_cor_classical, aes(x=Chan_basal_all_score_1, y=Chan_classical_all_score_1, col=log10(subtype_aveexp))) +
  geom_point(size=2) + 
  scale_color_viridis_b() + 
  geom_text_repel(aes(label = rownames(subtype_aveexp_cor_classical)), color = "black")

ggplot(subtype_aveexp_cor_basal, aes(x=Chan_basal_all_score_1, y=Chan_classical_all_score_1, col=log10(subtype_aveexp))) +
  geom_point(size=4) + 
  scale_color_viridis_b() + 
  geom_text_repel(aes(label = rownames(subtype_aveexp_cor_basal)), color = "black")


classical_score_summary <- summary(subtype_cor_df$Chan_classical_all_score_1)
classical_score_summary
basal_score_summary <- summary(subtype_cor_df$Chan_basal_all_score_1)
basal_score_summary

## top 25 genes correlate to subtype score
table(subtype_aveexp_cor$Chan_classical_all_score_1 > classical_score_summary[5] )
table(subtype_aveexp_cor$Chan_basal_all_score_1 > basal_score_summary[5] )

rownames(subtype_cor_df)[subtype_aveexp_cor$Chan_classical_all_score_1 > classical_score_summary[5]] -> classical_corHigh
rownames(subtype_cor_df)[subtype_aveexp_cor$Chan_basal_all_score_1 > basal_score_summary[5]] -> basal_corHigh

subtype_list <- list(Chan_classical_all, Chan_basal_all, classical_corHigh, basal_corHigh)
names(subtype_list) <- c("Chan_classical_all", "Chan_basal_all", "classical_new1", "basal_new1")
# library(VennDiagram)
# venn.plot <- venn.diagram(subtype_list, filename = "NULL", fill= c(2:5), disable.logging = T)

##  64 genes correlate both scores
table(subtype_aveexp_cor$Chan_classical_all_score_1 > classical_score_summary[4] & subtype_aveexp_cor$Chan_basal_all_score_1 > basal_score_summary[4])
# FALSE  TRUE 
# 823    51    

#both_corHigh <- rownames(subtype_cor_df)[subtype_aveexp_cor$Chan_classical_all_score_1 > classical_score_summary[4] & subtype_aveexp_cor$Chan_basal_all_score_1 > basal_score_summary[4]]
both_corHigh <- IC_all

## Get modified subtype signatures #it only take genes which is unique for classical or basal
classical_signature <- setdiff(classical_corHigh, both_corHigh)
basal_signature <- setdiff(basal_corHigh, both_corHigh)

raghavan_sigs <- fread("~/Projects/Johnsen_group/gene_signatures/raghavan2021_scSubtype_scores.tsv")

Ductal2 <- AddModuleScore(Ductal2, features = list(classical_signature), name = "classical_unique_score_", seed = 1)
Ductal2 <- AddModuleScore(Ductal2, features = list(basal_signature), name = "basal_unique_score_", seed = 1)
Ductal2 <- AddModuleScore(Ductal2, features = list(both_corHigh), name = "both_corHigh_", seed = 1)
Ductal2 <- AddModuleScore(Ductal2, features = raghavan_sigs, name = "raghavan", seed = 1, search = T)

FeaturePlot(Ductal2, features = c("classical_unique_score_1", "basal_unique_score_1", "both_corHigh_1"), 
            order = T, reduction = "umap.scvi_samp", label = T) &
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n=11, name="RdBu")))

FeaturePlot(Ductal2, features = c("scClassical", "scBasal", "IC"), 
            order = T, reduction = "umap.scvi_samp", label = T) &
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n=11, name="RdBu")))

#saveRDS(Ductal2, file = "Ductal2_all_score.RDS", ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

get_gene_state <- function(sobj, gene) {
  gene_info <- FetchData(Ductal2, vars = gene, layer = "counts")
  gene_pos_cells <- rownames(gene_info)[gene_info[[gene]] > 0]
  gene_neg_cells <- rownames(gene_info)[gene_info[[gene]] == 0]
  
  # gene_pos <- paste0(gene, "+")
  # gene_neg <- paste0(gene, "-")
  
  sobj[[gene]] <- ifelse(Cells(sobj) %in% gene_pos_cells, 1, 0)
  return(sobj)
}

Ductal2 <- get_gene_state(Ductal2, "SNAI2")
Ductal2 <- get_gene_state(Ductal2, "FAT2")
Ductal2 <- get_gene_state(Ductal2, "KRT5")
Ductal2 <- get_gene_state(Ductal2, "TP63")
Ductal2 <- get_gene_state(Ductal2, "HNF4A")
Ductal2 <- get_gene_state(Ductal2, "GATA6")
Ductal2 <- get_gene_state(Ductal2, "LGALS4")
Ductal2 <- get_gene_state(Ductal2, "REG4")



tumor.meta.data <- Ductal2@meta.data

#make Basal score negative
tumor.meta.data$flipped_Basal <- -1*tumor.meta.data$basal_unique_score_1
#Tumor_Type_score (X axis)
tumor.meta.data$Tumor_Type_score <- tumor.meta.data$flipped_Basal + tumor.meta.data$classical_unique_score_1
#Undifferentiated score (Y axis)
tumor.meta.data$Undiff_Score <- tumor.meta.data$both_corHigh_1 - abs(tumor.meta.data$Tumor_Type_score)



#wishbone plot

ggplot(tumor.meta.data) +
  geom_point(aes(Tumor_Type_score, Undiff_Score), size = 0.3)+
  theme(axis.text.x = element_text(size = 20, color = "black"), 
        axis.text.y = element_text(size = 20, color = "black"),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line = element_line(size = 0.8),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        aspect.ratio = 1,
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.5, "cm"),
        legend.text = element_text(size = 12),
        panel.spacing.x=unit(1, "lines"),
        panel.border = element_rect(colour="black",size=1, fill = NA))


gene_boomerang1 <- function(seurat_obj, tumor_metadata, feature, basal=T) {
  
  require(Seurat)
  require(ggplot2)
  
  #Labeling gene of interest in the wishbone plot 
  # feat_plot <- FeaturePlot(seurat_obj, features = feature)
  # feat_cells <- feat_plot$data
  # 
  # feat_score <- paste(feature, "score", sep = "_")
  metadata <- tumor_metadata
  metadata <- metadata |> arrange(.data[[feature]])
  #metadata[[feat_score]] <- feat_cells[[feature]]
  
  #x_jit <- diff(range(metadata$Tumor_Type_score, na.rm = TRUE)) * 0.05
  #y_jit <- diff(range(metadata$Undiff_Score, na.rm = TRUE)) * 0.05

  plot <- ggplot(metadata)+
    geom_point(aes(Tumor_Type_score, Undiff_Score, color = .data[[feature]],
               alpha = .data[[feature]]), size = 0.3,
               position = position_jitter(width = x_jit, height = y_jit)) +
    scale_alpha(range = c(0.01, 1), guide = "none")+
    theme(axis.text.x = element_text(size = 20, color = "black"), 
          axis.text.y = element_text(size = 20, color = "black"),
          strip.text.x = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "right",
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #axis.line = element_line(size = 0.8),
          axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13),
          aspect.ratio = 1,
          axis.ticks = element_line(color = "black"),
          axis.ticks.length = unit(0.5, "cm"),
          legend.text = element_text(size = 12),
          panel.spacing.x=unit(1, "lines"),
          panel.border = element_rect(colour="black",size=1, fill = NA),
          plot.title = element_text(size = 16, hjust = 0.5,
                                    margin = margin(b = 20)))
  
  high_col <- "blue2"#ifelse(basal, "orangered3", "blue2")
  
  plot <- plot +
    scale_color_gradient(
      low = "lightgrey",
      high = high_col,
      name = paste0(feature, " Gene Expression")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate("text", x = -1.5, y = Inf, label = "Basal",
             color = "orangered3", vjust = -0.5) +
    annotate("text", x =  1.5, y = Inf, label = "Classical",
             color = "blue2", vjust = -0.5) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(1,1,1,1, "cm"))
  
  print(plot + NoLegend() + ggtitle(feature))
  
}

gene_boomerang <- function(seurat_obj, tumor_metadata, feature, basal=T) {
  
  require(Seurat)
  require(ggplot2)
  
  #Labeling gene of interest in the wishbone plot 
  feat_plot <- FeaturePlot(seurat_obj, features = feature)
  feat_cells <- feat_plot$data
  
  feat_score <- paste(feature, "score", sep = "_")
  metadata <- tumor_metadata
  metadata[[feat_score]] <- feat_cells[[feature]] 
  metadata <- metadata |> arrange(.data[[feat_score]])
  
  #x_jit <- diff(range(metadata$Tumor_Type_score, na.rm = TRUE)) * 0.05
  #y_jit <- diff(range(metadata$Undiff_Score, na.rm = TRUE)) * 0.05
  
  plot <- ggplot(metadata)+
    geom_point(aes(Tumor_Type_score, Undiff_Score, color = .data[[feat_score]],
                   alpha = .data[[feat_score]]), size = 0.3,
               position = position_jitter(width = x_jit, height = y_jit)) +
    scale_alpha(range = c(0.01, 1), guide = "none")+
    theme(axis.text.x = element_text(size = 20, color = "black"), 
          axis.text.y = element_text(size = 20, color = "black"),
          strip.text.x = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "right",
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #axis.line = element_line(size = 0.8),
          axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13),
          aspect.ratio = 1,
          axis.ticks = element_line(color = "black"),
          axis.ticks.length = unit(0.5, "cm"),
          legend.text = element_text(size = 12),
          panel.spacing.x=unit(1, "lines"),
          panel.border = element_rect(colour="black",linewidth = 1, fill = NA),
          plot.title = element_text(size = 16, hjust = 0.5,
                                    margin = margin(b = 20)))
  
  high_col <- "blue2"#ifelse(basal, "orangered3", "blue2")
  
  plot <- plot +
    scale_color_gradient(
      low = "lightgrey",
      high = high_col,
      name = paste0(feature, " Expression")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate("text", x = -1.5, y = Inf, label = "Basal",
             color = "orangered3", vjust = -0.5) +
    annotate("text", x =  1.5, y = Inf, label = "Classical",
             color = "blue2", vjust = -0.5) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(1,1,1,1, "cm"))
  
  print(plot + ggtitle(feature))
  
}

#tumor.meta.data <- transform(tumor.meta.data, SNAI2_state = ifelse(snai2 == "SNAI2-", 0, 1))
#tumor.meta.data <- arrange(tumor.meta.data, SNAI2_score)

gene_boomerang1(Ductal2, tumor.meta.data, "SNAI2")
gene_boomerang1(Ductal2, tumor.meta.data, "KRT5")
gene_boomerang1(Ductal2, tumor.meta.data, "TP63")
gene_boomerang1(Ductal2, tumor.meta.data, "FAT2")

gene_boomerang1(Ductal2, tumor.meta.data, "HNF4A", basal=F)
gene_boomerang1(Ductal2, tumor.meta.data, "GATA6", basal = F)
gene_boomerang1(Ductal2, tumor.meta.data, "LGALS4", basal = F)
gene_boomerang1(Ductal2, tumor.meta.data, "REG4", basal = F)



gene_boomerang(Ductal2, tumor.meta.data, "rna_SNAI2")
gene_boomerang(Ductal2, tumor.meta.data, "rna_KRT5")
gene_boomerang(Ductal2, tumor.meta.data, "rna_TP63")
gene_boomerang(Ductal2, tumor.meta.data, "rna_FAT2")

gene_boomerang(Ductal2, tumor.meta.data, "rna_HNF4A", basal=F)
gene_boomerang(Ductal2, tumor.meta.data, "rna_GATA6", basal = F)
gene_boomerang(Ductal2, tumor.meta.data, "rna_LGALS4", basal = F)
gene_boomerang(Ductal2, tumor.meta.data, "rna_REG4", basal = F)

gene_boomerang(Ductal2, tumor.meta.data, "Chan_basal_all_score_1")
gene_boomerang(Ductal2, tumor.meta.data, "Chan_classical_all_score_1", basal = F)
gene_boomerang(Ductal2, tumor.meta.data, "basal_unique_score_1")
gene_boomerang(Ductal2, tumor.meta.data, "classical_unique_score_1", basal = F)


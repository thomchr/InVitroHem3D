source(file.path(Sys.getenv("EHT_REPO", "."), "config", "paths.R"))  # DATA_DIR, OUT_DIR, REPO_DIR
cell_topic <- read.table(file.path(DATA_DIR, "snATACseq/cisTopic/cell_topic.txt"),sep='\t')
colnames(cell_topic) <- cell_topic[1,]
rownames(cell_topic) <- cell_topic[,1]
cell_topic <- cell_topic[-1,-1]

cell_data=read.table(file.path(DATA_DIR, "snATACseq/HPC_HE_seurat.metadata.txt"))
cell_data <- cell_data[which(rownames(cell_data) %in% colnames(cell_topic)),]
cell_data$group="HE"
cell_data$group[cell_data$orig.ident=="HPC_seurat_project"]="HPC"

cells <- rownames(cell_data)
group.by <- "Refined_cell_type"
groups.use <- cell_data[cells, group.by, drop = FALSE]
group.use <- groups.use[, 1, drop = TRUE]
group.use <- factor(x = group.use)
names(x = group.use) <- cells
group.use <- group.use[c(which(group.use =="1_Stromal"), which(group.use =="2_Endo"), which(group.use =="3_HE"), which(group.use =="4_HPC"), which(group.use =="5_LMPP_Myeloid"), which(group.use =="6_LMPP_Lymphoid"),  which(group.use =="7_EMP"),  which(group.use =="8_EMP_Myeloid"),  which(group.use =="9_EMP_Mk"), which(group.use =="10_EMP_Ery"))]
order <- names(x = group.use)
cell_topic <- cell_topic[,order]

cell_topic <- t(apply(cell_topic, 1, as.numeric))

library(ComplexHeatmap)
colors <- scales::hue_pal()(length(unique(groups.use[,1])))
names(colors) <- sort(unique(groups.use[,1]))
mat <- cell_topic
colorPal <- grDevices::colorRampPalette(c('floralwhite', 'red', 'darkred'))
cellColor <- setNames(colors[groups.use[,1]], rownames(group.use))
cellColor <- cellColor[order]
colVars <- list()
colVars[['Cell_type']] <- colors
library(RColorBrewer)
c1 <- brewer.pal(4, 'Reds')[2]
c2 <- brewer.pal(3, 'Blues')[3]
colVars[['Sample']]  <- c(c1, c2)
names(colVars[['Sample']]) <- c("HPC","HE")
celltype <- cbind(as.data.frame(groups.use[order,1]),cell_data[order, 'group'])
colnames(celltype) <- c('Cell_type', 'Sample')
annotation <- ComplexHeatmap::HeatmapAnnotation(df = celltype, col = colVars, which='column')
heatmap <- ComplexHeatmap::Heatmap(data.matrix(mat), col=colorPal(20), cluster_columns = TRUE, cluster_rows = TRUE, 
				   show_row_dend = TRUE, show_column_names=FALSE, show_row_names = TRUE, 
				   top_annotation = annotation, name='Standarized AUC', row_names_gp = gpar(fontsize = 10), 
				   heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"), 
							       title_position='topcenter'), column_title_gp = gpar(fontface = 'bold'))

pdf(file.path(DATA_DIR, "snATACseq/topic_heatmap.pdf"))
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()






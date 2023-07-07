# Welcome to CompBio1!!
# Course instructor: Fahd Qadir Dragonmasterx87 on Github, for questions create an issue in the course repository
# First we will install packages required to run Seurat

# STEP1 Installation
# Install the remotes, dplyr, patchwork and dectools packages
install.packages('remotes')
install.packages("dplyr")
install.packages("patchwork")
install.packages("devtools")
install.packages("qs")
devtools::install_github("gaospecial/ggVennDiagram")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
BiocManager::install("dittoSeq")

# Install version 4.3.0 of Seurat
remotes::install_github("cran/spatstat.core")
remotes::install_version(package = 'Seurat', version = package_version('4.3.0'))

# Install a tutorial data library, if you cant access dont worry I have a backup
devtools::install_github("satijalab/seurat-data", ref = 'develop', force = TRUE)

# STEP2 Now load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(SeuratData)
  library(qs)
  library(ggplot2)
  library(harmony)
  set.seed(1234)
})

# Check
packageVersion("Seurat")
packageVersion("harmony")

# STEP3 can be bypassed and you can directly upload data from here ~Fahd_shared_with_participants\data proceed to STEP4
# STEP3 install dataset
InstallData("ifnb")

# load and process dataset
LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# Split two objects into 2
ctrl <- ifnb.list[["CTRL"]]
stim <- ifnb.list[["STIM"]]

# Lets randomly split data into n = 3 donors, remember this is just for simulation purposes
# the actual data is just n = 1
k <- 4
ctrl$k.assign <- sample(x = 1:k, size = ncol(ctrl), replace = TRUE)
ctrl.split <- SplitObject(ctrl, split.by = 'k.assign')
length(ctrl.split)
ctrl.split
names(ctrl.split) <- c("ctrl.d1", "ctrl.d2", "ctrl.d3", "ctrl.d4")

k <- 4
stim$k.assign <- sample(x = 1:k, size = ncol(stim), replace = TRUE)
stim.split <- SplitObject(stim, split.by = 'k.assign')
length(stim.split)
stim.split
names(stim.split) <- c("stim.d1", "stim.d2", "stim.d3", "stim.d4")

# Save files for tutorial
{
  qsave(ifnb.list[["ctrl.d1"]], r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\ctrl.d1.qs)")
  qsave(ifnb.list[["ctrl.d2"]], r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\ctrl.d2.qs)")
  qsave(ifnb.list[["ctrl.d3"]], r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\ctrl.d3.qs)")
  qsave(ifnb.list[["ctrl.d4"]], r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\ctrl.d4.qs)")

  qsave(ifnb.list[["stim.d1"]], r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\stim.d1.qs)")
  qsave(ifnb.list[["stim.d2"]], r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\stim.d2.qs)")
  qsave(ifnb.list[["stim.d3"]], r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\stim.d3.qs)")
  qsave(ifnb.list[["stim.d4"]], r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\stim.d4.qs)")
}

# STEP4 Load prepared Seurat files
{
  ctrl.d1 <- qread(r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\ctrl.d1.qs)")
  ctrl.d2 <- qread(r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\ctrl.d2.qs)")
  ctrl.d3 <- qread(r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\ctrl.d3.qs)")
  ctrl.d4 <- qread(r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\ctrl.d4.qs)")
  
  stim.d1 <- qread(r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\stim.d1.qs)")
  stim.d2 <- qread(r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\stim.d2.qs)")
  stim.d3 <- qread(r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\stim.d3.qs)")
  stim.d4 <- qread(r"(C:\Users\mqadir\Box\Courses-Workshops\CompBioW1\Fahd_shared_with_participants\data\stim.d4.qs)")
}

# STEP5 Create a unified list, remember object name comes first in list notation
ifnb.list <- list("ctrl.d1" = ctrl.d1, "ctrl.d2" = ctrl.d2, "ctrl.d3" = ctrl.d3, "ctrl.d4" = ctrl.d4,
                  "stim.d1" = stim.d1, "stim.d2" = stim.d2, "stim.d3" = stim.d3, "stim.d4" = stim.d4)

# STEP6 Addition of donor metadata
ifnb.list[["ctrl.d1"]]$donor <- "d1"
ifnb.list[["ctrl.d2"]]$donor <- "d2"
ifnb.list[["ctrl.d3"]]$donor <- "d3"
ifnb.list[["ctrl.d4"]]$donor <- "d4"

ifnb.list[["stim.d1"]]$donor <- "d1"
ifnb.list[["stim.d2"]]$donor <- "d2"
ifnb.list[["stim.d3"]]$donor <- "d3"
ifnb.list[["stim.d4"]]$donor <- "d4"

# STEP7 Merge objects
pbmc <- merge(ifnb.list[["ctrl.d1"]], y = c(ifnb.list[["ctrl.d2"]], ifnb.list[["ctrl.d3"]], ifnb.list[["ctrl.d4"]],
                                            ifnb.list[["stim.d1"]], ifnb.list[["stim.d2"]], ifnb.list[["stim.d3"]], ifnb.list[["stim.d4"]]), 
              add.cell.ids = c("ctrl.d1", "ctrl.d2", "ctrl.d3", "ctrl.d4",
                               "stim.d1", "stim.d2", "stim.d3", "stim.d4"), project = "pbmc")

# STEP8 QC
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
grep ("^CCL", rownames(pbmc[["RNA"]]),value = T)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-") # this dataset doesnt contain MT DNA, otherwise we subset on <10% MT

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Subset data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & nCount_RNA < 6000)

# Lets visualize new QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# STEP9 Data processing
# normalize and identify variable features for each dataset independently
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# select highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc, assay = "RNA")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose = TRUE)

# STEP10 Linear dimensionality reduction
# Linear dimensional reduction
pbmc <- RunPCA(pbmc, pc.genes = pbmc@assays$RNA@var.features, npcs = 20, verbose = TRUE)

# Visualize PCA
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)

# Visualize elbow plot of PC
ElbowPlot(pbmc)

# STEP11 Batch correction using Harmony
#Run Harmony batch correction with library and tissue source covariates
pbmc <- RunHarmony(pbmc,
                   assay.use = "RNA",
                   reduction = "pca",
                   dims.use = 1:20,
                   group.by.vars = c("donor", "stim"),
                   kmeans_init_nstart=20, kmeans_init_iter_max=100,
                   plot_convergence = TRUE)

# STEP12 Non linear multidimensional projection using UMAP
# Run UMAP, on PCA NON-batch corrected data
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20, return.model = TRUE)
DimPlot(pbmc, reduction = 'umap', label = FALSE, pt.size = 2, raster=TRUE)

# Now run Harmony
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20, return.model = TRUE)
DimPlot(pbmc, reduction = 'umap', label = FALSE, pt.size = 2, raster=TRUE)

# STEP13 Clustering
# algorithm 3 is the smart local moving (SLM) algorithm https://link.springer.com/article/10.1140/epjb/e2013-40829-0
pbmc <- pbmc %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=3,resolution = c(0.5), method = 'igraph') #25 res

Idents(pbmc) <- "seurat_annotations"
DimPlot(pbmc, reduction = "umap", label = TRUE)

# Observe gene exprssion
FeaturePlot(pbmc, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                               "CCL2", "PPBP"), min.cutoff = "q9")

# Rename clusters
Idents(pbmc) <- "RNA_snn_res.0.5"
table(pbmc@meta.data[["RNA_snn_res.0.5"]])
DimPlot(pbmc, reduction = "umap", label = TRUE)

new.cluster.ids <- c("CD14 Mono", "CD4 Naive T", "CD4 Memory T", "CD16 Mono", 
                     "B", "CD8 T", "T activated", "NK", "DC", "B Activated",
                     "Mk", "pDC", "Mono/Mk Doublets", "Eryth")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype <- Idents(pbmc)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# STEP13 Differential gene testing
# STEP13A: Single cell gene testing
# In order to fine DE genes its important to annotate cells,
# this is where metadata comes to be important
# Our first method employs single cell differential gene expression using LR
pbmc$celltype.stim <- paste(pbmc$celltype, pbmc$stim, sep = "_")
Idents(pbmc) <- "celltype.stim"
b.interferon.response <- FindMarkers(pbmc, ident.1 = "B_STIM", ident.2 = "B_CTRL",
                                     slot = 'data',
                                     test.use = "LR",
                                     min.pct = 0.1,
                                     latent.vars = c("donor"),
                                     logfc.threshold = 0.5849, #~1.5FC
                                     only.pos = TRUE,
                                     verbose = FALSE)

head(b.interferon.response, n = 15)
b.interferon.response <- dplyr::filter(b.interferon.response, p_val_adj < 5e-2)
plots <- VlnPlot(pbmc, features = c("ISG15", "ISG20"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# # STEP13B: Pseudobulk testing
DefaultAssay(pbmc) <- "RNA"
pbmc$celltype.stim.donor <- paste(pbmc$celltype, pbmc$stim, pbmc$donor, sep = "_")
Idents(pbmc) <- "celltype.stim.donor"
table(pbmc[["celltype.stim.donor"]])
combined_pbmc <- AggregateExpression(pbmc, 
                                     assays = c("RNA"), 
                                     features = NULL, return.seurat = TRUE,  
                                     group.by = "celltype.stim.donor",
                                     slot = "counts", verbose = FALSE)

combined_pbmc$celltype.stim.donor <- Cells(combined_pbmc)

# Metadata organization and addition to aggregated object
{
  Idents(combined_pbmc) <- 'celltype.stim.donor'
  combined_pbmc$celltype <- combined_pbmc@meta.data[["orig.ident"]]
  metadat <- combined_pbmc@meta.data
  metadat$celltype <- metadat[c('celltype')] <- str_split_i(metadat$celltype.stim.donor, "_", -3)
  metadat$stim <- metadat[c('stim')] <- str_split_i(metadat$celltype.stim.donor, '_', -2)
  metadat$donor <- metadat[c('donor')] <- str_split_i(metadat$celltype.stim.donor, '_', -1)
  combined_pbmc@meta.data = metadat
}

table(combined_pbmc@meta.data[["celltype"]])
table(combined_pbmc@meta.data[["donor"]])
table(combined_pbmc@meta.data[["celltype.stim.donor"]])

combined_pbmc$celltype.stim <- paste(combined_pbmc$celltype, combined_pbmc$stim, sep = "_")
table(combined_pbmc@meta.data[["celltype.stim"]])

Idents(combined_pbmc) <- "celltype.stim"
b.interferon.response.aggr <- FindMarkers(combined_pbmc, ident.1 = "B_STIM", ident.2 = "B_CTRL",
                                          slot = 'data',
                                          test.use = "DESeq2",
                                          min.pct = 0.1,
                                          latent.vars = c("donor"),
                                          logfc.threshold = 0.5849, #~1.5FC
                                          only.pos = TRUE,
                                          verbose = FALSE)

head(b.interferon.response.aggr, n = 15)

# Why is my adj pval 0?
.Machine$double.xmin

# Replace o with low pval, subset for 0.05 FDR and FC > 1.5
b.interferon.response.aggr$p_val_adj[b.interferon.response.aggr$p_val_adj == 0] <- 2e-302
b.interferon.response.aggr <- dplyr::filter(b.interferon.response.aggr, p_val_adj < 5e-2)
b.interferon.response.aggr <- dplyr::filter(b.interferon.response.aggr, avg_log2FC > 0.5849)

intersect(rownames(b.interferon.response), rownames(b.interferon.response.aggr))

# Lets look at genes in a venn diagram'
LR.genes <- rownames(b.interferon.response)
DESeq2.genes <- rownames(b.interferon.response.aggr)

# Make a list
x <- list("LR.genes" = LR.genes, "DESeq2.genes" = DESeq2.genes)

# Make a Venn object
venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, "(", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 3) +
  scale_fill_gradient(low = "white", high = "dodgerblue3") + # change color based on celltype
  # scale_color_manual(values = c("bmvsbf" = "black",
  #                               "mvsf" ="black", 
  #                               "wmvswf" = 'black'),
  # scale_color_manual(values = c("beta_m" = "black",
  #                               "beta_f" ="black",
  #                               "alpha_m" = 'black',
  #                               "alpha_f" = 'black')) +
  scale_color_manual(values = c("LR" = "black",
                                "DESeq2" ="black"),
                     labels = c('D' = 'D = bdiv_human')) +
  theme_void()

# Look at all sets of genes forming overlaps
# https://github.com/yanlinlin82/ggvenn/issues/21
mylist <- data@region[["item"]]
names(mylist)
names(mylist) <- data@region[["name"]]
mylist


# Visualization
Idents(pbmc) <- factor(Idents(pbmc), levels = c("Mono/Mk Doublets", "pDC", "Eryth", 
                                                "Mk", "DC", "CD14 Mono", "CD16 Mono", 
                                                "B Activated", "B", "CD8 T", "NK", "T activated",
                                                "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
Idents(pbmc) <- "celltype"
DotPlot(pbmc, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()

# Change colors on UMAP
DimPlot(pbmc, #switch here to plot
        #split.by = "Diabetes Status", 
        group.by = "celltype", 
        label = FALSE, 
        ncol = 1, 
        raster = FALSE,
        pt.size = 0.05,
        cols = c("dodgerblue3",      
                 "turquoise2",       
                 "lightseagreen",    
                 "darkseagreen2",    
                 "khaki2",            
                 "springgreen4",     
                 "chartreuse3",      
                 "burlywood3",       
                 "darkorange2",      
                 "salmon3",          
                 "orange",           
                 "salmon",           
                 "red",              
                 "magenta3",         
                 "orchid1",          
                 "red4",             
                 "grey30"            
        )
)

# We can also look at cellular number across various cell types
pbmc$donor.stim <- paste(pbmc$donor, pbmc$stim, sep = "_")
table(pbmc$donor.stim)
dittoBarPlot(pbmc, "celltype", 
             retain.factor.levels = TRUE,
             scale = "count",
             color.panel = c("dodgerblue3",      
                             "turquoise2",       
                             "lightseagreen",   
                             "darkseagreen2",   
                             "khaki2",           
                             "springgreen4",     
                             "chartreuse3",      
                             "burlywood3",      
                             "darkorange2",      
                             "salmon3",         
                             "orange",           
                             "salmon",           
                             "red",              
                             "magenta3",         
                             "orchid1",          
                             "red4",             
                             "grey30"),          
             group.by = "donor.stim") + coord_flip()

# Or even cellular proportion
dittoBarPlot(pbmc, "celltype", 
             retain.factor.levels = TRUE,
             scale = "percent",
             color.panel = c("dodgerblue3",      
                             "turquoise2",       
                             "lightseagreen",   
                             "darkseagreen2",   
                             "khaki2",           
                             "springgreen4",     
                             "chartreuse3",      
                             "burlywood3",      
                             "darkorange2",      
                             "salmon3",         
                             "orange",           
                             "salmon",           
                             "red",              
                             "magenta3",         
                             "orchid1",          
                             "red4",             
                             "grey30"),          
             group.by = "donor.stim") + coord_flip()


# Heatmap
label_genes <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                 "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                 "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
genes.to.plot <- pbmc@assays[["RNA"]]@var.features
dittoHeatmap(
  combined_pbmc,
  genes = genes.to.plot,
  # metas = NULL,
  # cells.use = NULL,
  annot.by = c("celltype", "donor", "stim"),
  #annot.by = c("lib", "sex", "source"),
  order.by = c("celltype"),
  # main = NA,
  # cell.names.meta = NULL,
  # assay = .default_assay(object),
  # slot = .default_slot(object),
  # swap.rownames = NULL,
  heatmap.colors = colorRampPalette(c("dodgerblue", "white", "red3"))(50),
  breaks=seq(-2, 2, length.out=50),
  scaled.to.max = FALSE,
  # heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
  # annot.colors = c(dittoColors(), dittoColors(1)[seq_len(7)]),
  # annotation_col = NULL,
  annotation_colors = list(celltype = c("CD14 Mono" = "salmon3",
                                        "CD4 Naive T" = "orange",
                                        "CD4 Memory T"= "lightseagreen",
                                        "CD16 Mono" = "dodgerblue3",
                                        "B" = "turquoise2",
                                        "CD8 T" = "burlywood3",
                                        "T activated" = "darkseagreen2",
                                        "NK" = "chartreuse3",
                                        "DC" = "darkorange2",
                                        "B Activated" = "red",
                                        "Mk" = "khaki2",
                                        "pDC" = "springgreen4",
                                        "Mono/Mk Doublets" = "orchid1",
                                        "Eryth" = "magenta3"),
                           donor = c("d1" = "dodgerblue",
                                     "d2" = "red2",
                                     "d3" = "green3",
                                     "d4" = "purple2"),
                           stim = c("CTRL" = "red4",
                                    "STIM" = "deepskyblue3")),
                           # ancestry = c("white" = "deepskyblue3",
                           #              "black" = "black",
                           #              "hispanic" = "darkorange"),
                           # source = c("nPOD" = "dodgerblue",
                           #            "Tulane" = "springgreen4",         
                           #            "UPENN" = "red4")),
  # # data.out = FALSE,
  # highlight.features = NULL,
  # highlight.genes = NULL,
  # show_colnames = isBulk(object),
  # show_rownames = TRUE,
  # scale = "row",
  #cluster_cols = TRUE,
  # border_color = NA,
  # legend_breaks = NA,
  # drop_levels = FALSE,
  # breaks = NA,
  # complex = FALSE
  #gaps_col = c(460),
  complex = TRUE,
  use_raster = TRUE,
  raster_quality = 5
) + rowAnnotation(mark = anno_mark(at = match(label_genes, 
                                              rownames(pbmc[genes.to.plot,])), 
                                   labels = label_genes, 
                                   which = "row",
                                   labels_gp = list(cex=1),
                                   #link_width = unit(4, "mm"), link_height = unit(4, "mm"),
                                   padding = 0.1))







---
title: Normalization using VST and SC Transform
layout: default
nav_order: 2
parent: Day 2
---
## Normalization using VST
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

#
FeaturePlot(pbmc, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                               "CCL2", "PPBP"), min.cutoff = "q9")
```
----

[Just the Docs]: https://just-the-docs.github.io/just-the-docs/
[GitHub Pages]: https://docs.github.com/en/pages
[README]: https://github.com/just-the-docs/just-the-docs-template/blob/main/README.md
[Jekyll]: https://jekyllrb.com
[GitHub Pages / Actions workflow]: https://github.blog/changelog/2022-07-27-github-pages-custom-github-actions-workflows-beta/
[use this template]: https://github.com/just-the-docs/just-the-docs-template/generate

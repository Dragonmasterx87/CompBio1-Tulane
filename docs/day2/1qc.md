---
title: Data loading/QC/filtering
layout: default
nav_order: 1
parent: Day 2
---

## Data loading/QC/filtering
## Here we will load a dataset availabale for download as part of Seurat's `SeuratData` tutorial data collection
```r
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
  set.seed(1234)
})

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

# STEP5 Addition of donor metadata
ifnb.list <- c(ctrl.split, stim.split)
ifnb.list[["ctrl.d1"]]$donor <- "d1"
ifnb.list[["ctrl.d2"]]$donor <- "d2"
ifnb.list[["ctrl.d3"]]$donor <- "d3"
ifnb.list[["ctrl.d4"]]$donor <- "d4"

ifnb.list[["stim.d1"]]$donor <- "d1"
ifnb.list[["stim.d2"]]$donor <- "d2"
ifnb.list[["stim.d3"]]$donor <- "d3"
ifnb.list[["stim.d4"]]$donor <- "d4"

# STEP6 Create a unified list, remember object name comes first in list notation
ifnb.list <- list("ctrl.d1" = ctrl.d1, "ctrl.d2" = ctrl.d2, "ctrl.d3" = ctrl.d3, "ctrl.d4" = ctrl.d4,
                  "stim.d1" = stim.d1, "stim.d2" = stim.d2, "stim.d3" = stim.d3, "stim.d4" = stim.d4)

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

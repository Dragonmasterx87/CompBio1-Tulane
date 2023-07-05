---
title: Differential testing framework
layout: default
nav_order: 4
parent: Day 2
---
## Differential testing framework
### STEP13 Differential gene testing
### STEP13A: Single cell gene testing
```r
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
```

### STEP13B: Pseudobulk testing
```r
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
.Machine

# Replace o with low pval, subset for 0.05 FDR and FC > 1.5
b.interferon.response.aggr$p_val_adj[b.interferon.response.aggr$p_val_adj == 0] <- 2e-302
b.interferon.response.aggr <- dplyr::filter(b.interferon.response.aggr, p_val_adj < 5e-2)
b.interferon.response.aggr <- dplyr::filter(b.interferon.response.aggr, avg_log2FC > 0.5849)

intersect(rownames(b.interferon.response), rownames(b.interferon.response.aggr))
```

----

[Just the Docs]: https://just-the-docs.github.io/just-the-docs/
[GitHub Pages]: https://docs.github.com/en/pages
[README]: https://github.com/just-the-docs/just-the-docs-template/blob/main/README.md
[Jekyll]: https://jekyllrb.com
[GitHub Pages / Actions workflow]: https://github.blog/changelog/2022-07-27-github-pages-custom-github-actions-workflows-beta/
[use this template]: https://github.com/just-the-docs/just-the-docs-template/generate

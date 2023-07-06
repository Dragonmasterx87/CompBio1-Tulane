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
Warning message:
glm.fit: fitted probabilities numerically 0 or 1 occurred 

head(b.interferon.response, n = 15)
                p_val avg_log2FC pct.1 pct.2     p_val_adj
ISG15   1.886870e-251  4.5777825 0.998 0.244 2.651618e-247
ISG20   6.775871e-233  2.9286888 1.000 0.674 9.522132e-229
IFIT3   2.249518e-226  4.5012181 0.963 0.052 3.161247e-222
IFI6    7.675338e-222  4.2569268 0.965 0.077 1.078615e-217
IFIT1   1.235784e-197  4.1209722 0.908 0.032 1.736647e-193
MX1     1.447618e-160  3.2484413 0.908 0.116 2.034337e-156
TNFSF10 1.295544e-148  3.7724262 0.781 0.022 1.820627e-144
LY6E    1.309127e-148  3.1014142 0.894 0.151 1.839716e-144
IFIT2   4.554090e-142  3.6335959 0.786 0.037 6.399863e-138
B2M     2.530957e-119  0.6158652 1.000 1.000 3.556754e-115
CXCL10  4.799024e-113  5.3217240 0.644 0.010 6.744069e-109
PLSCR1  2.659608e-112  2.8015323 0.788 0.119 3.737548e-108
IRF7    3.667981e-108  2.5623663 0.837 0.193 5.154613e-104
HERC5    6.244788e-98  2.8136999 0.611 0.022  8.775800e-94
UBE2L6   7.032047e-91  2.1202128 0.857 0.301  9.882136e-87

b.interferon.response <- dplyr::filter(b.interferon.response, p_val_adj < 5e-2)
plots <- VlnPlot(pbmc, features = c("ISG15", "ISG20"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
The default behaviour of split.by has changed.
Separate violin plots are now plotted side-by-side.
To restore the old behaviour of a single split violin,
set split.plot = TRUE.
      
This message will be shown once per session.
wrap_plots(plots = plots, ncol = 1)
```
![](../../assets/images/vlnplt3.JPG)

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

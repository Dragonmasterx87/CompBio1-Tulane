---
title: Differential testing framework
layout: default
nav_order: 4
parent: Day 2
---
## Differential testing framework
### STEP13 Differential gene testing
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

----

[Just the Docs]: https://just-the-docs.github.io/just-the-docs/
[GitHub Pages]: https://docs.github.com/en/pages
[README]: https://github.com/just-the-docs/just-the-docs-template/blob/main/README.md
[Jekyll]: https://jekyllrb.com
[GitHub Pages / Actions workflow]: https://github.blog/changelog/2022-07-27-github-pages-custom-github-actions-workflows-beta/
[use this template]: https://github.com/just-the-docs/just-the-docs-template/generate

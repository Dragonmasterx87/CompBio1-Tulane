---
title: Visualization
layout: default
nav_order: 5
parent: Day 2
---
## Visualization
```r
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


# Note a how 70% of genes overlap, but who are they?


# Look at all sets of genes forming overlaps
# https://github.com/yanlinlin82/ggvenn/issues/21
mylist <- data@region[["item"]]
names(mylist)
names(mylist) <- data@region[["name"]]
mylist
$LR.genes
 [1] "RSAD2"    "OASL"     "PML"      "ETV7"     "DHX58"    "STAT2"    "DDX60"    "DEK"     
 [9] "ZBP1"     "ZNFX1"    "GSDMD"    "MTHFD2"   "HSP90AA1" "ELF1"     "B3GNT2"   "NDUFA9"  
[17] "LACTB"    "AP3B1"    "DTX3L"    "TMEM140"  "GBP7"     "GLRX"     "SP140"    "TMX1"    
[25] "BAG1"     "TAP2.1"   "ZNF267"   "CNDP2"    "STAP1"    "CD53"     "HIF1A"    "MYCBP2"  
[33] "CHI3L2"   "GNB4"     "SYNE2"    "RARRES3"  "GPBP1"    "SNX6"     "CTSS"     "TANK"    
[41] "IDH3A"    "SOD2"     "MCL1"     "CAST"     "HLA-F"    "EVL"      "RNASEH2B" "UBE2F"   
[49] "ANXA2R"   "IRF2"     "SQRDL"    "TSPAN13" 

$DESeq2.genes
[1] "HLA-B" "H3F3B" "HLA-C" "LAMP3"

$LR.genes..DESeq2.genes
  [1] "ISG15"     "ISG20"     "IFIT3"     "IFI6"      "IFIT1"     "MX1"       "TNFSF10"   "LY6E"     
  [9] "IFIT2"     "B2M"       "CXCL10"    "PLSCR1"    "IRF7"      "HERC5"     "UBE2L6"    "IFI44L"   
 [17] "EPSTI1"    "OAS1"      "GBP1"      "IFITM2"    "SAMD9L"    "NT5C3A"    "IFI35"     "PSMB9"    
 [25] "MX2"       "DYNLT1"    "BST2"      "IFITM3"    "CMPK2"     "SAT1"      "EIF2AK2"   "PPM1K"    
 [33] "GBP4"      "DDX58"     "PSMA2.1"   "LAP3"      "SAMD9"     "XAF1"      "IFI16"     "COX5A"    
 [41] "SOCS1"     "MYL12A"    "SP110"     "PARP14"    "PSME2"     "TMSB10"    "CHST12"    "FBXO6"    
 [49] "MT2A"      "PLAC8"     "TRIM22"    "DRAP1"     "SUB1"      "TNFSF13B"  "NMI"       "XRN1"     
 [57] "NEXN"      "RBCK1"     "CLEC2D"    "MNDA"      "RNF213"    "IFI44"     "GBP5"      "NPC2"     
 [65] "STAT1"     "WARS"      "OAS2"      "SELL"      "TAP1"      "DDX60L"    "IRF8"      "OAS3"     
 [73] "RTCB"      "IFITM1"    "KIAA0040"  "CXCL11"    "CARD16"    "PSMA4"     "DNAJA1"    "IFIH1"    
 [81] "TYMP"      "HLA-E"     "LGALS9"    "NUB1"      "C19orf66"  "GBP2"      "PSMB8"     "GNG5"     
 [89] "HAPLN3"    "PMAIP1"    "IFIT5"     "PARP9"     "CD38"      "GMPR"      "C5orf56"   "EAF2"     
 [97] "HERC6"     "CD48"      "RTP4"      "RABGAP1L"  "USP30-AS1" "TREX1"     "IGFBP4"    "INPP1"    
[105] "CCL8"      "CREM"      "CD164"     "CLIC1"     "APOL6"     "CCL2"      "SMCHD1"    "ODF2L"    
[113] "EHD4"      "NAPA"      "SP100"     "PHF11"     "FAM46A"    "PNPT1"     "ADAR"      "POMP"     
[121] "UNC93B1"   "DCK"       "IRF1"      "CCR7"      "TMEM123"   "CASP4"     "HSH2D"     "SLFN5"    
[129] "CFLAR"     "CASP1"     "VAMP5"     "PSME1"     "XBP1"      "MRPL44" 



# Visualization
Idents(pbmc) <- factor(Idents(pbmc), levels = c("Mono/Mk Doublets", "pDC", "Eryth", 
                                                "Mk", "DC", "CD14 Mono", "CD16 Mono", 
                                                "B Activated", "B", "CD8 T", "NK", "T activated",
                                                "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
DotPlot(pbmc, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()
```
----

[Just the Docs]: https://just-the-docs.github.io/just-the-docs/
[GitHub Pages]: https://docs.github.com/en/pages
[README]: https://github.com/just-the-docs/just-the-docs-template/blob/main/README.md
[Jekyll]: https://jekyllrb.com
[GitHub Pages / Actions workflow]: https://github.blog/changelog/2022-07-27-github-pages-custom-github-actions-workflows-beta/
[use this template]: https://github.com/just-the-docs/just-the-docs-template/generate

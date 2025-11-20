## 07_Revision.R
## Mapping pregnancy neutrophil scRNA-seq to Ostuni reference
## and NP vs GDM comparisons within LD and HD fractions.

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(stringr)
  library(openxlsx2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(fgsea)
  library(EnhancedVolcano)
})

rm(list = ls())
graphics.off()
gc()
set.seed(123)

## ---------------------------------------------------------------------------
## Section 0: Paths and basic setup
## ---------------------------------------------------------------------------

# Project directories (relative to repo root)
dir_data    <- "data"
dir_figures <- "figures"
dir_results <- "results"

# Ensure output directories exist
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)

# Input objects (adjust relative paths if needed)
path_fig6 <- file.path(dir_data, "reference", "Fig6_adjusted.rds")
path_seu  <- file.path(dir_data, "seurat",   "04_GDMvsNP_allCellType.rds")

# Output RDS paths
path_seu_mapped     <- file.path(dir_data, "seurat",   "07_seu_mapped.rds")
path_fig6_processed <- file.path(dir_data, "reference", "Fig6_processed.rds")

## ---------------------------------------------------------------------------
## Section 1: Map scRNA-seq data to reference
## ---------------------------------------------------------------------------

# 1) Load reference and query
fig6 <- readRDS(path_fig6)
seu  <- readRDS(path_seu)

# Use published cluster annotation as identity on the reference
Idents(fig6) <- fig6$clustering_paper

# Work on SCT assay for reference; re-run SCT for query below
DefaultAssay(fig6) <- "SCT"

# Re-run SCT on query (assumes RNA assay is present)
DefaultAssay(seu) <- "RNA"
if (!"percent.mt" %in% colnames(seu@meta.data)) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
}
cc_genes <- Seurat::cc.genes.updated.2019
seu <- NormalizeData(seu, verbose = FALSE)
seu <- CellCycleScoring(
  seu,
  s.features   = intersect(cc_genes$s.genes, rownames(seu)),
  g2m.features = intersect(cc_genes$g2m.genes, rownames(seu)),
  set.ident    = FALSE
)
vars_to_regress <- intersect(
  c("percent.mt", "percent_mito", "S.Score", "G2M.Score"),
  colnames(seu@meta.data)
)
seu <- SCTransform(
  seu,
  vst.flavor = "v2",
  return.only.var.genes = FALSE,
  vars.to.regress = vars_to_regress
)
DefaultAssay(seu) <- "SCT"

# 2) SCT feature selection + PrepSCTIntegration (do not overwrite originals)
features  <- SelectIntegrationFeatures(object.list = list(fig6, seu), nfeatures = 3000)
prep_list <- PrepSCTIntegration(object.list = list(fig6, seu), anchor.features = features)
fig6_prep <- prep_list[[1]]
seu_prep  <- prep_list[[2]]
DefaultAssay(fig6_prep) <- "SCT"
DefaultAssay(seu_prep)  <- "SCT"

# Carry labels and identities
fig6_prep$clustering_paper <- fig6$clustering_paper
Idents(fig6_prep) <- fig6_prep$clustering_paper

# Reuse reference PCA; fit new UMAP with projection model
dims_use <- 1:20
fig6_prep@reductions$pca <- fig6@reductions$pca
fig6_prep <- RunUMAP(
  fig6_prep,
  reduction    = "pca",
  dims         = dims_use,
  return.model = TRUE,
  seed.use     = 123
)

# 3) Transfer anchors (SCT) + MapQuery to reference UMAP
anchors <- FindTransferAnchors(
  reference            = fig6_prep,
  query                = seu_prep,
  normalization.method = "SCT",
  reference.reduction  = "pca",
  dims                 = dims_use
)

seu_mapped <- MapQuery(
  anchorset           = anchors,
  reference           = fig6_prep,
  query               = seu_prep,
  refdata             = list(paper_cluster = "clustering_paper"),
  reference.reduction = "pca",
  reduction.model     = "umap"
)

# 4) Add meta factors for NP/GDM x immature/mature, if available
if ("CellType_l2" %in% colnames(seu_mapped@meta.data)) {
  seu_mapped$Type <- ifelse(
    seu_mapped$CellType_l2 %in% c("LD-DEFA3", "LD-OLFM4", "LD-OLR1", "LD-MMP9", "LD-S100A4"),
    "Immature",
    "Mature"
  )
}
if ("Sample_Name" %in% colnames(seu_mapped@meta.data)) {
  seu_mapped$Group <- ifelse(
    seu_mapped$Sample_Name %in% c("GDM_HD", "GDM_LD"),
    "GDM",
    "NP"
  )
}
if (all(c("Type", "Group") %in% colnames(seu_mapped@meta.data))) {
  seu_mapped$Group_Type <- paste(seu_mapped$Group, seu_mapped$Type, sep = "_")
}

# 5) Reference vs query on reference UMAP
ref_levels <- levels(Idents(fig6_prep))

seu_mapped$predicted.paper_cluster <- factor(
  seu_mapped$predicted.paper_cluster,
  levels = ref_levels
)

# Fig.6-like palette (clusters 1..27)
pal_fig6 <- c(
  "1"  = "#5BA6D1",
  "2"  = "#40C2C7",
  "3"  = "#8E83B8",
  "4"  = "#2EB9B2",
  "5"  = "#2AA6B0",
  "6"  = "#B79FC9",
  "7"  = "#7569A9",
  "8"  = "#60C1A3",
  "9"  = "#7FC97A",
  "10" = "#4FA7D4",
  "11" = "#57B95F",
  "12" = "#5CB6C0",
  "13" = "#CF7DAE",
  "14" = "#3E86C2",
  "15" = "#4DBE49",
  "16" = "#CFA939",
  "17" = "#E28A2B",
  "18" = "#D96227",
  "19" = "#E07479",
  "20" = "#79BE69",
  "21" = "#D8B546",
  "22" = "#E2C65F",
  "23" = "#E6B24A",
  "24" = "#6FB05A",
  "25" = "#C57BB8",
  "26" = "#A9719E",
  "27" = "#E6A139"
)
pal_use <- pal_fig6[levels(Idents(fig6_prep))]

p_ref <- DimPlot(
  fig6_prep,
  reduction = "umap",
  group.by  = "clustering_paper",
  label     = TRUE,
  cols      = pal_use,
  raster    = FALSE
) + coord_fixed(1)

p_qry <- DimPlot(
  seu_mapped,
  reduction = "ref.umap",
  group.by  = "predicted.paper_cluster",
  label     = TRUE,
  cols      = pal_use,
  raster    = FALSE
) + coord_fixed(1)

ref_umap_cluster <- p_ref + p_qry

ggsave(
  ref_umap_cluster,
  filename = file.path(dir_figures, "fig14_ref_umap_cluster.png"),
  width    = 500,
  height   = 300,
  units    = "mm",
  dpi      = 300
)
ggsave(
  ref_umap_cluster,
  filename = file.path(dir_figures, "fig14_ref_umap_cluster.pdf"),
  width    = 500,
  height   = 300,
  units    = "mm",
  dpi      = 300,
  device   = "pdf",
  bg       = "transparent"
)

# 6) Reference vs query split by group/sample
fig6$group <- factor(
  fig6$Tissue,
  levels = c("BM_Healthy", "BM_HSC-T", "PB_Healthy", "G-CSF", "PB_HSC-T", "PB_PDAC")
)

p_ref_split <- DimPlot(
  fig6,
  reduction = "umap",
  group.by  = "clustering_paper",
  split.by  = "group",
  cols      = pal_use,
  label     = FALSE,
  raster    = FALSE
) + coord_fixed(1)

ggsave(
  p_ref_split,
  filename = file.path(dir_figures, "fig15_ref_umap_cluster_split.png"),
  width    = 500,
  height   = 150,
  units    = "mm",
  dpi      = 300
)
ggsave(
  p_ref_split,
  filename = file.path(dir_figures, "fig15_ref_umap_cluster_split.pdf"),
  width    = 500,
  height   = 150,
  units    = "mm",
  dpi      = 300,
  device   = "pdf",
  bg       = "transparent"
)

p_qry_split <- DimPlot(
  seu_mapped,
  reduction = "ref.umap",
  group.by  = "predicted.paper_cluster",
  split.by  = "Sample_Name",
  cols      = pal_use,
  label     = FALSE,
  raster    = FALSE
) + coord_fixed(1)

ggsave(
  p_qry_split,
  filename = file.path(dir_figures, "fig16_qry_umap_cluster_split.png"),
  width    = 500,
  height   = 150,
  units    = "mm",
  dpi      = 300
)
ggsave(
  p_qry_split,
  filename = file.path(dir_figures, "fig16_qry_umap_cluster_split.pdf"),
  width    = 500,
  height   = 150,
  units    = "mm",
  dpi      = 300,
  device   = "pdf",
  bg       = "transparent"
)

# 7) Query clusters on reference UMAP
cols_celltype <- c(
  "LD-DEFA3"  = "#FA9645",
  "LD-OLFM4"  = "#4fc2dc",
  "LD-OLR1"   = "#F684C1",
  "LD-MMP9"   = "#AA8ED6",
  "LD-S100A4" = "#EB6F5D",
  "HD-NFKBIA" = "#e0c52e",
  "HD-CXCL1"  = "#4cb842",
  "HD-CXCR4"  = "#8C564B"
)

p_qry_celltype <- DimPlot(
  seu_mapped,
  reduction = "ref.umap",
  group.by  = "CellType_l2",
  cols      = cols_celltype,
  label     = FALSE,
  raster    = FALSE
) + coord_fixed(1)

ggsave(
  p_qry_celltype,
  filename = file.path(dir_figures, "fig17_qry_umap_our_cluster.png"),
  width    = 200,
  height   = 150,
  units    = "mm",
  dpi      = 300
)
ggsave(
  p_qry_celltype,
  filename = file.path(dir_figures, "fig17_qry_umap_our_cluster.pdf"),
  width    = 200,
  height   = 150,
  units    = "mm",
  dpi      = 300,
  device   = "pdf",
  bg       = "transparent"
)

# 8) Cluster abundance (reference + query)
source_cluster <- tibble::tibble(
  Cluster = c(fig6$clustering_paper, seu_mapped$predicted.paper_cluster),
  Group   = c(fig6$group,           seu_mapped$Sample_Name)
) %>%
  dplyr::group_by(Cluster, Group) %>%
  dplyr::summarise(no.cell = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  ) %>%
  dplyr::select(Cluster, Group, perc)

p_abundance <- ggplot(
  source_cluster,
  aes(x = Group, y = perc, fill = Cluster)
) +
  geom_col(colour = "black") +
  scale_fill_manual(values = pal_use) +
  coord_fixed(ratio = 1 / 20) +
  theme_bw() +
  xlab("") +
  ylab("%")

ggsave(
  p_abundance,
  filename = file.path(dir_figures, "fig18_ref_cluster_abundance.png"),
  width    = 500,
  height   = 150,
  units    = "mm",
  dpi      = 300
)
ggsave(
  p_abundance,
  filename = file.path(dir_figures, "fig18_ref_cluster_abundance.pdf"),
  width    = 500,
  height   = 150,
  units    = "mm",
  dpi      = 300,
  device   = "pdf",
  bg       = "transparent"
)

# Export cluster abundance
wb_ref <- wb_workbook()
wb_ref <- wb_add_worksheet(wb_ref, sheet = "Ref_mapping_cluster")
wb_ref <- wb_add_data(wb_ref, sheet = "Ref_mapping_cluster", x = as.data.frame(source_cluster))
wb_save(wb_ref, file.path(dir_results, "step7_ref_mapping.xlsx"), overwrite = TRUE)

# Save mapped objects
dir.create(dirname(path_seu_mapped),     showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(path_fig6_processed), showWarnings = FALSE, recursive = TRUE)
saveRDS(seu_mapped, file = path_seu_mapped)
saveRDS(fig6_prep,  file = path_fig6_processed)

## ---------------------------------------------------------------------------
## Section 2: NP vs GDM within LD and within HD
## ---------------------------------------------------------------------------

# Prepare meta factors on original object
if (!"group" %in% colnames(seu@meta.data)) {
  seu$group <- factor(dplyr::recode(
    seu$Sample_Name,
    "NP_HD"  = "NP",
    "NP_LD"  = "NP",
    "GDM_HD" = "GDM",
    "GDM_LD" = "GDM"
  ), levels = c("NP", "GDM"))
}
if (!"gran_den" %in% colnames(seu@meta.data)) {
  seu$gran_den <- factor(dplyr::recode(
    seu$Sample_Name,
    "NP_HD"  = "HD",
    "GDM_HD" = "HD",
    "NP_LD"  = "LD",
    "GDM_LD" = "LD"
  ), levels = c("LD", "HD"))
}

DefaultAssay(seu) <- "SCT"

# Helper: clean feature set for DEG
hist_genes <- grep("Hist", rownames(seu@assays$SCT@counts), value = TRUE)
hb_genes   <- grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$SCT@counts), value = TRUE)
ensg_genes <- grep("^ENSG[0-9]+", rownames(seu@assays$SCT@counts), value = TRUE)
bad_features <- unique(c(
  hist_genes,
  hb_genes,
  ensg_genes,
  grep(
    "^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
    rownames(seu@assays$SCT@counts),
    value = TRUE
  )
))
clean_features <- setdiff(rownames(seu[["SCT"]]@data), bad_features)

# Safe wrapper for fgsea
run_fgsea_safe <- function(ranks, category, subcategory = NULL) {
  ranks <- ranks[is.finite(ranks)]
  if (length(ranks) < 20 || is.null(names(ranks))) return(NULL)

  mdf <- msigdbr(
    species    = "Homo sapiens",
    category   = category,
    subcategory = subcategory
  )
  if (!nrow(mdf)) return(NULL)

  pathways <- split(mdf$gene_symbol, mdf$gs_name)
  if (!length(pathways)) return(NULL)

  res <- tryCatch(
    {
      suppressWarnings(
        fgsea(
          pathways,
          ranks,
          minSize = 10,
          maxSize = 500,
          nproc   = 1
        )
      )
    },
    error = function(e) {
      message("fgsea error: ", e$message)
      NULL
    }
  )
  if (is.null(res)) return(NULL)

  res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ", "), "")
  res[order(res$padj, -abs(res$NES)), ]
}

## -----------------------
## A) LD: NP vs GDM
## -----------------------
cat("\n=== LD: NP vs GDM ===\n")

ld <- subset(seu, subset = gran_den == "LD" & group %in% c("NP", "GDM"))
stopifnot(ncol(ld) > 0)
tab_ld <- table(ld$group)
print(tab_ld)
if (length(tab_ld) < 2) stop("LD contains only one group (NP or GDM); cannot compare.")

Idents(ld) <- "group"
DefaultAssay(ld) <- "SCT"

deg_ld_sct <- FindMarkers(
  ld,
  assay           = "SCT",
  test.use        = "wilcox",
  ident.1         = "GDM",
  ident.2         = "NP",
  features        = clean_features,
  logfc.threshold = 0,
  densify         = FALSE
)
deg_ld_sct <- deg_ld_sct[order(deg_ld_sct$p_val_adj, -abs(deg_ld_sct$avg_log2FC)), ]

DefaultAssay(ld) <- "SCT"
universe_ld <- rownames(ld[["SCT"]]@data)

up_ld   <- rownames(subset(deg_ld_sct, p_val_adj <= 0.05 & avg_log2FC >= 0.25))
down_ld <- rownames(subset(deg_ld_sct, p_val_adj <= 0.05 & avg_log2FC <= -0.25))

ego_up_ld <- if (length(intersect(up_ld, universe_ld)) >= 10) enrichGO(
  gene          = intersect(up_ld, universe_ld),
  universe      = universe_ld,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  minGSSize     = 50,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  readable      = TRUE
) else NULL

ego_dn_ld <- if (length(intersect(down_ld, universe_ld)) >= 10) enrichGO(
  gene          = intersect(down_ld, universe_ld),
  universe      = universe_ld,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  minGSSize     = 50,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  readable      = TRUE
) else NULL

ranks_ld <- deg_ld_sct$avg_log2FC
names(ranks_ld) <- rownames(deg_ld_sct)
ranks_ld <- sort(ranks_ld[is.finite(ranks_ld)], decreasing = TRUE)

gsea_h_ld    <- run_fgsea_safe(ranks_ld, "H")
gsea_keg_ld  <- run_fgsea_safe(ranks_ld, "C2", "CP:KEGG")
gsea_rea_ld  <- run_fgsea_safe(ranks_ld, "C2", "CP:REACTOME")
gsea_gobp_ld <- run_fgsea_safe(ranks_ld, "C5", "GO:BP")

# Save LD results
out_ld <- file.path(dir_results, "step7_GDMvsNP_LD_major.xlsx")
wb_ld <- wb_workbook()
wb_ld <- wb_add_worksheet(wb_ld, "DEG_SCT")
wb_ld <- wb_add_data(wb_ld, "DEG_SCT", data.frame(gene = rownames(deg_ld_sct), deg_ld_sct))

if (!is.null(ego_up_ld)) {
  wb_ld <- wb_add_worksheet(wb_ld, "GO_BP_up")
  wb_ld <- wb_add_data(wb_ld, "GO_BP_up", as.data.frame(ego_up_ld))
}
if (!is.null(ego_dn_ld)) {
  wb_ld <- wb_add_worksheet(wb_ld, "GO_BP_down")
  wb_ld <- wb_add_data(wb_ld, "GO_BP_down", as.data.frame(ego_dn_ld))
}
if (!is.null(gsea_h_ld)) {
  wb_ld <- wb_add_worksheet(wb_ld, "GSEA_Hallmark")
  wb_ld <- wb_add_data(wb_ld, "GSEA_Hallmark", as.data.frame(gsea_h_ld))
}
if (!is.null(gsea_keg_ld)) {
  wb_ld <- wb_add_worksheet(wb_ld, "GSEA_KEGG")
  wb_ld <- wb_add_data(wb_ld, "GSEA_KEGG", as.data.frame(gsea_keg_ld))
}
if (!is.null(gsea_rea_ld)) {
  wb_ld <- wb_add_worksheet(wb_ld, "GSEA_REACTOME")
  wb_ld <- wb_add_data(wb_ld, "GSEA_REACTOME", as.data.frame(gsea_rea_ld))
}
if (!is.null(gsea_gobp_ld)) {
  wb_ld <- wb_add_worksheet(wb_ld, "GSEA_GO_BP")
  wb_ld <- wb_add_data(wb_ld, "GSEA_GO_BP", as.data.frame(gsea_gobp_ld))
}
wb_save(wb_ld, out_ld, overwrite = TRUE)

saveRDS(
  list(
    deg_sct   = deg_ld_sct,
    ego_up    = ego_up_ld,
    ego_down  = ego_dn_ld,
    gsea_h    = gsea_h_ld,
    gsea_keg  = gsea_keg_ld,
    gsea_rea  = gsea_rea_ld,
    gsea_gobp = gsea_gobp_ld
  ),
  file.path(dir_results, "step7_GDMvsNP_LD_major_DEG.rds")
)

## -----------------------
## B) HD: NP vs GDM
## -----------------------
cat("\n=== HD: NP vs GDM ===\n")

hd <- subset(seu, subset = gran_den == "HD" & group %in% c("NP", "GDM"))
stopifnot(ncol(hd) > 0)
tab_hd <- table(hd$group)
print(tab_hd)
if (length(tab_hd) < 2) stop("HD contains only one group (NP or GDM); cannot compare.")

Idents(hd) <- "group"
DefaultAssay(hd) <- "SCT"

deg_hd_sct <- FindMarkers(
  hd,
  assay           = "SCT",
  test.use        = "wilcox",
  ident.1         = "GDM",
  ident.2         = "NP",
  features        = clean_features,
  logfc.threshold = 0,
  densify         = FALSE
)
deg_hd_sct <- deg_hd_sct[order(deg_hd_sct$p_val_adj, -abs(deg_hd_sct$avg_log2FC)), ]

DefaultAssay(hd) <- "SCT"
universe_hd <- rownames(hd[["SCT"]]@data)

up_hd   <- rownames(subset(deg_hd_sct, p_val_adj <= 0.05 & avg_log2FC >= 0.25))
down_hd <- rownames(subset(deg_hd_sct, p_val_adj <= 0.05 & avg_log2FC <= -0.25))

ego_up_hd <- if (length(intersect(up_hd, universe_hd)) >= 10) enrichGO(
  gene          = intersect(up_hd, universe_hd),
  universe      = universe_hd,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  minGSSize     = 50,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  readable      = TRUE
) else NULL

ego_dn_hd <- if (length(intersect(down_hd, universe_hd)) >= 10) enrichGO(
  gene          = intersect(down_hd, universe_hd),
  universe      = universe_hd,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  minGSSize     = 50,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  readable      = TRUE
) else NULL

ranks_hd <- deg_hd_sct$avg_log2FC
names(ranks_hd) <- rownames(deg_hd_sct)
ranks_hd <- sort(ranks_hd[is.finite(ranks_hd)], decreasing = TRUE)

gsea_h_hd    <- run_fgsea_safe(ranks_hd, "H")
gsea_keg_hd  <- run_fgsea_safe(ranks_hd, "C2", "CP:KEGG")
gsea_rea_hd  <- run_fgsea_safe(ranks_hd, "C2", "CP:REACTOME")
gsea_gobp_hd <- run_fgsea_safe(ranks_hd, "C5", "GO:BP")

# Save HD results
out_hd <- file.path(dir_results, "step7_GDMvsNP_HD_major.xlsx")
wb_hd <- wb_workbook()
wb_hd <- wb_add_worksheet(wb_hd, "DEG_SCT")
wb_hd <- wb_add_data(wb_hd, "DEG_SCT", data.frame(gene = rownames(deg_hd_sct), deg_hd_sct))

if (!is.null(ego_up_hd)) {
  wb_hd <- wb_add_worksheet(wb_hd, "GO_BP_up")
  wb_hd <- wb_add_data(wb_hd, "GO_BP_up", as.data.frame(ego_up_hd))
}
if (!is.null(ego_dn_hd)) {
  wb_hd <- wb_add_worksheet(wb_hd, "GO_BP_down")
  wb_hd <- wb_add_data(wb_hd, "GO_BP_down", as.data.frame(ego_dn_hd))
}
if (!is.null(gsea_h_hd)) {
  wb_hd <- wb_add_worksheet(wb_hd, "GSEA_Hallmark")
  wb_hd <- wb_add_data(wb_hd, "GSEA_Hallmark", as.data.frame(gsea_h_hd))
}
if (!is.null(gsea_keg_hd)) {
  wb_hd <- wb_add_worksheet(wb_hd, "GSEA_KEGG")
  wb_hd <- wb_add_data(wb_hd, "GSEA_KEGG", as.data.frame(gsea_keg_hd))
}
if (!is.null(gsea_rea_hd)) {
  wb_hd <- wb_add_worksheet(wb_hd, "GSEA_REACTOME")
  wb_hd <- wb_add_data(wb_hd, "GSEA_REACTOME", as.data.frame(gsea_rea_hd))
}
if (!is.null(gsea_gobp_hd)) {
  wb_hd <- wb_add_worksheet(wb_hd, "GSEA_GO_BP")
  wb_hd <- wb_add_data(wb_hd, "GSEA_GO_BP", as.data.frame(gsea_gobp_hd))
}
wb_save(wb_hd, out_hd, overwrite = TRUE)

saveRDS(
  list(
    deg_sct   = deg_hd_sct,
    ego_up    = ego_up_hd,
    ego_down  = ego_dn_hd,
    gsea_h    = gsea_h_hd,
    gsea_keg  = gsea_keg_hd,
    gsea_rea  = gsea_rea_hd,
    gsea_gobp = gsea_gobp_hd
  ),
  file.path(dir_results, "step7_GDMvsNP_HD_major_DEG.rds")
)

## ---------------------------------------------------------------------------
## Section 3: DEG plots for NP vs GDM within LD and HD
## ---------------------------------------------------------------------------

# Volcano plot helper
make_volcano_from_FindMarkers <- function(
  deg_df,
  title_txt,
  p_col        = "p_val_adj",
  fc_col       = "avg_log2FC",
  p_cut        = 0.05,
  fc_cut       = 0.25,
  n_label_sig  = 15,
  n_label_updn = 10
) {
  stopifnot(all(c(p_col, fc_col) %in% colnames(deg_df)))
  dd <- deg_df %>%
    as.data.frame() %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::filter(!is.na(.data[[p_col]]), !is.na(.data[[fc_col]]))

  lab_sig <- dd %>%
    dplyr::arrange(.data[[p_col]]) %>%
    head(n_label_sig) %>%
    dplyr::pull(symbol)

  lab_up <- dd %>%
    dplyr::filter(.data[[p_col]] < p_cut) %>%
    dplyr::arrange(dplyr::desc(.data[[fc_col]])) %>%
    head(n_label_updn) %>%
    dplyr::pull(symbol)

  lab_down <- dd %>%
    dplyr::filter(.data[[p_col]] < p_cut) %>%
    dplyr::arrange(.data[[fc_col]]) %>%
    head(n_label_updn) %>%
    dplyr::pull(symbol)

  selectLab <- unique(c(lab_sig, lab_up, lab_down))

  EnhancedVolcano(
    dd,
    lab         = dd$symbol,
    x           = fc_col,
    y           = p_col,
    xlab        = bquote(~ Log[2] ~ "fold change"),
    title       = title_txt,
    pCutoff     = p_cut,
    FCcutoff    = fc_cut,
    pointSize   = 1.0,
    labSize     = 5.0,
    labCol      = "black",
    labFace     = "bold",
    selectLab   = selectLab,
    boxedLabels = TRUE,
    colAlpha    = 4 / 5,
    legendPosition = "top",
    legendLabSize  = 12,
    legendIconSize = 4,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    lengthConnectors = 2,
    arrowheads     = FALSE,
    colConnectors  = "black"
  )
}

plots_vol <- list()
if (exists("deg_ld_sct") && nrow(deg_ld_sct) > 0) {
  plots_vol$LD <- make_volcano_from_FindMarkers(deg_ld_sct, "LD: GDM vs NP (SCT)")
}
if (exists("deg_hd_sct") && nrow(deg_hd_sct) > 0) {
  plots_vol$HD <- make_volcano_from_FindMarkers(deg_hd_sct, "HD: GDM vs NP (SCT)")
}

if (length(plots_vol) > 0) {
  volcano_combined <- wrap_plots(plots_vol, ncol = length(plots_vol))
  ggsave(
    volcano_combined,
    filename = file.path(dir_figures, "fig19_volcano_GDMvsNP_LD_HD.png"),
    width    = 360,
    height   = 250,
    units    = "mm",
    dpi      = 300
  )
  ggsave(
    volcano_combined,
    filename = file.path(dir_figures, "fig19_volcano_GDMvsNP_LD_HD.pdf"),
    width    = 360,
    height   = 250,
    units    = "mm",
    dpi      = 300,
    device   = "pdf",
    bg       = "transparent"
  )
}

# GSEA GO:BP barplot helper
make_gsea_gobp_bar <- function(
  gsea_df,
  title_txt,
  top_n  = 5,
  col_up = "#CF6A17",
  col_dn = "#1B9E77"
) {
  if (is.null(gsea_df) || !nrow(gsea_df)) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(paste0(title_txt, "\n(no significant terms)"))
    )
  }

  df <- as.data.frame(gsea_df) %>%
    dplyr::filter(!is.na(padj), padj < 0.05)

  if (!nrow(df)) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(paste0(title_txt, "\n(no padj < 0.05)"))
    )
  }

  up_terms <- df %>%
    dplyr::arrange(dplyr::desc(NES)) %>%
    dplyr::slice(1:min(top_n, n()))

  down_terms <- df %>%
    dplyr::arrange(NES) %>%
    dplyr::slice(1:min(top_n, n()))

  plot_df <- dplyr::bind_rows(down_terms, up_terms) %>%
    dplyr::mutate(
      Direction = ifelse(NES > 0, "Enriched in GDM", "Enriched in NP"),
      pathway   = stringr::str_replace_all(pathway, "_", " "),
      pathway   = stringr::str_wrap(pathway, width = 45)
    ) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(pathway = factor(pathway, levels = unique(pathway)))

  ggplot(plot_df, aes(x = NES, y = pathway, fill = Direction)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
    scale_fill_manual(values = c("Enriched in NP" = col_dn, "Enriched in GDM" = col_up)) +
    labs(title = title_txt, x = "NES (GSEA, GO:BP)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      axis.text.y     = element_text(size = 10)
    )
}

bar_ld <- make_gsea_gobp_bar(gsea_gobp_ld, "LD: GSEA GO:BP (GDM vs NP, top 5)")
ggsave(
  bar_ld,
  filename = file.path(dir_figures, "fig20_GSEA_GOBP_bar_LD.png"),
  width    = 250,
  height   = 150,
  units    = "mm",
  dpi      = 300
)
ggsave(
  bar_ld,
  filename = file.path(dir_figures, "fig20_GSEA_GOBP_bar_LD.pdf"),
  width    = 250,
  height   = 150,
  units    = "mm",
  dpi      = 300,
  device   = "pdf",
  bg       = "transparent"
)

bar_hd <- make_gsea_gobp_bar(gsea_gobp_hd, "HD: GSEA GO:BP (GDM vs NP, all significant)")
ggsave(
  bar_hd,
  filename = file.path(dir_figures, "fig21_GSEA_GOBP_bar_HD.png"),
  width    = 250,
  height   = 100,
  units    = "mm",
  dpi      = 300
)
ggsave(
  bar_hd,
  filename = file.path(dir_figures, "fig21_GSEA_GOBP_bar_HD.pdf"),
  width    = 250,
  height   = 100,
  units    = "mm",
  dpi      = 300,
  device   = "pdf",
  bg       = "transparent"
)

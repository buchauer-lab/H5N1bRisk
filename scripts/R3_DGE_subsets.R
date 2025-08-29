library("tidyverse")
library("DESeq2")
library("EnhancedVolcano")
library("cowplot")

## SETUP #######################################################################
# Load Data --------------------------------------------------------------------
projdir <- file.path("/data/cephfs-1/work/groups/buchauer/h5n1risk")
outdir <- file.path("/data/cephfs-1/work/groups/buchauer/tmp_h5n1risk/") # adapted due to writing limitation in original project
setwd(projdir)

cts <- readRDS(file.path(projdir, "data/cts_combined_filtered.rds"))
coldata <- readRDS(file.path(projdir, "data/coldata.rds"))

coldata <- coldata[colnames(cts), ] # match row/col order of data and metadata
all(rownames(coldata) == colnames(cts)) # check if that worked 👆

# optionally subset to omit viral reads
viral_genes <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
cts <- cts[!rownames(cts) %in% viral_genes, ]

# DESeq ########################################################################
# Setup ========================================================================
# Settings
plim <- 0.05 # limit for adjusted p value
lfclim <- 1.5 # limit for relevant effect size (after lfc shrinking)
lfclim <- log2(lfclim) # fold changes are handled in log2 in deseq
shrink_type = "ashr"

res <- list() # collect all results here

# function to extract and consistently save results
get_deseq_res <- function(dds, contrast, filename = "condition_test-ctrl", 
                          subtitle = NULL, plim, lfclim, shrink_type) {

  res1 <- results(dds, contrast = contrast, alpha=plim, lfcThreshold = lfclim) %>% 
    .[order(.$padj),]
  

  
  shr <- lfcShrink(dds, contrast = contrast, type=shrink_type) %>% 
    .[order(.$padj),]
  
  vulc <- EnhancedVolcano(shr, 
                          lab = rownames(shr),  
                          x = 'log2FoldChange',
                          y = 'padj',
                          ylab = bquote(~-Log[10] ~ italic(P~adj)),
                          pCutoff = plim,
                          FCcutoff = lfclim,
                          # title is filename after the underscore with spaces around the hyphen
                          title = gsub("-", " - ", sub("^[^_]*_", "", filename)), 
                          subtitle = sub("_.*", "", filename),
                          caption = NULL,
                          legendLabels= c("Not Sig.", 
                                          bquote(Log[2] ~ FC > "|" * .(2^lfclim) * "|"),
                                          bquote(adj.~pvalue < .(plim)),
                                          bquote(adj.~pvalue ~ and ~ Log[2] ~ FC)),
                          legendPosition = 'right'
  )

  
  # Save all DGE results for running gsea later
  shr %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(desc(abs(log2FoldChange))) %>%
    write.csv(file.path(outdir, "data/gsea/deg", paste0(filename, "_all_deg.csv")
    ), row.names = FALSE)
  
  # Save only significant DEGs for Analysis
  shr %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(desc(abs(log2FoldChange))) %>%
    filter(padj < plim) %>%
    write.csv(file.path(outdir, "data/dge", paste0(filename, "_sig_deg.csv")
    ), row.names = FALSE)
  
  ggsave(filename = file.path(outdir, "data/dge", paste0(filename, "_volcano.png")),
         plot = vulc, width = 12, height = 7, units = "cm", scale = 2)
  
  return(list(res = res1, res_shrunk = shr, volcano = vulc))
}

# All Types ====================================================================
# !!! This creates the full DDS object, but runs on the subsets only
# analysis on the full DDS object is done in R3_DGE_all

# Setup ------------------------------------------------------------------------
# constructing the full DESeq dataset including all types and interaction 
dds_full <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ type + treatment + type:treatment
)
dds_full$treatment <- relevel(dds_full$treatment, ref = "control")
dds_full$type <- relevel(dds_full$type, ref = "HuLu")

# ALI ==========================================================================
# Setup ------------------------------------------------------------------------
# Subsetting the dds to one type, allows us to use a simpler design to exclude 
# treatment:type interactions and donor effects from the analysis
dds_ali <- dds_full[, dds_full$type == "ALI"] # Subset to only ALI 

dds_ali$type <- droplevels(dds_ali$type)  # Clean up unused factor levels
dds_ali$donor <- droplevels(dds_ali$donor)
dds_ali$treatment <- droplevels(dds_ali$treatment)

dds_ali <- DESeqDataSet(dds_ali, design = ~ donor + treatment)
dds_ali <- DESeq(dds_ali)

# get the model matrix to define contrasts
mod_mat_ali <- model.matrix(design(dds_ali), colData(dds_ali))

# define coefficient vectors for each group
coeffs_ali <- list()
coeffs_ali$A_SH <- colMeans(mod_mat_ali[dds_ali$treatment == "A_SH", ])
coeffs_ali$A_Tx <- colMeans(mod_mat_ali[dds_ali$treatment == "A_Tx", ])
coeffs_ali$KAN_1 <- colMeans(mod_mat_ali[dds_ali$treatment == "KAN_1", ])
coeffs_ali$PanIII <- colMeans(mod_mat_ali[dds_ali$treatment == "PanIII", ])
coeffs_ali$control <- colMeans(mod_mat_ali[dds_ali$treatment == "control", ])

# Results - All Virus vs Control -----------------------------------------------
res$atx_ctrl$ali <- get_deseq_res(
  dds_ali, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_ali$A_Tx - coeffs_ali$control,
  filename = "ALI_A_Tx-control")

res$ash_ctrl$ali <- get_deseq_res(
  dds_ali, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_ali$A_SH - coeffs_ali$control,
  filename = "ALI_A_SH-control")

res$kan_ctrl$ali <- get_deseq_res(
  dds_ali, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_ali$KAN_1 - coeffs_ali$control,
  filename = "ALI_KAN_1-control")

res$pan_ctrl$ali <- get_deseq_res(
  dds_ali, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_ali$PanIII - coeffs_ali$control,
  filename = "ALI_PanIII-control")

# Results - A/Tx vs All Others -------------------------------------------------
res$atx_ash$ali <- get_deseq_res(
  dds_ali, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_ali$A_Tx - coeffs_ali$A_SH,
  filename = "ALI_A_Tx-A_SH")

res$atx_kan$ali <- get_deseq_res(
  dds_ali, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_ali$A_Tx - coeffs_ali$KAN_1,
  filename = "ALI_A_Tx-KAN_1")

res$atx_pan$ali <- get_deseq_res(
  dds_ali, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_ali$A_Tx - coeffs_ali$PanIII,
  filename = "ALI_A_Tx-PanIII")


# bLO ==========================================================================
# Setup ------------------------------------------------------------------------
# Subsetting the dds to one type, allows us to use a simpler design to exclude 
# treatment:type interactions and donor effects from the analysis
dds_blo <- dds_full[, dds_full$type == "bLO"] # Subset to only bLO 

dds_blo$type <- droplevels(dds_blo$type)  # Clean up unused factor levels
dds_blo$donor <- droplevels(dds_blo$donor)  
#dds_blo$treatment <- droplevels(dds_blo$treatment)  

dds_blo <- DESeqDataSet(dds_blo, design = ~ donor + treatment)
#dds_blo$sizeFactor <- NULL
dds_blo <- DESeq(dds_blo)

# get the model matrix to define contrasts
mod_mat_blo <- model.matrix(design(dds_blo), colData(dds_blo))

# define coefficient vectors for each group
coeffs_blo <- list()
coeffs_blo$A_SH <- colMeans(mod_mat_blo[dds_blo$treatment == "A_SH", ])
coeffs_blo$A_Tx <- colMeans(mod_mat_blo[dds_blo$treatment == "A_Tx", ])
coeffs_blo$KAN_1 <- colMeans(mod_mat_blo[dds_blo$treatment == "KAN_1", ])
coeffs_blo$PanIII <- colMeans(mod_mat_blo[dds_blo$treatment == "PanIII", ])
coeffs_blo$control <- colMeans(mod_mat_blo[dds_blo$treatment == "control", ])

# Results - All Virus vs Control -----------------------------------------------
res$atx_ctrl$blo <- get_deseq_res(
  dds_blo, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_blo$A_Tx - coeffs_blo$control,
  filename = "bLO_A_Tx-control")

res$ash_ctrl$blo <- get_deseq_res(
  dds_blo, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_blo$A_SH - coeffs_blo$control,
  filename = "bLO_A_SH-control")

res$kan_ctrl$blo <- get_deseq_res(
  dds_blo, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_blo$KAN_1 - coeffs_blo$control,
  filename = "bLO_KAN_1-control")

res$pan_ctrl$blo <- get_deseq_res(
  dds_blo, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_blo$PanIII - coeffs_blo$control,
  filename = "bLO_PanIII-control")

# Results - A/Tx vs All Others -------------------------------------------------
res$atx_ash$blo <- get_deseq_res(
  dds_blo, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_blo$A_Tx - coeffs_blo$A_SH,
  filename = "bLO_A_Tx-A_SH")

res$atx_kan$blo <- get_deseq_res(
  dds_blo, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_blo$A_Tx - coeffs_blo$KAN_1,
  filename = "bLO_A_Tx-KAN_1")

res$atx_pan$blo <- get_deseq_res(
  dds_blo, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_blo$A_Tx - coeffs_blo$PanIII,
  filename = "bLO_A_Tx-PanIII")


# HuLu ==========================================================================
# Setup ------------------------------------------------------------------------
# Subsetting the dds to one type, allows us to use a simpler design to exclude 
# treatment:type interactions and donor effects from the analysis
dds_hulu <- dds_full[, dds_full$type == "HuLu"] # Subset to only HuLu 

dds_hulu$type <- droplevels(dds_hulu$type)  # Clean up unused factor levels
dds_hulu$donor <- droplevels(dds_hulu$donor)  

dds_hulu <- DESeqDataSet(dds_hulu, design = ~ donor + treatment)
#dds_hulu$sizeFactor <- NULL
dds_hulu <- DESeq(dds_hulu)

# get the model matrix to define contrasts
mod_mat_hulu <- model.matrix(design(dds_hulu), colData(dds_hulu))

# define coefficient vectors for each group
coeffs_hulu <- list()
coeffs_hulu$A_SH <- colMeans(mod_mat_hulu[dds_hulu$treatment == "A_SH", ])
coeffs_hulu$A_Tx <- colMeans(mod_mat_hulu[dds_hulu$treatment == "A_Tx", ])
coeffs_hulu$KAN_1 <- colMeans(mod_mat_hulu[dds_hulu$treatment == "KAN_1", ])
coeffs_hulu$PanIII <- colMeans(mod_mat_hulu[dds_hulu$treatment == "PanIII", ])
coeffs_hulu$control <- colMeans(mod_mat_hulu[dds_hulu$treatment == "control", ])

# Results - All Virus vs Control -----------------------------------------------
res$atx_ctrl$hulu <- get_deseq_res(
  dds_hulu, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_hulu$A_Tx - coeffs_hulu$control,
  filename = "HuLu_A_Tx-control")

res$ash_ctrl$hulu <- get_deseq_res(
  dds_hulu, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_hulu$A_SH - coeffs_hulu$control,
  filename = "HuLu_A_SH-control")

res$kan_ctrl$hulu <- get_deseq_res(
  dds_hulu, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_hulu$KAN_1 - coeffs_hulu$control,
  filename = "HuLu_KAN_1-control")

res$pan_ctrl$hulu <- get_deseq_res(
  dds_hulu, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_hulu$PanIII - coeffs_hulu$control,
  filename = "HuLu_PanIII-control")

# Results - A/Tx vs All Others -------------------------------------------------
res$atx_ash$hulu <- get_deseq_res(
  dds_hulu, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_hulu$A_Tx - coeffs_hulu$A_SH,
  filename = "HuLu_A_Tx-A_SH")

res$atx_kan$hulu <- get_deseq_res(
  dds_hulu, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_hulu$A_Tx - coeffs_hulu$KAN_1,
  filename = "HuLu_A_Tx-KAN_1")

res$atx_pan$hulu <- get_deseq_res(
  dds_hulu, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_hulu$A_Tx - coeffs_hulu$PanIII,
  filename = "HuLu_A_Tx-PanIII")

# ==============================================================================
saveRDS(res, file = file.path(outdir, "data", "dge_results.rds"))
res <- readRDS(file.path(outdir, "data", "dge_results.rds"))


# Create Composite Volcano Plots ===============================================
leg <- get_plot_component(res$atx_ctrl$full$volcano +
                            theme(legend.position = "bottom",
                                  legend.direction = "horizontal",
                                  legend.margin = margin(10,20,10,20)), #t,r,b,l
                          "guide-box-bottom")

# A/Tx vs Control --------------------------------------------------------------
title <- ggdraw() + 
  draw_label("A/Tx vs Control", fontface = 'bold', x = 0, hjust = 0, size = 18)+
  theme(plot.margin = margin(0, 0, 0, 15))


p_atx_ctrl <- plot_grid(title, 
          plot_grid(#res$atx_ctrl$full$volcano + ggtitle(NULL) +
                     # theme(legend.position="none",),
                           res$atx_ctrl$ali$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"),
                           res$atx_ctrl$blo$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"),
                           res$atx_ctrl$hulu$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"), 
                    ncol = 2, align = "vh"),
          leg, ncol = 1, rel_heights = c(0.1, 1.5, 0.1))

# A/SH vs Control --------------------------------------------------------------
title <- ggdraw() + 
  draw_label("A/SH vs Control", fontface = 'bold', x = 0, hjust = 0, size = 18)+
  theme(plot.margin = margin(0, 0, 0, 15))


p_ash_ctrl <- plot_grid(title, 
          plot_grid(#res$ash_ctrl$full$volcano + ggtitle(NULL) +
                     # theme(legend.position="none",),
                    res$ash_ctrl$ali$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"),
                    res$ash_ctrl$blo$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"),
                    res$ash_ctrl$hulu$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"), 
                    ncol = 2, align = "vh"),
          leg, ncol = 1, rel_heights = c(0.1, 1.5, 0.1))

# KAN-1 vs Control --------------------------------------------------------------
title <- ggdraw() + 
  draw_label("KAN-1 vs Control", fontface = 'bold', x = 0, hjust = 0, size = 18)+
  theme(plot.margin = margin(0, 0, 0, 15))


p_kan_ctrl <- plot_grid(title, 
          plot_grid(#res$kan_ctrl$full$volcano + ggtitle(NULL) +
                    #  theme(legend.position="none",),
                    res$kan_ctrl$ali$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"),
                    res$kan_ctrl$blo$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"),
                    res$kan_ctrl$hulu$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"), 
                    ncol = 2, align = "vh"),
          leg, ncol = 1, rel_heights = c(0.1, 1.5, 0.1))

# PanIII vs Control --------------------------------------------------------------
title <- ggdraw() + 
  draw_label("PanIII vs Control", fontface = 'bold', x = 0, hjust = 0, size = 18)+
  theme(plot.margin = margin(0, 0, 0, 15))


p_pan_ctrl <- plot_grid(title, 
          plot_grid(#res$pan_ctrl$full$volcano + ggtitle(NULL) +
                    #  theme(legend.position="none",),
                    res$pan_ctrl$ali$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"),
                    res$pan_ctrl$blo$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"),
                    res$pan_ctrl$hulu$volcano + ggtitle(NULL) + 
                      theme(legend.position="none"), 
                    ncol = 2, align = "vh"),
          leg, ncol = 1, rel_heights = c(0.1, 1.5, 0.1))

# A/Tx vs A/SH --------------------------------------------------------------
title <- ggdraw() + 
  draw_label("A/Tx vs A/SH", fontface = 'bold', x = 0, hjust = 0, size = 18)+
  theme(plot.margin = margin(0, 0, 0, 15))


p_atx_ash <- plot_grid(title, 
                        plot_grid(#res$atx_ash$full$volcano + ggtitle(NULL) +
                                  #  theme(legend.position="none",),
                                  res$atx_ash$ali$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"),
                                  res$atx_ash$blo$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"),
                                  res$atx_ash$hulu$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"), 
                                  ncol = 2, align = "vh"),
                        leg, ncol = 1, rel_heights = c(0.1, 1.5, 0.1))

# A/Tx vs KAN-1 ----------------------------------------------------------------
title <- ggdraw() + 
  draw_label("A/Tx vs KAN-1", fontface = 'bold', x = 0, hjust = 0, size = 18)+
  theme(plot.margin = margin(0, 0, 0, 15))


p_atx_kan <- plot_grid(title, 
                        plot_grid(#res$atx_kan$full$volcano + ggtitle(NULL) +
                                  #  theme(legend.position="none",),
                                  res$atx_kan$ali$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"),
                                  res$atx_kan$blo$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"),
                                  res$atx_kan$hulu$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"), 
                                  ncol = 2, align = "vh"),
                        leg, ncol = 1, rel_heights = c(0.1, 1.5, 0.1))

# A/Tx vs PanIII --------------------------------------------------------------
title <- ggdraw() + 
  draw_label("A/Tx vs PanIII", fontface = 'bold', x = 0, hjust = 0, size = 18)+
  theme(plot.margin = margin(0, 0, 0, 15))


p_atx_pan <- plot_grid(title, 
                        plot_grid(#res$atx_pan$full$volcano + ggtitle(NULL) +
                                  #  theme(legend.position="none",),
                                  res$atx_pan$ali$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"),
                                  res$atx_pan$blo$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"),
                                  res$atx_pan$hulu$volcano + ggtitle(NULL) + 
                                    theme(legend.position="none"), 
                                  ncol = 2, align = "vh"),
                        leg, ncol = 1, rel_heights = c(0.1, 1.5, 0.1))

# Save Plots -------------------------------------------------------------------
ggsave(filename = file.path(outdir, "data/dge", "A_Tx_ctrl_grid_volcano.png"),
       plot = p_atx_ctrl, width = 12, height = 12, units = "cm", scale = 2)
ggsave(filename = file.path(outdir, "data/dge", "A_SH_ctrl_grid_volcano.png"),
       plot = p_ash_ctrl, width = 12, height = 12, units = "cm", scale = 2)
ggsave(filename = file.path(outdir, "data/dge", "KAN_1_ctrl_grid_volcano.png"),
       plot = p_kan_ctrl, width = 12, height = 12, units = "cm", scale = 2)
ggsave(filename = file.path(outdir, "data/dge", "PanIII_ctrl_grid_volcano.png"),
       plot = p_pan_ctrl, width = 12, height = 12, units = "cm", scale = 2)
ggsave(filename = file.path(outdir, "data/dge", "A_Tx_A_SH_grid_volcano.png"),
       plot = p_atx_ash, width = 12, height = 12, units = "cm", scale = 2)
ggsave(filename = file.path(outdir, "data/dge", "A_Tx_KAN_1_grid_volcano.png"),
       plot = p_atx_kan, width = 12, height = 12, units = "cm", scale = 2)
ggsave(filename = file.path(outdir, "data/dge", "A_Tx_PanIII_grid_volcano.png"),
       plot = p_atx_pan, width = 12, height = 12, units = "cm", scale = 2)




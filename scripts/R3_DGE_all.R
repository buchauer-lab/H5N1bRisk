library("tidyverse")
library("DESeq2")
library("EnhancedVolcano")
library("cowplot")

## SETUP #######################################################################
# Load Data --------------------------------------------------------------------
projdir <- file.path("/data/cephfs-1/work/groups/buchauer/h5n1risk")
outdir <- file.path(projdir, "data/dge")
setwd(projdir)

cts <- readRDS(file.path(projdir, "data/cts_combined_filtered.rds"))
coldata <- readRDS(file.path(projdir, "data/coldata.rds"))

coldata <- coldata[colnames(cts), ] # match row/col order of data and metadata
all(rownames(coldata) == colnames(cts)) # check if that worked 👆

# optionally subset to omit viral reads
viral_genes <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
cts <- cts[!rownames(cts) %in% viral_genes, ]

# Compute counts per million ###################################################
lib_sizes <- colSums(cts)
cpm_matrix <- t(t(cts) / lib_sizes) * 1e6
write.csv(cpm_matrix, "data/CPM_without_viral_genes.csv")


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
  # Note: We continue with shr from here on instead of res1
  
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
# Setup ------------------------------------------------------------------------
# constructing the full DESeq dataset including all types and interaction 
dds_full <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ type + treatment + type:treatment
)
dds_full$treatment <- relevel(dds_full$treatment, ref = "control")
dds_full$type <- relevel(dds_full$type, ref = "HuLu")

dds_full <- DESeq(dds_full)


# work with model matrix based contrasts, see here:
# https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
# This allows to easily account for sample imbalances (ALI has less samples) 
# and lets us define the desired contrasts in a (relatively) comprehensible way 

# get the model matrix to define contrasts
mod_mat_full <- model.matrix(design(dds_full), colData(dds_full))

# define coefficient vectors for each group
coeffs_full <- list()
coeffs_full$A_SH <- colMeans(mod_mat_full[dds_full$treatment == "A_SH", ])
coeffs_full$A_Tx <- colMeans(mod_mat_full[dds_full$treatment == "A_Tx", ])
coeffs_full$KAN_1 <- colMeans(mod_mat_full[dds_full$treatment == "KAN_1", ])
coeffs_full$PanIII <- colMeans(mod_mat_full[dds_full$treatment == "PanIII", ])
coeffs_full$control <- colMeans(mod_mat_full[dds_full$treatment == "control", ])

# for the comparison of different tissues in A/Tx:
coeffs_full$A_Tx_ALI   <- colMeans(mod_mat_full[dds_full$treatment == "A_Tx" & dds_full$type == "ALI", ])
coeffs_full$A_Tx_bLO   <- colMeans(mod_mat_full[dds_full$treatment == "A_Tx" & dds_full$type == "bLO", ])
coeffs_full$A_Tx_HuLu  <- colMeans(mod_mat_full[dds_full$treatment == "A_Tx" & dds_full$type == "HuLu", ])

coeffs_full$control_ALI  <- colMeans(mod_mat_full[dds_full$treatment == "control" & dds_full$type == "ALI", ])
coeffs_full$control_bLO  <- colMeans(mod_mat_full[dds_full$treatment == "control" & dds_full$type == "bLO", ])
coeffs_full$control_HuLu <- colMeans(mod_mat_full[dds_full$treatment == "control" & dds_full$type == "HuLu", ])

# Results - All Virus vs Control -----------------------------------------------
res$atx_ctrl$full <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_full$A_Tx - coeffs_full$control,
  filename = "AllTypes_A_Tx-control")

res$ash_ctrl$full <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_full$A_SH - coeffs_full$control,
  filename = "AllTypes_A_SH-control")

res$kan_ctrl$full <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_full$KAN_1 - coeffs_full$control,
  filename = "AllTypes_KAN_1-control")

res$pan_ctrl$full <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_full$PanIII - coeffs_full$control,
  filename = "AllTypes_PanIII-control")

# Results - A/Tx vs All Others -------------------------------------------------
res$atx_ash$full <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_full$A_Tx - coeffs_full$A_SH,
  filename = "AllTypes_A_Tx-A_SH")

res$atx_kan$full <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_full$A_Tx - coeffs_full$KAN_1,
  filename = "AllTypes_A_Tx-KAN_1")

res$atx_pan$full <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = coeffs_full$A_Tx - coeffs_full$PanIII,
  filename = "AllTypes_A_Tx-PanIII")

# Results - A/Tx Type Comparison -----------------------------------------------
# Create the contrasts in this pattern:
# (A_Tx_ALI - A_Tx_bLO) - (control_ALI - control_bLO)
# this extracts the differences between ALI and bLO that we see in A_Tx and 
# removes the effects that we also observe in control

res$atx$ali_blo <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = (coeffs_full$A_Tx_ALI - coeffs_full$A_Tx_bLO) - 
    (coeffs_full$control_ALI - coeffs_full$control_bLO),
  filename = "A_Tx_ALI-bLO")

res$atx$ali_hulu <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = (coeffs_full$A_Tx_ALI - coeffs_full$A_Tx_HuLu) - 
    (coeffs_full$control_ALI - coeffs_full$control_HuLu),
  filename = "A_Tx_ALI-HuLu")

res$atx$blo_hulu <- get_deseq_res(
  dds_full, plim = plim, lfclim = lfclim, shrink_type = shrink_type,
  contrast = (coeffs_full$A_Tx_bLO - coeffs_full$A_Tx_HuLu) - 
    (coeffs_full$control_bLO - coeffs_full$control_HuLu),
  filename = "A_Tx_bLO-HuLu")


# fix Volcano Titles & Subtitles
res$atx$ali_blo$volcano$labels$title <- "ALI - bLO"
res$atx$ali_hulu$volcano$labels$title <- "ALI - HuLu"
res$atx$blo_hulu$volcano$labels$title <- "bLO - HuLu"
res$atx$ali_blo$volcano$labels$subtitle <- "A/Tx effect, Control effect removed"
res$atx$ali_hulu$volcano$labels$subtitle <- "A/Tx effect, Control effect removed"
res$atx$blo_hulu$volcano$labels$subtitle <- "A/Tx effect, Control effect removed"
ggsave(filename = file.path(outdir, "data/dge", "A_Tx_ALI-bLO_volcano.png"),
       plot = res$atx$ali_blo$volcano, width = 12, height = 7, 
       units = "cm", scale = 2)
ggsave(filename = file.path(outdir, "data/dge", "A_Tx_ALI-HuLu_volcano.png"),
       plot = res$atx$ali_hulu$volcano, width = 12, height = 7, 
       units = "cm", scale = 2)
ggsave(filename = file.path(outdir, "data/dge", "A_Tx_bLO-HuLu_volcano.png"),
       plot = res$atx$blo_hulu$volcano, width = 12, height = 7, 
       units = "cm", scale = 2)

# The analysis on subsets is done in R3_DGE_subsets.R
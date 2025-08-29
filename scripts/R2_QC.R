library("tidyverse")
library("broom")
library("DESeq2")
library("pheatmap")
library("cowplot")
library("msigdbr")

# Load Data --------------------------------------------------------------------
projdir <- file.path("/data/cephfs-1/work/groups/buchauer/h5n1risk")
figdir <- file.path(projdir, "data/qc/")
setwd(projdir)

cts_path   <- file.path(projdir, "data", "cts_combined_filtered.rds")
cts        <- readRDS(cts_path)
fileprefix <- if (grepl("proteincoding", cts_path)) "prot_" else ""

coldata <- readRDS(file.path(projdir, "data/coldata.rds"))

coldata <- coldata[colnames(cts), ] # match row/col order of data and metadata
all(rownames(coldata) == colnames(cts)) # check if that worked 👆

# construct DESeqDataSets ------------------------------------------------------
run_deseq <- function(cts, coldata, designFormula,filter_fn = NULL) {
  # subset like this:
  # type == "HuLu" and treatment == "control"
  # --> filter_fn = expr(type == "HuLu" & treatment == "control"))
  
  if (!is.null(filter_fn)) {
    coldata_sub <- coldata %>% filter(!!filter_fn)
  } else {
    coldata_sub <- coldata
  }
  cts_sub <- cts[, rownames(coldata_sub)]
  
  dds <- DESeqDataSetFromMatrix(
    countData = cts_sub,
    colData = coldata_sub,
    design = designFormula
  )
  return(DESeq(dds))
}

# All tissues, All viruses -----------------------------------------------------
dds <- run_deseq(cts, coldata, ~ type + treatment)
dds$treatment <- relevel(dds$treatment, ref = "control")
dds$type <- relevel(dds$type, ref = "HuLu")
vsd <- vst(dds, blind = FALSE)

# HuLu -------------------------------------------------------------------------
hulu_dds <- run_deseq(cts, coldata,
                      ~ donor + treatment,
                      expr(type == "HuLu"))
hulu_vsd <- vst(hulu_dds, blind = FALSE)

# ALI --------------------------------------------------------------------------
ali_dds <- run_deseq(cts, coldata,
                     ~ donor + treatment,
                     expr(type == "ALI"))
ali_vsd <- vst(ali_dds, blind = FALSE)

# bLO --------------------------------------------------------------------------
blo_dds <- run_deseq(cts, coldata,
                     ~ donor + treatment,
                     expr(type == "bLO"))
blo_vsd <- vst(blo_dds, blind = FALSE)


### QC #########################################################################
# # Plot counts of "Most Significant" Gene to give Overview
# plotCounts(dds[dds$type == "HuLu",], gene="RPE65", intgroup=c("type","treatment"))
# plotCounts(dds, gene="RPE65", intgroup=c("type","treatment"))


viral_genes <-c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

count_bars_df <- data.frame(
  sample = colnames(dds),
  counts = colSums(counts(dds)),
  vircounts = colSums(counts(dds)[rownames(counts(dds)) %in% viral_genes, ]),
  type = colData(dds)$type,
  treatment = colData(dds)$treatment,
  donor = colData(dds)$donor ) %>%
  left_join(read_csv(file.path(projdir, "readcounts.csv")), by = "sample")

count_bars_df <- count_bars_df[order(count_bars_df$type, count_bars_df$treatment), ]
count_bars_df$sample <- factor(count_bars_df$sample, levels = count_bars_df$sample)

# Total Counts
p_tot_cts <- ggplot(count_bars_df, aes(x = sample, y = counts, fill = treatment)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ type, scales = "free_x", nrow = 1) +  # Grouping by 'type'
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Total Counts per Sample",
    x = "Sample",
    y = "Total Counts",
    fill = "Treatment"
  )
p_tot_cts

# Viral Counts
p_vir_cts <- ggplot(count_bars_df, aes(x = sample, y = vircounts, fill = treatment)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ type, scales = "free_x", nrow = 1) +  # Grouping by 'type'
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Viral Counts per Sample",
    x = "Sample",
    y = "Total Counts",
    fill = "Treatment"
  )
p_vir_cts

# Viral Counts (log)
p_vir_log <- ggplot(count_bars_df, aes(x = sample, y = vircounts, fill = treatment)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ type, scales = "free_x", nrow = 1) +  # Grouping by 'type'
  scale_y_continuous(trans='log10') +
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Viral Counts per Sample (Log Scaled!)",
    x = "Sample",
    y = "Total Counts",
    fill = "Treatment"
  )


# Sequencing Read Counts
seq_reads <- count_bars_df %>%
  pivot_longer(
    cols = c(run1, run2),
    names_to = "run",
    values_to = "reads"
  ) %>%
  mutate(run = factor(run, levels = c("run2", "run1")))

p_seq_cts <- ggplot(seq_reads, aes(x = sample, y = reads, fill = treatment, alpha = run)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ type, scales = "free_x", nrow = 1) +  # Group by 'type'
  scale_alpha_manual(values = c("run1" = 1.0, "run2" = 0.6)) +
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Total Sequencing Reads per Sample",
    x = "Sample",
    y = "Read Counts",
    fill = "Treatment",
    alpha = "Run"
  ) 
p_seq_cts


# Heatmap of count Matrix
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:25]
df <- as.data.frame(colData(dds)[, c("type", "treatment")])
df <- df[order(df$type, df$treatment), ] # Order columns by type and treatment
mat <- assay(vsd)[select, rownames(df)] # Subset the vsd matrix accordingly
p_heat_cts <- pheatmap(mat, cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=df,
         main = "Genes with Highest Normalized Expression")
p_heat_cts

# Sample Distances
sampleDists <- dist(vir(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$type, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$treatment, vsd$type, sep="-")
p_heat_dist <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,)
p_heat_dist


# PCA --------------------------------------------------------------------------
options(repr.plot.width=12, repr.plot.height=7) # plot size (in cm?)
# PCA - All tissues
pcaData <- plotPCA(vsd, intgroup=c("treatment", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_all_pca <- ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=type)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(15, 16, 17, 18, 3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("PCA - All Tissues") +
  coord_fixed() +
  theme_minimal_grid(12)


# PCA - HuLu
pcaData <- plotPCA(hulu_vsd, intgroup=c("treatment", "donor"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_hulu_pca <- ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=donor)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(15, 16, 17, 18, 3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("PCA - HuLu") +
  coord_fixed() +
  theme_minimal_grid(12)


# PCA - ALI
pcaData <- plotPCA(ali_vsd, intgroup=c("treatment", "donor"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_ali_pca <- ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=donor)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(15, 16, 17, 18, 3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("PCA - ALI") +
  coord_fixed() +
  theme_minimal_grid(12)


# PCA - bLO
pcaData <- plotPCA(blo_vsd, intgroup=c("treatment", "donor"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_blo_pca <- ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=donor)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(15, 16, 17, 18, 3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("PCA - bLO") +
  coord_fixed() +
  theme_minimal_grid(12)
p_blo_pca
  

p_combined_pca <- plot_grid(p_all_pca, p_hulu_pca, p_blo_pca, p_ali_pca)


# Save Plots -------------------------------------------------------------------
ggsave(file.path(figdir, paste0(fileprefix, "viral_counts.png")), plot = p_vir_cts, 
       width = 12, height = 7, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "viral_log_counts.png")), plot = p_vir_log, 
       width = 12, height = 7, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "total_counts.png")), plot = p_tot_cts, 
       width = 12, height = 7, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "seq_read_counts.png")), plot = p_seq_cts, 
       width = 12, height = 7, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "heat_counts.png")), plot = p_heat_cts, 
       width = 12, height = 7, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "heat_dists.png")), plot = p_heat_dist, 
       width = 12, height = 7, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "pca_combined.png")), plot = p_combined_pca, 
       width = 12, height = 12, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "pca_hulu.png")), plot = p_hulu_pca, 
       width = 12, height = 7, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "pca_ali.png")), plot = p_ali_pca, 
       width = 12, height = 7, units = "cm",scale = 2)
ggsave(file.path(figdir, paste0(fileprefix, "pca_blo.png")), plot = p_blo_pca, 
       width = 12, height = 7, units = "cm",scale = 2)



### Correlation of Viral Reads to Interferon Gamma expression ##################
## Setup =======================================================================
# scoring function -------------------------------------------------------------
gene_set_scores <- function(mat, genes) {
  missing <- setdiff(genes, rownames(mat))
  if (length(missing) > 0) {
    warning("Dropping ", length(missing), " genes not found in matrix: ",
            paste(missing, collapse = ", "))
    genes <- intersect(genes, rownames(mat))
  }
  # keep only genes in the set, drop = FALSE to preserve matrix structure
  submat <- mat[genes, , drop = FALSE]
  
  mean <- colMeans(submat)
  median <- apply(submat, 2, median)
  # z-score, then average across genes
  zmean <- colSums(t(scale(t(submat))))/length(genes) 

  data.frame(
    sample      = colnames(mat),
    mean_score  = mean,
    median_score = median,
    zmean_score  = zmean,
    row.names   = NULL,
    check.names = FALSE
    )
}

# normalized/ variance stabilized counts --------------------------------------- 
# norm_cts <- counts(dds, normalized = TRUE)
vst_cts <- assay(vsd) # use this!

# gene sets of interest --------------------------------------------------------

# interferon gamma response (only in HuLu because immune cell dependent)
ifng_genes <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  dplyr::pull(gene_symbol)

il1b_genes <- msigdbr(species = "Homo sapiens", 
                         category = "C2", subcollection = "CP:REACTOME"
                         ) %>% 
  filter(gs_name == "REACTOME_INTERLEUKIN_1_SIGNALING") %>% 
  dplyr::pull(gene_symbol) %>%
  unique()


viral_genes <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

## Results =====================================================================
all_scores <- coldata %>%
  rownames_to_column("sample") %>%                            
  mutate( # scores using ALL samples
    full_viral_score = gene_set_scores(vst_cts,  viral_genes)$zmean_score,
    full_ifng_score  = gene_set_scores(vst_cts,  ifng_genes )$zmean_score,
    full_il1b_score  = gene_set_scores(vst_cts,  il1b_genes )$zmean_score
  ) %>% 
  group_by(type) %>%                                          
  mutate( # scores for each tissue individually
    type_viral_score = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                       viral_genes)$zmean_score,
    type_ifng_score  = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                       ifng_genes )$zmean_score,
    type_il1b_score  = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                       il1b_genes )$zmean_score
  ) %>% 
  group_by(treatment) %>%
  mutate( # scores for each tissue individually
    treatment_viral_score = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                     viral_genes)$zmean_score,
    treatment_ifng_score  = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                   ifng_genes )$zmean_score,
    treatment_il1b_score  = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                            il1b_genes )$zmean_score
    ) %>% 
  group_by(type, treatment) %>%
  mutate( # scores for each tissue individually
    group_viral_score = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                            viral_genes)$zmean_score,
    group_ifng_score  = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                            ifng_genes )$zmean_score,
    group_il1b_score  = gene_set_scores(vst_cts[ , sample, drop = FALSE], 
                                        il1b_genes )$zmean_score
  ) %>% 
  ungroup() 

# Overall ----------------------------------------------------------------------
cor_all <- cor.test(all_scores$full_viral_score, all_scores$full_il1b_score,
                     method = "spearman")

p_cor_all <- ggplot(all_scores, aes(x = full_viral_score, y = full_il1b_score))+
  geom_point(size = 2, aes(colour = treatment, shape = type)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  labs(x = "Viral Genes mean z‑score",
       y = "REACTOME IL1 Signaling mean z‑score",
       title = "Correlation Viral Genes and IL1 Signaling (all samples)",
       subtitle = paste0("Spearman ρ = ",
                      round(cor_all$estimate, 2),
                      ",  p = ",
                      signif(cor_all$p.value, 3))) +
  theme_minimal_grid()

ggsave(file.path(figdir, paste0(fileprefix, "corr_vir-il1b_all.png")), plot = p_cor_all, 
       width = 10, height = 7, units = "cm",scale = 2)
# print(cor_all)
p_cor_all


# By Type ----------------------------------------------------------------------
cors_by_type <- all_scores %>%
  group_by(type) %>%
  summarise( tidy(
    cor.test(full_viral_score, full_il1b_score, method = "spearman")
  ) )
cors_by_type

t <- "HuLu" # ALI bLO HuLu
p_cor_type <- all_scores %>% 
  filter(type == t) %>%                   
  ggplot(aes(x = type_viral_score,
             y = type_il1b_score)) +           
  geom_point(size = 2, aes(colour = treatment)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  labs(
    x = "Viral Genes mean z‑score",
    y = "REACTOME IL1 Signaling mean z‑score",
    title = str_glue("Correlation Viral Genes and IL1 Signaling ({t} only)"),
    subtitle = paste0("Spearman ρ = ",
                      round(cors_by_type$estimate[cors_by_type$type == t], 2),
                      ",  p = ",
                      signif(cors_by_type$p.value[cors_by_type$type == t], 3))
  ) +
  theme_minimal_grid()

ggsave(file.path(figdir, paste0(fileprefix, "corr_vir-il1b_",t ,".png")), plot = p_cor_type, 
       width = 10, height = 7, units = "cm",scale = 2)
p_cor_type


# By Treatment --------------------------------------------------------------
cors_by_treatment <- all_scores %>%
  group_by(treatment) %>%
  summarise( tidy(
    cor.test(treatment_viral_score, treatment_il1b_score, method = "spearman")
  ) )
cors_by_treatment

vir <- "KAN_1" # A_SH A_Tx control KAN_1 PanIII
p_cor_treatment <- all_scores %>%
  filter(treatment == vir) %>%
  ggplot(aes(x = treatment_viral_score,
             y = treatment_il1b_score)) +
  geom_point(size = 2, aes(colour = type)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  labs(
    x = "Viral Genes mean z‑score",
    y = "REACTOME IL1 Signaling mean z‑score",
    title = str_glue("Correlation Viral Genes and IL1 Signaling ({vir} only)"),
    subtitle = paste0("Spearman ρ = ",
                      round(cors_by_treatment$estimate[cors_by_treatment$treatment == vir], 2),
                      ",  p = ",
                      signif(cors_by_treatment$p.value[cors_by_treatment$treatment == vir], 3))
  ) +
  theme_minimal_grid()

ggsave(file.path(figdir, paste0(fileprefix, "corr_vir-il1b_",vir ,".png")), plot = p_cor_treatment, 
       width = 10, height = 7, units = "cm",scale = 2)
p_cor_treatment


cors_by_treatment









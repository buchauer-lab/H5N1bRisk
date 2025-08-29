library("tidyverse")
library("biomaRt")
library("rtracklayer")

# Samples:
# ctrl: no Virus
# atx: H5N1b
# ash: H5N8
# kan: H5N1
# pan: H3N2

projdir <- file.path("/data/cephfs-1/work/groups/buchauer/h5n1risk/")
setwd(file.path(projdir, "/data/"))

# Metadata #####################################################################
coldata <- read.csv(file.path(projdir, "h5n1meta.csv"),
                    header = TRUE, row.names = "sample") %>%
  dplyr::select(genome, type, treatment, donor) %>%
  mutate(type = recode(type,"bLO-derived layer" = "bLO")) %>%
  mutate(across(where(is.character), ~ gsub("[-/ ]", "_", .))) %>%
  mutate(across(everything(), as.factor))

ali_samples <- rownames(coldata[coldata$type == "ALI",])
blo_samples <- rownames(coldata[coldata$type == "bLO",])
hulu_samples <- rownames(coldata[coldata$type == "HuLu",])

# Create Count Matrix ##########################################################
# combined human and viral counts ----------------------------------------------
cts <- purrr::reduce(list(  
  rownames_to_column(read.table("counts/h5n1b_vir_GRCh38_counts.txt", 
                                header=TRUE, na.strings = c("","N/A","NaN"),
                                row.names="Geneid")[,-1:-5], "geneid"), 
  rownames_to_column(read.table("counts/h5n8_vir_GRCh38_counts.txt", 
                                header=TRUE, na.strings = c("","N/A","NaN"),
                                row.names="Geneid")[,-1:-5], "geneid"), 
  rownames_to_column(read.table("counts/h5n1_vir_GRCh38_counts.txt", 
                                header=TRUE, na.strings = c("","N/A","NaN"),
                                row.names="Geneid")[,-1:-5], "geneid"), 
  rownames_to_column(read.table("counts/h3n2_vir_GRCh38_counts.txt", 
                                header=TRUE, na.strings = c("","N/A","NaN"),
                                row.names="Geneid")[,-1:-5], "geneid")), 
  full_join, by="geneid") %>%  
  column_to_rownames("geneid") %>%
  as.matrix() %>% 
  { 
    . <- .[rowSums(. >= 10) >= 3, ] # Prefilter min 10 counts for 3 samples
    . <- .[ , sort(colnames(.)) ]   # Sort columns alphabetically
    .
  }

colnames(cts) <- str_split_i(colnames(cts), '\\.', 12)


# add gene names ---------------------------------------------------------------
# create dict to map ensembl ids to gene names
HumMart <- useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl")

hgncSymbols <- data.frame(getBM(filters= "ensembl_gene_id",
                                attributes= c("ensembl_gene_id","hgnc_symbol"),
                                values=rownames(cts),
                                mart= HumMart))

# create map of genes and ids
genemap <- tibble(geneid = rownames(cts)) %>%
  left_join(hgncSymbols, join_by(geneid==ensembl_gene_id), multiple = "any") %>%
  mutate(hgnc_symbol = case_when(hgnc_symbol == '' ~ geneid,
                                 is.na(hgnc_symbol) ~ geneid,
                                 .default = hgnc_symbol))


# extract the duplicates for documentation
dups_ids <- genemap[
  genemap$hgnc_symbol %in% genemap$hgnc_symbol[duplicated(genemap$hgnc_symbol)], 
                    "geneid"]

cts_dups <- cts[rownames(cts) %in% dups_ids$geneid, ] %>%
  as.data.frame() %>%
  rownames_to_column("geneid") %>%
  left_join(genemap, by = "geneid") %>%
  dplyr::select(geneid, hgnc_symbol, everything())

# Assign gene names
all(rownames(cts) == genemap$geneid)
rownames(cts) <- genemap$hgnc_symbol[match(rownames(cts), genemap$geneid)]

# check size and duplicates before merging duplicates
dim(cts)
unique(rownames(cts)[duplicated(rownames(cts))])

# Group by gene name and sum across rows with same name
cts <- cts %>%
  as.data.frame(cts) %>%
  mutate(gene=rownames(cts)) %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames("gene") %>%
  as.matrix() 

# check size and duplicates after merging
dim(cts)
unique(rownames(cts)[duplicated(rownames(cts))])


# Extract Protein Coding Genes #################################################
# also keep viral gene segments
viral_genes <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

# Import GTF file
gtf <- import(file.path(projdir, "genomes/individualGenomes/Homo_sapiens.GRCh38.113.gtf"))

# Filter to keep only genes
genes <- gtf[gtf$type == "gene"]

# Look at available metadata columns
mcols(genes)

# Filter protein-coding genes
protein_coding_genes <- genes[genes$gene_biotype == "protein_coding" ]

# Extract gene IDs
protein_coding_ids <- unique(protein_coding_genes$gene_id)

# keep only viral and protein coding genes
keepers <- union(protein_coding_genes$gene_id, viral_genes)
keepers <- genemap[genemap$geneid %in% keepers,]# have to go through ensembl ids
cts_prot <- cts[rownames(cts) %in% keepers$hgnc_symbol,]


# Stricter Prefiltering ########################################################
cts <- readRDS(file.path(projdir, "data/cts_combined_raw.rds"))
cts_pc <- readRDS(file.path(projdir, "data/cts_combined_proteincoding_raw.rds"))

ali_samples <- rownames(coldata[coldata$type == "ALI",])
blo_samples <- rownames(coldata[coldata$type == "bLO",])
hulu_samples <- rownames(coldata[coldata$type == "HuLu",])

cts_filtered <- cts %>% 
  as_tibble(rownames= NA) %>%
  filter(
    rowSums(dplyr::select(., all_of(ali_samples)))   >= 10,   # total ALI counts ≥ 10
    rowSums(dplyr::select(., all_of(blo_samples)))   >= 10,   # total bLO counts ≥ 10
    rowSums(dplyr::select(., all_of(hulu_samples)))  >= 10    # total HuLu counts ≥ 10
  ) %>%
  as.matrix()
dim(cts)
dim(cts_filtered)

cts_pc_filtered <- cts_pc %>% 
  as_tibble(rownames= NA) %>%
  filter(
    rowSums(dplyr::select(., all_of(ali_samples)))   >= 10,   # total ALI counts ≥ 10
    rowSums(dplyr::select(., all_of(blo_samples)))   >= 10,   # total bLO counts ≥ 10
    rowSums(dplyr::select(., all_of(hulu_samples)))  >= 10    # total HuLu counts ≥ 10
  ) %>%
  as.matrix()
dim(cts_pc)
dim(cts_pc_filtered)




# Save files  ------------------------------------------------------------------
saveRDS(cts, file.path(projdir, "data/cts_combined_raw.rds"))
saveRDS(cts_prot, file.path(projdir, "data/cts_combined_proteincoding_raw.rds"))
saveRDS(cts_filtered, file.path(projdir, "data/cts_combined_filtered.rds"))
saveRDS(cts_pc_filtered, file.path(projdir, "data/cts_combined_proteincoding_filtered.rds"))
saveRDS(genemap, file.path(projdir, "data/genenames.rds"))
write_csv(cts_dups, file.path(projdir, "data/qc/genename_duplicates.csv"))
saveRDS(coldata, file.path(projdir, "data/coldata.rds"))


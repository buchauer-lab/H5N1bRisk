library("tidyverse")
library("DESeq2")
library("EnhancedVolcano")
library("cowplot")

## SETUP #######################################################################
# Load Data --------------------------------------------------------------------
projdir <- file.path("/data/cephfs-1/work/groups/buchauer/h5n1risk")
setwd(projdir)

cts <- readRDS(file.path(projdir, "data/cts_combined_proteincoding_filtered.rds"))
coldata <- readRDS(file.path(projdir, "data/coldata.rds"))

coldata <- coldata[colnames(cts), ] # match row/col order of data and metadata
all(rownames(coldata) == colnames(cts)) # check if that worked 👆

# optionally subset to omit viral reads
viral_genes <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
cts <- cts[!rownames(cts) %in% viral_genes, ]


# Define Functions -------------------------------------------------------------
# compute fold change between two conditions
# sample2 supposed to be control
get_FC <- function(df, sample1, sample2){
  return(df[[sample1]]/df[[sample2]])
}

# basic approach
donors <- unique(coldata[['donor']])
result_df <- data.frame(matrix(nrow = nrow(cts), ncol = 4*length(donors)))

colnames(result_df) <- paste(rep(donors, each=4), rep(c("A_Tx", "A_SH", "KAN_1", "PanIII"), 4), sep = "__")
rownames(result_df) <- rownames(cts)

for (donor in donors) {
  for (treatment in c("A_Tx", "A_SH", "KAN_1", "PanIII")) {
    colnam <- paste(donor, treatment, sep = "__")
    print(colnam)
    smpl1 <- rownames(coldata[(coldata$donor == donor) & (coldata$treatment == treatment),])
    smpl2 <- rownames(coldata[(coldata$donor == donor) & (coldata$treatment == "control"),])
    print(smpl1)
    print(smpl2)
    result_df[,colnam] <- cts[,smpl1] / cts[,smpl2]
  }
}

log_df <- log2(result_df)

write.csv(result_df, "data/FC_donors.csv")
write.csv(log_df, "data/logFC_donors.csv")

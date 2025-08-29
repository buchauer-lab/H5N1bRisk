# BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library("tidyverse")
library("clusterProfiler")
library("enrichplot")
library("cowplot")
library("org.Hs.eg.db", character.only = TRUE)


## SETUP #######################################################################
# Load Data --------------------------------------------------------------------
projdir <- file.path("/data/cephfs-1/work/groups/buchauer/h5n1risk")
outdir <- file.path(projdir, "data/gsea")
indir <- file.path(outdir, "deg")
setwd(projdir)
filenames = dir(indir, pattern="*.csv")
genemap <- read_rds(file.path(projdir, "data/genenames.rds"))


for( i in 1:length(filenames)){
  file <- filenames[i]
  name <- gsub("_all_deg.csv","", file)
  print(paste(i,name, sep = ": "))
  degs <- read_csv(file.path(outdir, "deg", file), show_col_types = FALSE) %>%
    as_tibble() %>%
    left_join(genemap, by = join_by(gene == hgnc_symbol)) %>% # return ensembl ids
    dplyr::select(geneid, log2FoldChange) %>%
    na.omit() %>%
    arrange(desc(log2FoldChange)) %>%
    deframe()
  
  print("GO")
  # GO #########################################################################
  gsea_GO <- gseGO(degs,
                   ont = "ALL",
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   minGSSize = 25,
                   maxGSSize = 500,
                   verbose = FALSE
                   )
  
  write_csv(as.data.frame(gsea_GO), 
            file.path(outdir, paste0(name, "_gseaGO.csv")))

  
  # probably more useful in top down case:
  # gsea_GO <- pairwise_termsim(gsea_GO) # necessary for some cluster vis
  # emapplot(simpleGO, showCategory = 20)
  
  
  if (!is.null(gsea_GO) && nrow(gsea_GO@result) > 0) {
    gsea_GO <- simplify(gsea_GO)
    
    p_dot_GO <- dotplot(gsea_GO, showCategory = 20, 
            label_format = 80,
            title = paste0(name, " - GSEA (GO)") , split=".sign") + 
      facet_grid(.~.sign)
    
    ggsave(filename = file.path(outdir, paste0(name, "_dotplot_gseaGO.png")),  
           plot = p_dot_GO, width = 12, height = 7, units = "cm", scale = 2)
    
  }
  

 
  # KEGG #######################################################################
  print("KEGG")
  # gene IDs need to be converted for gseKEGG
  ids<-bitr(names(degs), fromType = "ENSEMBL", 
            toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  # remove duplicate IDS 
  ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  keggdegs <- degs %>%
    tibble(ENSEMBL = names(.), log2FoldChange = .) %>% 
    left_join(ids) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    dplyr::select(ENTREZID, log2FoldChange) %>%
    na.omit() %>%
    arrange(desc(log2FoldChange)) %>%
    deframe()
  
  gsea_KEGG <- gseKEGG(keggdegs,
                       organism = "hsa",
                       keyType = "ncbi-geneid",
                       minGSSize = 25,
                       maxGSSize = 500,
                       verbose = FALSE
                       )

  
  
  if (!is.null(gsea_KEGG) && nrow(gsea_KEGG@result) > 0) {
      
    p_dot_KEGG <- dotplot(gsea_KEGG, showCategory = 20, 
                     label_format = 80,
                     title = paste0(name, " - GSEA (KEGG)") , split=".sign") + 
      facet_grid(.~.sign)
    
    ggsave(filename = file.path(outdir, paste0(name, "_dotplot_gseaKEGG.png")),  
           plot = p_dot_KEGG, width = 12, height = 7, units = "cm", scale = 2)
    
  }
  
  write_csv(as.data.frame(gsea_KEGG), 
            file.path(outdir, paste0(name, "_gseaKEGG.csv")))

}

  
  
  
  
  



All required software and packages are listed in /environment.yml

1. makegtf.ipynb:
    Create custom (minimal) GTFs from the supplied fasta files.

2. merge_seqs.sh:
    Merge sequencing runs prior to alignment

3. combinefasta.sh:
    Combine genomes by concatenating human and viral fasta.
    Check chromosome labels as sanity check.

4. combinegtf.sh:
    Combine human and viral gtf (remove header of second file).

5. star_buildIndex.sh:
    Create STAR indices from the combined gtf & fasta files

6. qc.sh:
    Run FastQC on sequencing data and summarize with multiQC

7. star_align_all.sh star_align_one.sh:
    align reads to one of the/ all human + virus genomes (for each Genome extract the corresponding read files from the metadata file)
    
8. star_align_summarize.sh:
    summarizes _Log.final.out files of all STAR runs

9. featureCounts.sh:
    count reads per genome using subreads featureCounts
    
10. R1_preprocess.R:
    Combine Host and Viral Count matrices and Assign gene names to the GeneIDs, produces:
      - cts_human_raw.rds & cts_viral_raw.rds: raw count matrices with only slight prefiltering
      - dds_both_geneNames.rds: DeSeq Object with combined counts of host and virus, and Gene Names assigned

11. R2_QC.R:
    Some QC Plots for readcounting and PCA
    Correlation between gene sets
    
12. R3_DGE_all.R R3_DGE_subsets.R:
    DGE Analysis comparing:
      - Each Virus vs Control:
        - All samples combined controlling for Type (ALI,bLO ,HuLu) and Type:Treatment Interaction
        - each Type (R3_DGE_subsets.R) individually controlling for Donor
      - A/Tx vs all other Viruses
        - All samples combined controlling for Type (ALI,bLO ,HuLu) and Type:Treatment Interaction
        - each Type (R3_DGE_subsets.R) individually controlling for Donor
      - A/Tx Effects between different Organoids (removing Control effects)
        
13. R4_GSEA.R:
    GSEA for all DGE comparisons

14. R5_FoldChanges_per_patient.R
    Calculate fold changes within individual Donors.

#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --nodes=1             # number of nodes to allocate
#SBATCH --ntasks=1            # number of processes (tasks)
#SBATCH --cpus-per-task=24    # cores per process/task
#SBATCH --mem=256G            # request ram 
#SBATCH --time=1-0            # max time for job (2d)
#SBATCH --output=slurmlogs/slurm-%j.out
#SBATCH --error=slurmlogs/slurm-%j.out

source /data/cephfs-1/home/users/${USER}/.bashrc # activate conda/mamba
mamba activate bulk # env with STAR

# h3n2_vir_GRCh38 h5n1_vir_GRCh38 h5n1b_vir_GRCh38 h5n8_vir_GRCh38
INPUT="h5n8_vir_GRCh38"

OUTPATH="/data/cephfs-1/work/groups/buchauer/h5n1risk/genomes/"
INPATH=${OUTPATH}"combinedGenomes/"${INPUT}

STAR --runThreadN 24 \
    --runMode genomeGenerate \
    --genomeDir ${OUTPATH}STAR_${INPUT/_vir/} \
    --genomeFastaFiles ${INPATH}.fa \
    --sjdbGTFfile ${INPATH}.gtf \
    --sjdbOverhang 149 
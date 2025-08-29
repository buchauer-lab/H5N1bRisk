#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=1-0
#SBATCH --output=slurmlogs/slurm-%j.out
#SBATCH --error=slurmlogs/slurm-%j.out

source /data/cephfs-1/home/users/${USER}/.bashrc
mamba activate bulk


# h3n2_GRCh38  h5n1_GRCh38  h5n1b_GRCh38  h5n8_GRCh38
GENOME="h5n8_GRCh38"

ROOTPATH="/data/cephfs-1/work/groups/buchauer/h5n1risk"
GENOMEPATH="${ROOTPATH}/genomes"
READPATH="${ROOTPATH}/raw"
ALIGNMENTPATH="${ROOTPATH}/alignment"
METADATA="${ROOTPATH}/h5n1meta.csv"
GENOME_DIR="${GENOMEPATH}/STAR_${GENOME}"

if [[ ! -d "$GENOME_DIR" ]]; then
    echo "Error: Genome directory ${GENOME_DIR} does not exist!" >&2
    exit 1
fi


while IFS=',' read -r sample genome_name _rest; do
    if [[ "$sample" == "sample" || "$genome_name" == "genome" || -z "$sample" || -z "$genome_name" ]]; then
        continue
    fi

    if [[ "$genome_name" == "$GENOME" ]]; then
        R1="${READPATH}/*${sample}*_R1_*.fastq.gz"
        R2="${R1/_R1_/_R2_}"

        eval R1=$(echo $R1)
        eval R2=$(echo $R2)

        if [[ ! -f "$R1" || ! -f "$R2" ]]; then
            echo "Warning: Paired files not found for sample ${sample}" >&2
            continue
        fi

        OUTDIR="${ALIGNMENTPATH}/${sample}"
        OUTPREFIX="${OUTDIR}/${sample}"

        mkdir -p "$OUTDIR"
        echo "Running STAR for ${sample}..."

        STAR --genomeDir "$GENOME_DIR" \
             --readFilesIn "$R1" "$R2" \
             --readFilesCommand zcat \
             --genomeLoad LoadAndKeep \
             --outReadsUnmapped Fastx \
             --runThreadN 24 \
             --outFileNamePrefix "${OUTPREFIX}_" \
             --limitBAMsortRAM 15000000000 \
             --outSAMtype BAM SortedByCoordinate
    fi
done < "$METADATA"

echo "Unloading STAR genome index for $GENOME..."
STAR --genomeDir "$GENOME_DIR" --genomeLoad Remove

#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --nodes=1             # number of nodes to allocate
#SBATCH --ntasks=1            # number of processes (tasks)
#SBATCH --cpus-per-task=24    # cores per process/task
#SBATCH --mem=128G            # request ram 
#SBATCH --time=1-0            # max time for job
#SBATCH --output=slurmlogs/slurm-%j.out
#SBATCH --error=slurmlogs/slurm-%j.out

source /data/cephfs-1/home/users/${USER}/.bashrc
mamba activate bulk

ROOTPATH="/data/cephfs-1/work/groups/buchauer/h5n1risk"
GENOMEPATH="${ROOTPATH}/genomes"
READPATH="${ROOTPATH}/raw"
ALIGNMENTPATH="${ROOTPATH}/alignment"
METADATA="${ROOTPATH}/h5n1meta.csv"

declare -A SAMPLE_TO_GENOME
declare -A GENOME_TO_SAMPLES

# Skip header and populate mappings
while IFS=',' read -r sample genome _rest; do
    if [[ "$sample" == "sample" || "$genome" == "genome" || -z "$sample" || -z "$genome" ]]; then
        continue
    fi
    SAMPLE_TO_GENOME["$sample"]="$genome"
    GENOME_TO_SAMPLES["$genome"]+="$sample "
done < "$METADATA"

# Process each genome group
for genome in "${!GENOME_TO_SAMPLES[@]}"; do
    echo "Loading STAR genome index for ${genome}..."
    GENOME_DIR="${GENOMEPATH}/STAR_${genome}"

    if [[ ! -d "$GENOME_DIR" ]]; then
        echo "Error: Genome directory ${GENOME_DIR} does not exist!" >&2
        exit 1
    fi

    for sample in ${GENOME_TO_SAMPLES[$genome]}; do
        R1="${READPATH}/*${sample}*_R1_*.fastq.gz"
        R2="${R1/_R1_/_R2_}"

        # Expand wildcards safely
        eval R1=$(echo $R1)
        eval R2=$(echo $R2)

        if [[ ! -f "$R1" || ! -f "$R2" ]]; then
            echo "Warning: Paired files not found for sample ${sample}" >&2
            continue
        fi

        OUTDIR="${ALIGNMENTPATH}/${sample}"
        OUTPREFIX="${OUTDIR}/${sample}"

        mkdir -p "$OUTDIR"
        echo "Running STAR for ${R1} and ${R2}..."

        STAR --genomeDir "$GENOME_DIR" \
             --readFilesIn "$R1" "$R2" \
             --readFilesCommand zcat \
             --genomeLoad LoadAndKeep \
             --outReadsUnmapped Fastx \
             --runThreadN 24 \
             --outFileNamePrefix "${OUTPREFIX}_" \
             --limitBAMsortRAM 15000000000 \
             --outSAMtype BAM SortedByCoordinate
    done

    echo "Unloading genome ${genome}..."
    STAR --genomeDir "$GENOME_DIR" --genomeLoad Remove
done

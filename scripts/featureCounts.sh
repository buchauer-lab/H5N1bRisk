#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --nodes=1             # number of nodes to allocate
#SBATCH --ntasks=1            # number of processes (tasks)
#SBATCH --cpus-per-task=16    # cores per process/task
#SBATCH --mem=32G             # request ram 
#SBATCH --time=1-0            # max time for job
#SBATCH --output=slurmlogs/slurm-%j.out
#SBATCH --error=slurmlogs/slurm-%j.out

source /data/cephfs-1/home/users/${USER}/.bashrc
mamba activate bulk

ROOTPATH="/data/cephfs-1/work/groups/buchauer/h5n1risk/"
GTFPATH="${ROOTPATH}/genomes/combinedGenomes/"
ALIGNMENTPATH="${ROOTPATH}/alignment/"
COUNTSPATH="${ROOTPATH}/data/counts/"
METADATA="${ROOTPATH}/h5n1meta.csv"

declare -A GENOME_TO_SAMPLES

# Skip header and populate mappings
while IFS=',' read -r sample genome _rest; do
    if [[ "$sample" == "sample" || "$genome" == "genome" || -z "$sample" || -z "$genome" ]]; then
        continue
    fi
    genome="${genome/_/_vir_}"
    GENOME_TO_SAMPLES["$genome"]+="$sample "
done < "$METADATA"


# Process each genome group
for genome in "${!GENOME_TO_SAMPLES[@]}"; do
    echo "Processing samples for genome: $genome"
    GTF_FILE="${GTFPATH}/${genome}.gtf"

    if [[ ! -f "$GTF_FILE" ]]; then
        echo "Error: GTF file ${GTF_FILE} does not exist!" >&2
        exit 1
    fi

    # Prepare a list of BAM files for this genome
    BAM_LIST=()
    for sample in ${GENOME_TO_SAMPLES[$genome]}; do
        BAM="${ALIGNMENTPATH}/${sample}/${sample}_Aligned.sortedByCoord.out.bam"
        if [[ -f "$BAM" ]]; then
            BAM_LIST+=("$BAM")
        else
            echo "Warning: BAM file not found for ${sample}" >&2
        fi
    done
    
    # Output file
    OUTFILE="${COUNTSPATH}/${genome}_counts.txt"
    
    echo "Running featureCounts for genome ${genome}..."
    featureCounts -T 8 -p -B -C -a "$GTF_FILE" -o "$OUTFILE" "${BAM_LIST[@]}"
    

done
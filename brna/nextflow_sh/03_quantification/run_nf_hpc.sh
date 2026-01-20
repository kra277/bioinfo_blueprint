#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=bp_brna_quant
#SBATCH --mail-type=END
#SBATCH --mail-user=kra277@nyu.edu
#SBATCH --account=torch_pr_505_general
#SBATCH --output=log/brna_quant.txt

module purge
module load nextflow/25.10.2

mkdir -p log

# Set required variables
export NF_GSE="gse312831" 
export SING_CONT="/scratch/projects/aouizeratlab/containers/2026/brna/brna_cont_v1.sif"
export NF_MAIN_DIR="$SCRATCH/bp_brna/"

echo "Running Pipeline for: $NF_MAIN_DIR / $NF_GSE"

nextflow run main.nf \
    -c nf.config \
    -profile hpc \
    --sif "$SING_CONT" \
    -resume

# Check progress squeue -u id
# To stop the process, use: scancel <JOBID> 
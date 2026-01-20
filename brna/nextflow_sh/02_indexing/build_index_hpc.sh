#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=04:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=index_build
#SBATCH --mail-type=END
#SBATCH --account=torch_pr_505_general
#SBATCH --mail-user=kra277@nyu.edu
#SBATCH --output=log/slurm_%j_log.txt

# Make sure the params refelect your file locations and your requirements
# change the num of threads based on your local computer config and needs
# set main local directory
MAIN_DIR=/scratch/kra277
# Location of the container
SING_CONT=${MAIN_DIR}/test_brna/brna_cont_v1.sif
# Location to store the index
INDEX_CONT=${MAIN_DIR}/gse236374/data/star_index
# Ensembl release
ENSEMBL_RELEASE=115
# could be homo_sapiens or mus_musculus
ORGANISM=mus_musculus
# capitalize the first letter
ORG_NAME=Mus_musculus
# for humans GRCh38, for mouse GRCm39
BUILD=GRCm39
# read length -1 -> OVER_HANG
OVER_HANG=149

ENSEMBL_BASE=https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${ORGANISM}/dna
ENSEMBL_GTF_BASE=https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${ORGANISM}
GTF_FILE=${ORG_NAME}.${BUILD}.${ENSEMBL_RELEASE}.gtf
FA_FILE=${ORG_NAME}.${BUILD}.dna.primary_assembly.fa

# Make required directories
mkdir -p log
mkdir -p ${INDEX_CONT}

# Download if missing
if [ ! -f "${INDEX_CONT}/${FA_FILE}" ]; then
  wget -P "${INDEX_CONT}" "${ENSEMBL_BASE}/${FA_FILE}.gz"
  gunzip -f "${INDEX_CONT}/${FA_FILE}.gz"
fi

if [ ! -f "${INDEX_CONT}/${GTF_FILE}" ]; then
  wget  -P "${INDEX_CONT}" "${ENSEMBL_GTF_BASE}/${GTF_FILE}.gz"
  gunzip -f "${INDEX_CONT}/${GTF_FILE}.gz"
fi

# Make Star Index
singularity exec ${SING_CONT} STAR \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --runMode genomeGenerate \
    --genomeDir ${INDEX_CONT} \
    --genomeFastaFiles ${INDEX_CONT}/${FA_FILE} \
    --sjdbGTFfile ${INDEX_CONT}/${GTF_FILE} \
    --sjdbOverhang ${OVER_HANG}
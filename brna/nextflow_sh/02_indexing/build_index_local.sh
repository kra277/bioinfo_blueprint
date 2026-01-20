#!/bin/bash
set -euo pipefail
mkdir -p "log"

# Make sure the params refelect your file locations and your requirements
# change the num of threads based on your local computer config and needs
# set main local directory
MAIN_DIR=/mnt/myvolume
# Location of the container
SING_CONT=${MAIN_DIR}/containers/2026/singularity/brna/brna_cont_v1.sif
# Location to store the index
INDEX_CONT=${MAIN_DIR}/blueprint/gse236374/data/star_index
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
# Threads to use on your workstation
THREADS=20

# ----------------------------
# Derived variables (leave as is)
# ----------------------------
ENSEMBL_BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${ORGANISM}/dna"
ENSEMBL_GTF_BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${ORGANISM}"
GTF_FILE="${ORG_NAME}.${BUILD}.${ENSEMBL_RELEASE}.gtf"
FA_FILE="${ORG_NAME}.${BUILD}.dna.primary_assembly.fa"

# ----------------------------
# Setup
# ----------------------------
mkdir -p "${INDEX_CONT}"

# ----------------------------
# Download FASTA + GTF if missing
# ----------------------------
if [ ! -f "${INDEX_CONT}/${FA_FILE}" ]; then
  wget -P "${INDEX_CONT}" "${ENSEMBL_BASE}/${FA_FILE}.gz"
  gunzip -f "${INDEX_CONT}/${FA_FILE}.gz"
fi

if [ ! -f "${INDEX_CONT}/${GTF_FILE}" ]; then
  wget -P "${INDEX_CONT}" "${ENSEMBL_GTF_BASE}/${GTF_FILE}.gz"
  gunzip -f "${INDEX_CONT}/${GTF_FILE}.gz"
fi

# ----------------------------
# Build STAR index
# ----------------------------
echo "Building STAR index..."

singularity exec \
  --bind "${MAIN_DIR}:${MAIN_DIR}" \
  "${SING_CONT}" STAR \
  --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${INDEX_CONT}" \
  --genomeFastaFiles "${INDEX_CONT}/${FA_FILE}" \
  --sjdbGTFfile "${INDEX_CONT}/${GTF_FILE}" \
  --sjdbOverhang "${OVER_HANG}" \
  > "log/star_index_build.log" 2>&1

echo "Done. STAR index is in: ${INDEX_CONT}"
echo "Log: log/star_index_build.log"
#!/bin/bash

mkdir -p log
logfile="log/brna_quant.txt"

# Set params via environment variables
# set the local main directory where the gseXXXXX is located
# Make sure the gse number is correct
# Set the container path
export NF_GSE="gse312831" 
export SING_CONT="/scratch/projects/aouizeratlab/containers/2026/brna/brna_cont_v1.sif"
export NF_MAIN_DIR="/mnt/myvolume/blueprint"

echo "Running Pipeline for: $NF_MAIN_DIR / $NF_GSE"

nohup nextflow run main.nf \
  -c nf.config \
  -profile local \
  --sif "$SING_CONT" \
  -resume \
  > "${logfile}" 2>&1 &

echo "Started Nextflow in background."
echo "Log: ${logfile}"
echo "Check progress with tail -f log/brna_quant.txt"

# To stop the process, use: kill <PID>
# To find the PID, use: ps aux | grep nextflow

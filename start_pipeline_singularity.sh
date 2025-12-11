#!/usr/bin/env bash

# ~10GB for local cache available
export TMPDIR=~/tmp # make sure this folder exists
export SINGULARITY_LOCALCADHEDIR=~/singularity_cache/

run_config_file=configs/run_template.config
projectDir=$(dirname "$0")

singularity exec ptato_1.2.0.sif /ptato/nextflow/nextflow run \
PTATO/ptato.nf \
-c ${run_config_file} \
-profile slurm -resume
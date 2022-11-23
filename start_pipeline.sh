#!/usr/bin/env bash

run_config_file=$1
projectDir=$(dirname "$0")

nextflow run \
${projectDir}/ptato.nf \
-c ${run_config_file} \
-profile slurm -resume

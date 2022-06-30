#!/usr/bin/env bash

run_config_file=$1
projectDir=$(dirname "$0")

source ${projectDir}/scripts/bash/params_check.sh

check_params ${projectDir} ${run_config_file}

nextflow run \
${projectDir}/ptato.nf \
-c ${run_config_file} \
-profile slurm -resume

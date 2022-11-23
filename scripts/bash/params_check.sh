#!/usr/bin/env bash

source ${projectDir}/scripts/bash/config_file_reader.sh
source ${projectDir}/scripts/bash/validate_path.sh
source ${projectDir}/scripts/bash/levenshtein.sh

check_params() {
  projectDir=$1
  run_config_file=$2

  expected_params=()
  for template_config_fname in ${projectDir}/configs/*.config; do
    validate_path ${projectDir} ${template_config_fname}
    p=$(read_config_file ${template_config_fname})
    expected_params+=(${p})
  done

  input_params=($(read_config_file ${run_config_file}))

  for param in "${input_params[@]}"; do
    if [[ ! " ${expected_params[*]} " =~ " ${param} " ]]; then
      echo "Found unexpected parameter: ${param}"
      ment_param=""
      score=99999
      for exp_param in "${expected_params[@]}"; do
        exp_score=$(levenshtein ${param} ${exp_param})
        if [ "${exp_score}" -lt "${score}" ]; then
          score=${exp_score}
          ment_param=${exp_param}
        fi
      done
      echo "Do you mean ${ment_param}?"
      exit
    fi
  done
}

#!/usr/bin/env bash

read_config_file () {
  config_fname=$1
  includeconfig_regex="includeConfig\s+('|\")(.+)('|\")"
  comment_regex="^(\s*)//.*$"
  emptyline_regex="^(\s|\t|\n)*$"
  keys=()
  key=""
  params=()
  while read line; do
    if [[ "${line}" =~ ${emptyline_regex} ]]; then
      continue
    elif [[ "${line}" =~ ${comment_regex} ]]; then
      continue
    elif [[ "${line}" =~ ${includeconfig_regex} ]]; then
      fname="${BASH_REMATCH[2]}"
      p=$(read_config_file ${fname})
      params+=($p)
    elif [[ "${line}" == *"="* ]]; then
      key_regex="^(.*)=(.*)$"
      if [[ ${line} =~ ${key_regex} ]]; then
        k="${BASH_REMATCH[1]}"
        k=${k/ /}
        a=${keys[*]}
        param=${a// /.}
        params_regex="^params"
        if [[ $param =~ ${params_regex} ]]; then
          if [[ ! ${k} == "dummy" ]]; then
            params+=("${param}.${k}")
          fi
        fi
      fi
    elif [[ "${line}" == *"{"* ]]; then
      key_regex="^(.*)\{$"
      if [[ ${line} =~ ${key_regex} ]]; then
        k="${BASH_REMATCH[1]}"
        k=${k/ /}
        key="${key}${k}"
        keys+=(${key})
        key=""
      fi
    elif [[ "${line}" == *"}"* ]]; then
      unset 'keys[${#keys[@]}-1]'
    fi
  done < ${config_fname}

  echo ${params[@]}
}

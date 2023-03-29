#!/usr/bin/env bash

validate_path() {
  projectDir=$1
  config_file=$2
  while read line; do
    myarray=($( echo ${line} | sed "s/=/ /g;s/\"/ /g;s/'/ /g"))
    for path in ${myarray[@]}; do
      if [[ "${path}" == *"/"* ]]; then
        if [[ ! "${path}" == *"out_dir"* ]]; then
          path=$(eval echo $path)
          if [ ! -f ${path} ]; then
            if [ ! -d ${path} ]; then
              echo "${path} is no such file or directory"
            fi
          fi
        fi
      fi
    done
  done < <(grep -P "/" ${config_file} | grep -vP "^\s*//")
}

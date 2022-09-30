process smurf {
  tag {"SMuRF ${germline_sample_id}"}
  label 'SMuRF'
  label 'SMuRF_3_0_0'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(germline_sample_id), path(germline_vcf), path( germline_tbi), val(bam_sample_ids), val(bam_files), val(bai_files), val(bulk_names) )

  output:
    tuple( val(donor_id), path("${germline_sample_id}.SMuRF.vcf"), path("${germline_sample_id}.SMuRF.filtered.vcf"), path("${germline_sample_id}.SMuRF.vafplot.pdf"), path("${germline_sample_id}_*.SMuRF.filtered.vcf"), emit: somatic_vcfs_dir )

  script:
    b = bam_files ? ' -b ' + bam_files.join(' -b ') : ''
    n = bulk_names ? ' -n ' + bulk_names.join(' -n ') : ''

    """
    host=\$(hostname)
    echo \${host}

    . /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/venv_3.6/bin/activate

    python /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/SMuRF.py \
    -i ${germline_vcf} \
    ${b} \
    ${n} \
    -t ${task.cpus} \
    -c ${params.smurf.config}

    bash /hpc/pmc_vanboxtel/tools/ToolsVanBox/SMuRF-3.0.0/scripts/split_in_single_sample_vcfs.sh ${germline_sample_id}.SMuRF.filtered.vcf

    for BULK in ${n}; do
      if [[ "\${BULK}" != "-n" ]]; then
        rm ${germline_sample_id}_\${BULK}.SMuRF.filtered.vcf*
      fi
    done
    """
}

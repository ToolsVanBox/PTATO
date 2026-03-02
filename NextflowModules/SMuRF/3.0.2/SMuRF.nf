process smurf {
  tag {"SMuRF ${germline_sample_id}"}
  label 'SMuRF'
  label 'SMuRF_3_0_2'
  shell = ['/bin/bash', '-euo', 'pipefail']
  //container = 'docker://vanboxtelbioinformatics/smurf:3.0.2'
//  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'docker://vanboxtelbioinformatics/smurf:3.0.2':
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/smurf@sha256:7db74359db85390f702bf049e937a0fef747f6928222ef4f26f7222504708c5a' }"
    
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/smurf:3.0.4':
        params.artifact_registry_path + '/smurf:3.0.5' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/smurf@sha256:c1a0437a29792cf5ab7fda7c18efefe5bdee9dcc4e5cb77b56950c6ae54ca793' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/smurf@sha256:bbc5bd7d50ed3ebe3ed3a38c8cf89970b244628fe2f0950ea7cbbb27846594ea' }"
  
  input:
    tuple( val(donor_id), val(germline_sample_id), path(germline_vcf), path( germline_tbi), val(bam_sample_ids), path(bam_files), path(bai_files), val(bulk_names) )
    path( config )
  output:
    tuple val(donor_id), path("*SMuRF.vcf"), emit: smurf_vcf
    tuple val(donor_id), path("${germline_sample_id}.SMuRF.filtered.vcf"), emit: smurf_filtered_vcf
    tuple val(donor_id), path("${germline_sample_id}_*.SMuRF.filtered.vcf"), emit: smurf_filtered_single_vcf
//    tuple( val(donor_id), path("${germline_sample_id}.SMuRF.vcf"), path("${germline_sample_id}.SMuRF.filtered.vcf"), path("${germline_sample_id}.SMuRF.vafplot.pdf"), path("${germline_sample_id}_*.SMuRF.filtered.vcf"), emit: somatic_vcfs_dir )

  script:
    b = bam_files ? ' -b ' + bam_files.join(' -b ') : ''
    n = bulk_names ? ' -n ' + bulk_names.join(' -n ') : ''

    """
    host=\$(hostname)
    echo \${host}

    ls -lh ./
    
    export projectDir=${projectDir}

    python /smurf/SMuRF.py \
    -i ${germline_vcf} \
    ${b} \
    ${n} \
    -t ${task.cpus} \
    -c ${config}

    bash /smurf/scripts/split_in_single_sample_vcfs.sh ${germline_sample_id}.SMuRF.filtered.vcf
    
    for BULK in ${n}; do
      if [[ "\${BULK}" != "-n" ]]; then
        rm ${germline_sample_id}_\${BULK}.SMuRF.filtered.vcf*
      fi
    done
    
    touch ${germline_sample_id}.SMuRF.vafplot.pdf
    """
}

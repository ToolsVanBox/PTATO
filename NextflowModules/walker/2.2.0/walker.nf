process walker {
  tag {"walker ${sample_id}"}
  label 'walker'
  label 'walker_2_2_0'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'docker://vanboxtelbioinformatics/walker:2.2.0'
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/walker:2.2.0':
        params.artifact_registry_path + '/walker:2.2.1' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/walker@sha256:002846d2e74e64295da3424f6fe851573e903d887f51c73f06bf54da33c08308' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/walker@sha256:21b524c98bf879d0e7f5ae26c37b70f67f23c8179538994f2661a164e04d1151' }"

  input:
    tuple( val(donor_id), val(germline_sample_id), path(germline_vcf), path( germline_tbi), val(bam_sample_ids), path(bam_files), path(bai_files), val(sample_id), path(somatic_vcf), path(somatic_tbi) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.walker.vcf"), path("${sample_id}.walker.bed"), path("${sample_id}.walker.txt"), emit: walker_out)

  script:
    b = bam_files ? ' -b ' + bam_files.join(' -b ') : ''

    """
    host=\$(hostname)
    echo \${host}
    
    ls -lh ./

    python /walker/walker.py \
    -g ${germline_vcf} \
    -s ${somatic_vcf} \
    -t ${task.cpus} \
    ${b} \
    -o ${sample_id} \
    -f vcf -f bed -f txt

#    touch ${sample_id}.walker.vcf
#    touch ${sample_id}.walker.bed
#    touch ${sample_id}.walker.txt
    
    """
}

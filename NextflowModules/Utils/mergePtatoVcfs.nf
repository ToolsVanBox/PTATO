process mergePtatoVcfs {
  tag {"mergePtatoVcfs ${donor_id}"}
  label 'mergePtatoVcfs'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        params.artifact_registry_path + '/ptato_r:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(sample_id), path(input_vcf), path(input_tbi), val(ptato_snvs_sample_ids), path(ptato_snvs_vcfs), path(ptato_snvs_tbis), val(ptato_indels_sample_ids), path(ptato_indels_vcfs), path(ptato_indels_tbis) )

  output:
    tuple( val(donor_id), val(donor_id), path("${donor_id}.ptato.merged.vcf"), emit: ptato_intersect_vcf)

  script:
    ptato_vcfs = ptato_snvs_vcfs.join(',')+","+ptato_indels_vcfs.join(',')
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/merge_ptato_vcfs.R --args ${input_vcf} ${ptato_vcfs} ${donor_id}.ptato.merged.vcf
    """
}

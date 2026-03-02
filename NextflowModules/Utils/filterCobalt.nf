process FilterCobalt {
  tag {"FilterCobalt ${normal_sample_id} ${tumor_sample_id}"}
  label 'FilterCobalt'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        params.artifact_registry_path + '/ptato_r:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(cobalt_tsv_file) )
    path( centromeres )
    path( cytoband )
    path( pon )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.readcounts.*"), emit: cobalt_filtered_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/filter_cobalt.R --args \
    ${cobalt_tsv_file} \
    ${centromeres} \
    ${cytoband} \
    ${pon} \
    ${params.ref_genome} \
    ${tumor_sample_id} \
    ${params.filtercobalt.optional}
    """
}

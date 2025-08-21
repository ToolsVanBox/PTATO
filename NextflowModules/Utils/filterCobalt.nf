process FilterCobalt {
  tag {"FilterCobalt ${normal_sample_id} ${tumor_sample_id}"}
  label 'FilterCobalt'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_rsha256:01608f437f46324b7ef4c6ae3ce9fe06da8d160ab6d46a44b04342e979f13be6' }"

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(cobalt_tsv_file) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.readcounts.*"), emit: cobalt_filtered_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/filter_cobalt.R --args \
    ${cobalt_tsv_file} \
    ${params.svs.centromeres} \
    ${params.svs.cytoband} \
    ${params.cobalt.pon} \
    ${params.ref_genome} \
    ${tumor_sample_id} \
    ${params.filtercobalt.optional}
    """
}

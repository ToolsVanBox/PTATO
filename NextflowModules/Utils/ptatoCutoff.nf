process ptatoCutoff {
  tag {"ptatoCutoff ${sample_id}"}
  label 'ptatoCutoff'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_rsha256:01608f437f46324b7ef4c6ae3ce9fe06da8d160ab6d46a44b04342e979f13be6' }"

  input:
    tuple( val(donor_id), val(sample_id), path(ptato_vcf), path(ptato_tbi), path(walker_vcf), path(walker_tbi))

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.ptatotable.txt"), path("${sample_id}.ptaprobcutoff.txt"), emit: ptatofilter_out)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/ptatoCutoff.R --args ${ptato_vcf} ${walker_vcf} ${params.ref_genome} ${params.ptatocutoff.optional} ${sample_id}.ptatotable.txt > ${sample_id}.ptaprobcutoff.txt
    """
}

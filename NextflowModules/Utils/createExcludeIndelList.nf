process createIndelExcludeIndelList {
  tag {"createIndelExcludeIndelList"}
  label 'createIndelExcludeIndelList'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_rsha256:01608f437f46324b7ef4c6ae3ce9fe06da8d160ab6d46a44b04342e979f13be6' }"


  input:
    tuple( val(label), path(beds) )

  output:
    tuple( val("excludeindellist"), val("excludeindellist"), path("excludeindellist.bed"), emit: excludeindellist_bed )

  script:
    b1 = beds ? ' ' + beds.join(' ') : ''
    """
    host=\$(hostname)
    echo \${host}

    cat ${b1} | cut -f 1,2,3 | sort | uniq > excludeindellist.bed
    """
}

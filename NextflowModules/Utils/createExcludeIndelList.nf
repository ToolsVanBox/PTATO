process createIndelExcludeIndelList {
  tag {"createIndelExcludeIndelList"}
  label 'createIndelExcludeIndelList'
  shell = ['/bin/bash', '-euo', 'pipefail']

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

process createIndelBlacklist {
  tag {"createIndelBlacklist"}
  label 'createIndelBlacklist'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(label), path(beds) )

  output:
    tuple( val("blacklist"), val("blacklist"), path("blacklist.bed"), emit: blacklist_bed )

  script:
    b1 = beds ? ' ' + beds.join(' ') : ''
    """
    host=\$(hostname)
    echo \${host}

    cat ${b1} | cut -f 1,2,3 | sort | uniq > blacklist.bed
    """
}

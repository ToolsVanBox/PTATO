process AutosomalCallableLoci {
  tag {"AutosomalCallableLoci"}
  label 'AutosomalCallableLoci'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(callable_bed_file), path(callable_txt_file) )
  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.callableloci.autosomal.txt"), emit: autosomal_callable_file )
  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/AutosomalCallableLoci.R --args ${callable_bed_file} ${sample_id}.callableloci.autosomal.txt
    """
}
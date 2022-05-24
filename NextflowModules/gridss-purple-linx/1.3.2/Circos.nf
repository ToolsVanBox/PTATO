
process Circos {

  tag {"Circos ${normal_sample_id} ${tumor_sample_id}"}
  label 'Circos'
  label 'circos'
  container = 'docker://gridss/gridss-purple-linx:1.3.2'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(circos_config_file), path(circos_txt_files) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.*g"), emit: circos_plots )

  script:
    """
    circos \
    -conf ${circos_config_file} \
    -outputfile ${tumor_sample_id} \
    ${params.circos.optional}
    """
}

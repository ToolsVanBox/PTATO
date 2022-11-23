process bgzip {
  tag {"bgzip ${sample_id}"}
  label 'bgzip'
  label 'bgzip_1_15'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/htslib:1.15--h9753748_0'

  input:
    tuple( val(donor_id), val(sample_id), path(vcf) )

  output:
    tuple( val(donor_id), val(sample_id), path("${vcf}.gz"), emit: vcf_gz)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bgzip -c ${vcf} ${params.bgzip.optional} > ${vcf}.gz
    """
}

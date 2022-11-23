process tabix {
  tag {"tabix ${sample_id}"}
  label 'tabix'
  label 'tabix_1_15'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/htslib:1.15--h9753748_0'

  input:
    tuple( val(donor_id), val(sample_id), path(vcf_gz) )

  output:
    tuple( val(donor_id), val(sample_id), path( vcf_gz ), path("${vcf_gz}.tbi"), emit: tbi )

  script:
    """
    host=\$(hostname)
    echo \${host}

    tabix ${params.tabix.optional} ${vcf_gz}
    """
}

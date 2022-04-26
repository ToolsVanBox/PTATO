process SnpSift {
  tag {"SnpSift ${donor_id}"}
  label 'SnpSift'
  label 'SnpSift_4_3_1t__1'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'docker://davelabhub/snpsift:4.3.1t--1'

  input:
    tuple( val(donor_id), val(sample_id), path(vcf_gz), path(tbi), val(normal_sample_ids) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.germline.vcf"), emit: phased_vcf)

  script:
    isHet = normal_sample_ids ? " | isHet( GEN['" + normal_sample_ids.join("'] )  | isHet( GEN['") : ''
    """
    host=\$(hostname)
    echo \${host}

    zcat ${vcf_gz} | \
    SnpSift filter \
    "((exists ID) & ( ID =~ 'rs' )) \
    ${isHet}'] )\
    ${params.snpsift.optional} \
    " > ${sample_id}.germline.vcf
    """
}

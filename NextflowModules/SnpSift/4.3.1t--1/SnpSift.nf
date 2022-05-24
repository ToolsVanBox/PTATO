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
    normals = normal_sample_ids.join(" ")
    """
    host=\$(hostname)
    echo \${host}

    isHET_ARRAY=()
    HEADER=\$(zcat ${vcf_gz} | grep "^#CHROM" | cut -f 10-)

    for NORMAL in ${normals}; do
      NORMAL_GEN_ID=\$(echo \${HEADER/\${NORMAL}//} | cut -d/ -f1 | wc -w | tr -d ' ')
      isHET="isHet(GEN[\${NORMAL_GEN_ID}])"
      isHET_ARRAY+=(\${isHET})
    done
    isHET_STRING=\${isHET_ARRAY[@]}
    isHet=\${isHET_STRING// / | }

    zcat ${vcf_gz} | \
    SnpSift filter \
    "((exists ID) & ( ID =~ 'rs' )) | \
    ( \${isHet} )\
    ${params.snpsift.optional} \
    " > ${sample_id}.germline.vcf

    """

}

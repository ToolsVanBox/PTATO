process shapeit {
  tag {"shapeit ${donor_id} ${chrom}"}
  label 'shapeit'
  label 'shapeit_4_2_2'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/shapeit4:4.2.2--h24bf969_1'

  input:
    tuple( val(donor_id), val(sample_id), path(vcf_gz), path(tbi), val(chrom) )

  output:
    tuple( val(donor_id), val(sample_id), val(chrom), path("${sample_id}.${chrom}.phased.vcf.gz"), emit: phased_vcf)

  script:
    """
    host=\$(hostname)
    echo \${host}

    shapeit4 \
    --input ${vcf_gz} \
    --region ${chrom} \
    --output ${sample_id}.${chrom}.phased.vcf.gz \
    --map ${params.shapeit.maps}/chr${chrom}.b38.gmap.gz \
    --thread ${task.cpus} \
    --sequencing \
    --reference ${params.shapeit.reference}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz \
    ${params.shapeit.optional}
    """
}

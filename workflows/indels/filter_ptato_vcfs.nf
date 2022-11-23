include { ptatoIndelFilter } from '../../NextflowModules/Utils/ptatoFilter.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)

workflow filter_ptato_vcfs {
  take:
    ptato_vcfs

    main:
      ptatoIndelFilter( ptato_vcfs )
      bgzip( ptatoIndelFilter.out )

      tabix( bgzip.out )

      ptato_filtered_vcfs = tabix.out
        .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
          vcf_name = vcf_gz.getName()
          tbi_name = vcf_tbi.getName()
          vcf_gz = vcf_gz.copyTo("${params.out_dir}/indels/${donor_id}/${sample_id}/${vcf_name}")
          vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/indels/${donor_id}/${sample_id}/${tbi_name}")
          [ donor_id, sample_id, vcf_gz, vcf_tbi ]
        }
  emit:
    ptato_filtered_vcfs
}

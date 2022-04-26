include { ptatoFilter } from '../../NextflowModules/Utils/ptatoFilter.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)

workflow filter_ptato_vcfs {
  take:
    ptato_vcfs
    walker_vcfs

  main:
    input_ptato_filter = ptato_vcfs
      .join( walker_vcfs, by: [0,1] )

    ptatoFilter( input_ptato_filter )
    bgzip( ptatoFilter.out
      .map{ donor_id, sample_id, ptato_filtered_vcf, ptato_table ->
        vcf_name = ptato_filtered_vcf.getName()
        table_name = ptato_table.getName()
        ptato_table.copyTo("${params.out_dir}/indels/${donor_id}/${sample_id}/${table_name}")
        [ donor_id, sample_id, ptato_filtered_vcf ]
      }
    )
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

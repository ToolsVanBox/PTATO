include { test_indel_rf } from '../../NextflowModules/Utils/randomForest.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)

include { extractPtatoVcfFromDir } from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow get_ptato_indel_vcfs {
  take:
    somatic_vcfs
    rf_tables
  main:
    if ( params.optional.indels.ptato_vcfs_dir ) {
      raw_ptato_vcfs = extractPtatoVcfFromDir( params.optional.indels.ptato_vcfs_dir )
      get_gzipped_vcfs( raw_ptato_vcfs )
      ptato_vcfs = get_gzipped_vcfs.out
    } else {
      input_test_indel_rf = somatic_vcfs.combine( rf_tables, by: [0,1] )

      test_indel_rf( input_test_indel_rf )

      bgzip( test_indel_rf.out )
      tabix( bgzip.out )
      ptato_vcfs = tabix.out
        .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
          vcf_name = vcf_gz.getName()
          tbi_name = vcf_tbi.getName()
          vcf_gz = vcf_gz.copyTo("${params.out_dir}/indels/${donor_id}/${vcf_name}")
          vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/indels/${donor_id}/${tbi_name}")
          [ donor_id, sample_id, vcf_gz, vcf_tbi ]
        }
    }
  emit:
    ptato_vcfs
}

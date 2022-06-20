include { intersectPON } from '../../NextflowModules/bedtools/2.30.0/intersect.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)

include { extractPonFilteredVcfFromDir } from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow get_pon_filtered_vcfs {
  take:
    somatic_vcfs
  main:
    if ( params.optional.indels.pon_filtered_vcfs_dir ) {
      raw_pon_filtered_vcfs = extractPonFilteredVcfFromDir( params.optional.indels.pon_filtered_vcfs_dir )
      get_gzipped_vcfs( raw_pon_filtered_vcfs )
      pon_filtered_vcfs = get_gzipped_vcfs.out
    } else {
      input_files = somatic_vcfs.combine( Channel.from( params.indels.pon ) )

      intersectPON( input_files )
      bgzip( intersectPON.out )
      tabix( bgzip.out )
      pon_filtered_vcfs = tabix.out
        .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
          vcf_name = vcf_gz.getName()
          tbi_name = vcf_tbi.getName()
          vcf_gz = vcf_gz.copyTo("${params.out_dir}/intermediate/indels/pon/${donor_id}/${vcf_name}")
          vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/intermediate/indels/pon/${donor_id}/${tbi_name}")
          [ donor_id, sample_id, vcf_gz, vcf_tbi ]
        }
    }
  emit:
    pon_filtered_vcfs
}

include { intersectExcludeIndelList } from '../../NextflowModules/bedtools/2.30.0/intersect.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)

include { extractExcludeIndelListFilteredVcfFromDir } from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow get_excludeindellist_filtered_vcfs {
  take:
    somatic_vcfs
  main:
    if ( params.optional.indels.excludeindellist_filtered_vcfs_dir ) {
      raw_excludeindellist_filtered_vcfs = extractexcludeindellistFilteredVcfFromDir( params.optional.indels.excludeindellist_filtered_vcfs_dir )
      get_gzipped_vcfs( raw_excludeindellist_filtered_vcfs )
      excludeindellist_filtered_vcfs = get_gzipped_vcfs.out
    } else {
      input_files = somatic_vcfs.combine( Channel.from( params.indels.excludeindellist ) )

      intersectExcludeIndelList( input_files )
      bgzip( intersectExcludeIndelList.out )
      tabix( bgzip.out )
      excludeindellist_filtered_vcfs = tabix.out
        .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
          vcf_name = vcf_gz.getName()
          tbi_name = vcf_tbi.getName()
          vcf_gz = vcf_gz.copyTo("${params.out_dir}/intermediate/indels/excludeindellist/${donor_id}/${vcf_name}")
          vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/intermediate/indels/excludeindellist/${donor_id}/${tbi_name}")
          [ donor_id, sample_id, vcf_gz, vcf_tbi ]
        }
    }
  emit:
    excludeindellist_filtered_vcfs
}

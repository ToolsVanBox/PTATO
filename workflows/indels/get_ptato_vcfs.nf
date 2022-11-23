include { excludeIndels } from '../../NextflowModules/Utils/excludeIndels.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)

include { extractPtatoVcfFromDir } from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow get_ptato_vcfs {
  take:
    somatic_vcfs
    context_beds

  main:
    if ( params.optional.indels.ptato_vcfs_dir ) {
      raw_ptato_vcfs = extractPtatoVcfFromDir( params.optional.indels.ptato_vcfs_dir )
      get_gzipped_vcfs( raw_ptato_vcfs )
      ptato_vcfs = get_gzipped_vcfs.out
    } else {
      // somatic_vcfs.view()
      // context_beds.view()
      // params.indels.excludeindellist.view()

      input_files = somatic_vcfs
        .combine( context_beds, by: [0,1] )
        .combine( Channel.from( params.indels.excludeindellist ) )


      excludeIndels( input_files )

      bgzip( excludeIndels.out )
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

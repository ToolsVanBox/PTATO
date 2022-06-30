include { mergePtatoVcfs } from '../../NextflowModules/Utils/mergePtatoVcfs.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)

workflow merge_ptato_vcfs {
  take:
    ptato_intersect_vcfs
    snvs_ptato_vcfs
    indels_ptato_vcfs

  main:
  input_merge_vcfs = ptato_intersect_vcfs
    .combine(
      snvs_ptato_vcfs
        .join( indels_ptato_vcfs )
        .groupTuple( by: [0] )
    , by: [0] )

  mergePtatoVcfs( input_merge_vcfs )
  bgzip( mergePtatoVcfs.out )
  tabix( bgzip.out )
  ptato_merged_vcfs = tabix.out
    .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
      vcf_name = vcf_gz.getName()
      tbi_name = vcf_tbi.getName()
      vcf_gz = vcf_gz.copyTo("${params.out_dir}/ptato_vcfs/${donor_id}/${vcf_name}")
      vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/ptato_vcfs/${donor_id}/${tbi_name}")
      [ donor_id, sample_id, vcf_gz, vcf_tbi ]
    }

  emit:
    ptato_merged_vcfs
}

include { intersectPTATO } from '../../NextflowModules/bedtools/2.30.0/intersect.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)

workflow intersect_ptato_vcfs {
  take:
    input_vcfs
    snvs_ptato_vcfs
    indels_ptato_vcfs

  main:
    input_intersect_vcfs = input_vcfs
      .combine(
        snvs_ptato_vcfs
          .join( indels_ptato_vcfs )
          .groupTuple( by: [0] )
      , by: [0] )

    intersectPTATO( input_intersect_vcfs )
    bgzip( intersectPTATO.out )
    tabix( bgzip.out )
    ptato_intersect_vcfs = tabix.out

  emit:
    ptato_intersect_vcfs
}

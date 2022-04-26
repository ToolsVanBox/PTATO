#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { get_germline_vcfs } from './workflows/germline.nf' params(params)
include { short_variants_train } from './workflows/short_variants.nf' params(params)
include {
  extractInputVcfGzFromDir;
  extractBamsFromDir
} from './NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow {
  main:
    input_vcfs = extractInputVcfGzFromDir( params.input_vcfs_dir )
    input_bams = extractBamsFromDir( params.bams_dir ).groupTuple( by: [0] )

    GetSampleName( input_bams.transpose() )
    bams = GetSampleName.out
      .groupTuple( by: [0] )

    if ( params.run.snvs || params.run.indels || params.run.svs ) {
      get_germline_vcfs( input_vcfs )
      germline_vcfs = get_germline_vcfs.out

      if ( params.run.snvs || params.run.indels ) {
        short_variants_train( input_vcfs, bams, germline_vcfs )
      }
    }

    // if ( params.QC ) {
    //   qc()
    // }
    // if ( params.snvs ) {
      // snvs( germline_vcfs )
    // }
    // if ( params.svs ) {
    //   svs()
    // }
}

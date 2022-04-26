#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { GetSampleName } from './NextflowModules/GATK/4.1.3.0/GetSampleName.nf' params(params)
include { get_germline_vcfs } from './workflows/germline.nf' params(params)
include { short_variants } from './workflows/short_variants.nf' params(params)
include { svs } from './workflows/svs.nf' params(params)
include { qc } from './workflows/qc.nf' params(params)

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
        short_variants( input_vcfs, bams, germline_vcfs )
      }

      if ( params.run.svs ) {
        svs( bams )
      }
    }
    if ( params.run.QC ) {
      qc( bams )
    }
}

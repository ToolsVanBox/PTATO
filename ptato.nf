#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include{ check_params } from './NextflowModules/Utils/checkInputParams.nf'

include { get_gzipped_vcfs } from './workflows/get_gzipped_vcfs.nf' params(params)
include { get_indexed_bams } from './workflows/get_indexed_bams.nf' params(params)

include { get_germline_vcfs } from './workflows/germline.nf' params(params)
include { short_variants } from './workflows/short_variants.nf' params(params)
include { svs } from './workflows/svs.nf' params(params)
include { cnvs } from './workflows/cnvs.nf' params(params)
include { qc } from './workflows/qc.nf' params(params)

include {
  extractInputVcfFromDir;
  extractBamsFromDir
} from './NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow {
  main:
    run_donor_ids = Channel.from( params.bulk_names )
      .map{
        donor_id, bulk_name ->
        [donor_id]
      }
      .unique()

    input_raw_vcfs = run_donor_ids.combine( extractInputVcfFromDir( params.input_vcfs_dir ), by: [0] )
    input_raw_bams = run_donor_ids.combine( extractBamsFromDir( params.bams_dir ), by: [0] )

    get_gzipped_vcfs( input_raw_vcfs )
    input_vcfs = get_gzipped_vcfs.out

    get_indexed_bams( input_raw_bams )
    input_bams = get_indexed_bams.out.groupTuple( by: [0] )

    if ( params.run.snvs || params.run.indels || params.run.svs || params.run.cnvs ) {

      get_germline_vcfs( input_vcfs )
      germline_vcfs = get_germline_vcfs.out

      if ( params.run.snvs || params.run.indels ) {
        short_variants( input_vcfs, input_bams, germline_vcfs )
      }

      if ( params.run.svs || params.run.cnvs ) {
        cnvs( input_bams, germline_vcfs )
        filtered_cnv_files = cnvs.out


        if ( params.run.svs ) {
          svs( input_bams, germline_vcfs, filtered_cnv_files )
        }
      }
    }
    if ( params.run.QC ) {
      qc( input_bams )
    }
}

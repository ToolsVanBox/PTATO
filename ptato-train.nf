#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { get_gzipped_vcfs } from './workflows/get_gzipped_vcfs.nf' params(params)

include { get_germline_vcfs } from './workflows/germline.nf' params(params)
include { short_variants_train } from './workflows/short_variants.nf' params(params)
include {
  extractInputVcfGzFromDir
} from './NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow {
  main:
    run_donor_ids = Channel.from( params.bulk_names )
      .map{
        donor_id, bulk_name ->
        [donor_id]
      }

    input_raw_vcfs = run_donor_ids.combine( extractInputVcfGzFromDir( params.input_vcfs_dir ), by: [0] )
    get_gzipped_vcfs( input_raw_vcfs )
    input_vcfs = get_gzipped_vcfs.out

    if ( params.run.snvs || params.run.indels || params.run.svs ) {
      get_germline_vcfs( input_vcfs )
      germline_vcfs = get_germline_vcfs.out

      if ( params.run.snvs || params.run.indels ) {
        short_variants_train( input_vcfs, germline_vcfs )
      }
    }
}

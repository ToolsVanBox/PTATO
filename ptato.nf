#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include{ check_params } from './NextflowModules/Utils/checkInputParams.nf'

include { get_gzipped_vcfs } from './workflows/get_gzipped_vcfs.nf' params(params)
include { get_gzipped_vcfs as get_gzipped_vcfs2 } from './workflows/get_gzipped_vcfs.nf' params(params)
include { get_indexed_bams } from './workflows/get_indexed_bams.nf' params(params)

include { get_germline_vcfs } from './workflows/germline.nf' params(params)
include { short_variants } from './workflows/short_variants.nf' params(params)

include { svs } from './workflows/svs.nf' params(params)
include { cnvs } from './workflows/cnvs.nf' params(params)
include { RunCallableLoci } from './workflows/QC/RunCallableLoci.nf' params(params)

include { 
  qc; 
  post_ptato_qc 
} from './workflows/qc.nf' params(params)

include {
  extractInputVcfFromDir;
  extractBamsFromDir;
  extractWalkerVcfFromDir;
  extractCombinedPtatoVcfFromDir; 
  extractPtatoTableFromDir
  extractInputVcfFromCloudDir;
  extractBamsFromCloudDir;
  extractBaisFromCloudBatchDir;
} from './NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow {
  main:
    run_donor_ids = Channel.from( params.bulk_names )
      .map{
        donor_id, bulk_name ->
        [donor_id]
      }
      .unique()

    def donor_id = params.bulk_names[0][0]
    
//    input_raw_vcfs = run_donor_ids.combine( extractInputVcfFromDir( params.input_vcfs_dir ), by: [0] )
//    input_raw_bams = run_donor_ids.combine( extractBamsFromDir( params.bams_dir ), by: [0] )
    
    Channel.fromPath( params.input_samplesheet ).
        splitCsv( header:true )
        .branch { row ->
            ch_bam: row.file_type == "bam"
                    return tuple( row.donor_id, row.sample_id, row.file, row.file_index )
            ch_vcf: row.file_type == "vcf"
                    return tuple( row.donor_id, row.sample_id, row.file, row.file_index )
           
        }
        .set { inputs }

    
    input_raw_vcfs = run_donor_ids.combine( inputs.ch_vcf, by: [0] )
    input_raw_bams = run_donor_ids.combine( inputs.ch_bam, by: [0] )

    input_raw_vcfs.view()
    
    input_raw_bams.view()

//    input_raw_bams = run_donor_ids.combine( extractBamsFromCloudDir( params.bams_dir, donor_id ), by: [0] )
//    input_bais = run_donor_ids.combine( extractBaisFromCloudBatchDir( "${params.input_vcfs_dir}/../../../bams/", donor_id ), by: [0] )
    
//    input_raw_bams = input_bais.combine(input_raw_bams, by: [0,1]).view()
    
    get_indexed_bams( input_raw_bams )
    input_bams = get_indexed_bams.out.groupTuple( by: [0] )

    // Define variables
    def fasta = file( params.genome_fasta, checkIfExists: true )
    def fai = file( params.genome_fai, checkIfExists: true )
    def dict = file( params.genome_dict, checkIfExists: true )
    
    ch_fasta = Channel.value( fasta )
      .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }
      
    ch_fai = Channel.value( fai )
      .map{ genome_fai -> [ [ id:'fai' ], genome_fai ] }
    
    ch_dict = Channel.value( dict )
      .map( genome_dict -> [ [ id:'dict'], genome_dict ] )
      
    if ( params.run.QC || params.run.postqc) {
      RunCallableLoci( input_bams, ch_fasta, ch_fai, ch_dict )
      if ( params.run.QC) {
        qc( input_bams, ch_fasta, ch_fai, ch_dict )
      }
    }

    if ( params.run.snvs || params.run.indels || params.run.svs || params.run.cnvs || params.run.postqc) {

      get_gzipped_vcfs( input_raw_vcfs )
      input_vcfs = get_gzipped_vcfs.out

      get_germline_vcfs( input_vcfs )
      germline_vcfs = get_germline_vcfs.out

      if( params.run.postqc ) {
        if ( params.optional.postqc.ptato_vcfs_dir && params.optional.walker_vcfs_dir ) {
          ptato_combined_vcfs = extractCombinedPtatoVcfFromDir( params.optional.postqc.ptato_vcfs_dir )
          raw_walker_vcfs = extractWalkerVcfFromDir( params.optional.walker_vcfs_dir )
          get_gzipped_vcfs2( raw_walker_vcfs )
          walker_vcfs = get_gzipped_vcfs2.out
          ptato_tables = extractPtatoTableFromDir(params.optional.postqc.ptato_vcfs_dir ) 

          postqc_combined_input = ptato_combined_vcfs.combine(
            ptato_tables, by: [0,1] ).combine(
            walker_vcfs, by: [0,1] )

          if ( params.run.snvs || params.run.indels ) {
            short_variants( input_vcfs, input_bams, germline_vcfs )
          }
        } else {
          // You need to run short_variant before you can run postQC 
          short_variants( input_vcfs, input_bams, germline_vcfs )
          postqc_combined_input = short_variants.out
        }
        post_ptato_qc( postqc_combined_input, RunCallableLoci.out )
      } else {
        if ( params.run.snvs || params.run.indels ) {
         short_variants( input_vcfs, input_bams, germline_vcfs )
        }
      }

      if ( params.run.svs || params.run.cnvs ) {
        cnvs( input_bams, germline_vcfs, ch_fasta, ch_fai, ch_dict )
        filtered_cnv_files = cnvs.out

        if ( params.run.svs ) {
          svs( input_bams, germline_vcfs, filtered_cnv_files, ch_fasta, ch_fai, ch_dict )
        }
      }
    }
}

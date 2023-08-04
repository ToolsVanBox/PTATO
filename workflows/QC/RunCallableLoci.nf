#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { extractCallableLociBedFromDir; 
  extractAutosomalCallableLociFromDir } from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)
include { CallableLoci } from '../../NextflowModules/GATK/3.8.1/CallableLoci.nf' params(params)
include { AutosomalCallableLoci } from '../../NextflowModules/Utils/AutosomalCallable.nf' params(params)


workflow RunCallableLoci {
  take:
    bams
  main:
    if ( params.optional.callableloci_dir ) {
      callableloci_files = extractCallableLociBedFromDir( params.optional.callableloci_dir )
    } else {
      CallableLoci( bams.transpose() )
      callableloci_files = CallableLoci.out
        .map{ donor_id, sample_id, callableloci_bed, callableloci_txt ->
          bed_filename = callableloci_bed.getName()
          txt_filename = callableloci_txt.getName()
          callableloci_bed = callableloci_bed.copyTo("${params.out_dir}/QC/CallableLoci/${donor_id}/${bed_filename}")
          callableloci_txt = callableloci_txt.copyTo("${params.out_dir}/QC/CallableLoci/${donor_id}/${txt_filename}")
          [ donor_id, sample_id, callableloci_bed, callableloci_txt ]
        }
    }

    // Create the automal call for R script 
    if ( params.optional.autosomal_callable_dir ) {
      autosomal_callable_files = extractAutosomalCallableLociFromDir( params.optional.autosomal_callable_dir )
    } else {
      AutosomalCallableLoci( callableloci_files )
      autosomal_callable_files = AutosomalCallableLoci.out
        .map{ donor_id, sample_id, autosomal_callable_file -> 
          auto_filename = autosomal_callable_file.getName()
          autosomal_callable_file = autosomal_callable_file.copyTo("${params.out_dir}/QC/AutosomalCallableLoci/${donor_id}/${auto_filename}")
          [donor_id, sample_id, autosomal_callable_file]
        }
    }
    // Combine and emit output 
    callable_combined_files = callableloci_files.combine(autosomal_callable_files, by: [0,1])
    emit: 
      callable_combined_files
}

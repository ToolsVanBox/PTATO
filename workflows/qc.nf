#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
  extractWGSMetricsFromDir;
  extractAlignmentSummaryMetricsFromDir;
  extractCallableLociBedFromDir
} from '../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { CollectWGSMetrics } from '../NextflowModules/GATK/4.2.6.1/CollectWGSMetrics.nf' params(params)
include { CollectAlignmentSummaryMetrics } from '../NextflowModules/GATK/4.2.6.1/CollectAlignmentSummaryMetrics.nf' params(params)
include { CallableLoci } from '../NextflowModules/GATK/3.8.1/CallableLoci.nf' params(params)
include { QCreport } from '../NextflowModules/Utils/QCreport.nf' params(params)
include { postQCreport } from '../NextflowModules/Utils/postQCreport.nf' params(params)

workflow qc {
  take:
    bams
  main :
    if ( params.optional.qc.alignment_summary_metrics_dir ) {
      alignment_summary_metrics_files = extractAlignmentSummaryMetricsFromDir( params.optional.qc.alignment_summary_metrics_dir )
    } else {
      CollectAlignmentSummaryMetrics( bams.transpose() )
      alignment_summary_metrics_files = CollectAlignmentSummaryMetrics.out
        .map{ donor_id, sample_id, alignment_summary_metrics_file ->
          filename = alignment_summary_metrics_file.getName()
          alignment_summary_metrics_file = alignment_summary_metrics_file.copyTo("${params.out_dir}/QC/alignment_summary_metrics/${donor_id}/${filename}")
          [ donor_id, sample_id, alignment_summary_metrics_file ]
        }
    }

    if ( params.optional.qc.wgs_metrics_dir ) {
      wgs_metrics_files = extractWGSMetricsFromDir( params.optional.qc.wgs_metrics_dir )
    } else {
      CollectWGSMetrics( bams.transpose() )
      wgs_metrics_files = CollectWGSMetrics.out
        .map{ donor_id, sample_id, wgs_metrics_file ->
          filename = wgs_metrics_file.getName()
          wgs_metrics_file = wgs_metrics_file.copyTo("${params.out_dir}/QC/wgs_metrics/${donor_id}/${filename}")
          [ donor_id, sample_id, wgs_metrics_file ]
        }
    }

    input_qc_report = alignment_summary_metrics_files
        .combine( wgs_metrics_files, by: [0,1] )
        .groupTuple( by: [0] )

    QCreport( input_qc_report )

    qc_reports = QCreport.out
      .map{ donor_id, qc_report_pdf, qc_report_txt ->
        pdf_name = qc_report_pdf.getName()
        txt_name = qc_report_txt.getName()
        qc_report_pdf = qc_report_pdf.copyTo("${params.out_dir}/QC/reports/${donor_id}/${pdf_name}")
        qc_report_txt = qc_report_txt.copyTo("${params.out_dir}/QC/reports/${donor_id}/${txt_name}")
        [ donor_id, qc_report_pdf, qc_report_txt ]
      }
}

workflow post_ptato_qc {
  take:
    postqc_combined_input
    callableloci_files
  main:
    // Rename the sample_id in order to combine them and select necessary values 
    postqc_combined_input = postqc_combined_input.map{ donor_id, sample_id, ptato_vcf, ptato_tbi, ptato_filt_vcf, ptato_filt_tbi, walker_vcf, walker_tbi, ptato_table -> 
      // This might be a tricky part of the postQC
      sample_id = sample_id.replaceAll(/.*_/, '')
      [ donor_id, sample_id, ptato_vcf, ptato_tbi, ptato_filt_vcf, ptato_filt_tbi, walker_vcf, walker_tbi, ptato_table ]
    }
    callableloci_files = callableloci_files.map{ donor_id, sample_id, callableloci_bed, callableloci_txt, autosomal_callable_files -> 
      [ donor_id, sample_id, autosomal_callable_files ]
    }
    
    // Combine the input 
    input_PostPTATO_qc = postqc_combined_input
        .combine( callableloci_files, by: [0,1] )
        .groupTuple( by: [0] )
    
    postQCreport( input_PostPTATO_qc )
    postQC_report = postQCreport.out
      .map{ donor_id, sample_id, postqc_report_pdf, postqc_report_txt -> 
        pdf_name = postqc_report_pdf.getName()
        txt_name = postqc_report_txt.getName()
        postqc_report_pdf = postqc_report_pdf.copyTo("${params.out_dir}/QC/reports/${donor_id}/${pdf_name}")
        postqc_report_txt = postqc_report_txt.copyTo("${params.out_dir}/QC/reports/${donor_id}/${txt_name}")
        [ donor_id, sample_id, postqc_report_pdf, postqc_report_txt ]
    }
}


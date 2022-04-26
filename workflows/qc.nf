#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
  extractWGSMetricsFromDir;
  extractAlignmentSummaryMetricsFromDir
} from '../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { CollectWGSMetrics } from '../NextflowModules/GATK/4.1.3.0/CollectWGSMetrics.nf' params(params)
include { CollectAlignmentSummaryMetrics } from '../NextflowModules/GATK/4.1.3.0/CollectAlignmentSummaryMetrics.nf' params(params)
include { QCreport } from '../NextflowModules/Utils/QCreport.nf' params(params)

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
        .collect()


    QCreport( input_qc_report )

    qc_reports = QCreport.out
      .map{ donor_id, qc_report_pdf ->
        pdf_name = qc_report_pdf.getName()
        qc_report_pdf = qc_report_pdf.copyTo("${params.out_dir}/QC/${pdf_name}")
        [ donor_id, qc_report_pdf ]
      }
}

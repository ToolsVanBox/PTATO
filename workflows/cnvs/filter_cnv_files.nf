
include {
  extractCobaltFilteredReadCounts;
  extractBafFilteredFiles;
  extractCobaltFilteredReadCountsSegments;
  extractBafSegmentsFiles;
  extractBafBinnedFiles;
} from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { FilterCobalt } from '../../NextflowModules/Utils/filterCobalt.nf' params(params)
include { FilterBAF } from '../../NextflowModules/Utils/filterBAF.nf' params(params)

workflow filter_cnv_files {
  take:
    cobalt_ratio_tsv_files
    germline_vcf_files
    normal_bams
    tumor_bams
  main:
    if ( params.optional.cnvs.cobalt_filtered_readcounts_dir ) {
      cobalt_filtered_readcounts_files = extractCobaltFilteredReadCounts( params.optional.cnvs.cobalt_filtered_readcounts_dir )
      cobalt_filtered_readcounts_segment_files = extractCobaltFilteredReadCountsSegments( params.optional.cnvs.cobalt_filtered_readcounts_dir )
      cobalt_filtered_readcounts_files = cobalt_filtered_readcounts_files.concat( cobalt_filtered_readcounts_segment_files )
    } else {
      FilterCobalt( cobalt_ratio_tsv_files )
      cobalt_filtered_readcounts_files = FilterCobalt.out
        .transpose()
        .map{
          donor_id, normal_sample_id, tumor_sample_id, cobalt_filtered_file ->
          file_name = cobalt_filtered_file.getName()
          cobalt_filtered_file.copyTo("${params.out_dir}/intermediate/cnvs/readCounts/${donor_id}/${normal_sample_id}/${file_name}")
          [ donor_id, normal_sample_id, tumor_sample_id, cobalt_filtered_file ]
        }
    }

    if ( params.optional.cnvs.baf_filtered_files_dir ) {
      baf_filtered_files = extractBafFilteredFiles( params.optional.cnvs.baf_filtered_files_dir )
      baf_filtered_segments_files = extractBafSegmentsFiles( params.optional.cnvs.baf_filtered_files_dir )
      // baf_filtered_segments_files.view()
      baf_binned_files = extractBafBinnedFiles( params.optional.cnvs.baf_filtered_files_dir )
      baf_filtered_files = baf_filtered_files.concat( baf_filtered_segments_files ).concat( baf_binned_files)
      baf_filtered_files.view()
    } else {
      input_baf_filter_files = normal_bams
        .combine( tumor_bams, by: [0] )
        .combine( germline_vcf_files, by: [0] )
        .map{
          donor_id, normal_sample_id, normal_bam, normal_bai, tumor_sample_id, tumor_bam, tumor_bai, sample_id, germline_vcf, germline_tbi ->
          [ donor_id, normal_sample_id, tumor_sample_id, germline_vcf, germline_tbi ]
        }

      FilterBAF( input_baf_filter_files )
      baf_filtered_files = FilterBAF.out
        .transpose()
        .map{
          donor_id, normal_sample_id, tumor_sample_id, baf_filtered_file ->
          file_name = baf_filtered_file.getName()
          baf_filtered_file.copyTo("${params.out_dir}/intermediate/cnvs/BAF/${donor_id}/${normal_sample_id}/${file_name}")
          [ donor_id, normal_sample_id, tumor_sample_id, baf_filtered_file ]
        }
    }
    filtered_cnv_files = cobalt_filtered_readcounts_files.concat( baf_filtered_files )
  emit:
    filtered_cnv_files
}

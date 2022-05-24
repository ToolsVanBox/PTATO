
include {
  extractCobaltFilteredReadCounts;
  extractBafFilteredFiles;
  extractCobaltFilteredReadCountsSegments;
  extractBafFilteredSegments;
  extractBafBinnedFiles;
  extractGripssFilteredFiles;
  extractIntegratedSvFiles;
} from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { FilterCobalt } from '../../NextflowModules/Utils/filterCobalt.nf' params(params)
include { FilterBAF } from '../../NextflowModules/Utils/filterBAF.nf' params(params)
include { FilterGripss } from '../../NextflowModules/Utils/filterGripss.nf' params(params)
include { IntegrateSvFiles } from '../../NextflowModules/Utils/integrateSvFiles.nf' params(params)
include { CreateSvPlots } from '../../NextflowModules/Utils/createSvPlots.nf' params(params)
include { CreateCircosConfig } from '../../NextflowModules/Utils/createCircosConfig.nf' params(params)
include { Circos } from '../../NextflowModules/gridss-purple-linx/1.3.2/Circos.nf' params(params)

workflow filter_sv_files {
  take:
    cobalt_ratio_tsv_files
    germline_vcf_files
    normal_bams
    tumor_bams
    gripss_somatic_filtered_vcfs
  main:
    if ( params.optional.svs.cobalt_filtered_readcounts_dir ) {
      cobalt_filtered_readcounts_files = extractCobaltFilteredReadCounts( params.optional.svs.cobalt_filtered_readcounts_dir )
      cobalt_filtered_readcounts_1kb_files = cobalt_filtered_readcounts_files.filter{ it[3] =~ /.readcounts.filtered.1kb.txt/ }
      cobalt_filtered_readcounts_100kb_files = cobalt_filtered_readcounts_files.filter{ it[3] =~ /.readcounts.filtered.100kb.txt/ }
      cobalt_filtered_readcounts_1mb_files = cobalt_filtered_readcounts_files.filter{ it[3] =~ /.readcounts.filtered.1mb.txt/ }

      cobalt_filtered_readcounts_segments = extractCobaltFilteredReadCountsSegments( params.optional.svs.cobalt_filtered_readcounts_dir )
      cobalt_filtered_readcounts_segments_bedpe_files = cobalt_filtered_readcounts_segments.filter{ it[3] =~ /.readcounts.segments.bedpe/ }
      cobalt_filtered_readcounts_segments_txt_files = cobalt_filtered_readcounts_segments.filter{ it[3] =~ /.readcounts.segments.txt/ }
    } else {
      FilterCobalt( cobalt_ratio_tsv_files )
      filter_cobalt_files = FilterCobalt.out
        .transpose()
        .map{
          donor_id, normal_sample_id, tumor_sample_id, cobalt_filtered_file ->
          file_name = cobalt_filtered_file.getName()
          cobalt_filtered_file.copyTo("${params.out_dir}/intermediate/svs/readCounts/${donor_id}/${normal_sample_id}/${file_name}")
          [ donor_id, normal_sample_id, tumor_sample_id, cobalt_filtered_file ]
        }
      cobalt_filtered_readcounts_segments_bedpe_files = filter_cobalt_files.filter{ it[3] =~ /.readcounts.segments.bedpe/ }
      cobalt_filtered_readcounts_segments_txt_files = filter_cobalt_files.filter{ it[3] =~ /.readcounts.segments.txt/ }
      cobalt_filtered_readcounts_1kb_files = filter_cobalt_files.filter{ it[3] =~ /.readcounts.filtered.1kb.txt/ }
      cobalt_filtered_readcounts_100kb_files = filter_cobalt_files.filter{ it[3] =~ /.readcounts.filtered.100kb.txt/ }
      cobalt_filtered_readcounts_1mb_files = filter_cobalt_files.filter{ it[3] =~ /.readcounts.filtered.1mb.txt/ }
    }

    if ( params.optional.svs.baf_filtered_files_dir ) {
      baf_filtered_files = extractBafFilteredFiles( params.optional.svs.baf_filtered_files_dir )
      baf_filtered_segments_files = extractBafFilteredSegments( params.optional.svs.baf_filtered_files_dir )
      baf_filtered_segments_txt_files = baf_filtered_segments_files.filter{ it[3] =~ /.baf.segments.txt/ }
      baf_filtered_segments_bedpe_files = baf_filtered_segments_files.filter{ it[3] =~ /.baf.segments.bedpe/ }
      baf_binned_files = extractBafBinnedFiles( params.optional.svs.baf_filtered_files_dir )
      baf_binned_100kb_files = baf_binned_files.filter{ it[3] =~ /.baf.binned.100kb.txt/ }
      baf_binned_1mb_files = baf_binned_files.filter{ it[3] =~ /.baf.binned.1mb.txt/ }
    } else {
      input_baf_filter_files = normal_bams
        .combine( tumor_bams, by: [0] )
        .combine( germline_vcf_files, by: [0] )
        .map{
          donor_id, normal_sample_id, normal_bam, normal_bai, tumor_sample_id, tumor_bam, tumor_bai, sample_id, germline_vcf, germline_tbi ->
          [ donor_id, normal_sample_id, tumor_sample_id, germline_vcf, germline_tbi ]
        }

      FilterBAF( input_baf_filter_files )
      filter_baf_files = FilterBAF.out
        .transpose()
        .map{
          donor_id, normal_sample_id, tumor_sample_id, baf_filtered_file ->
          file_name = baf_filtered_file.getName()
          baf_filtered_file.copyTo("${params.out_dir}/intermediate/svs/BAF/${donor_id}/${normal_sample_id}/${file_name}")
          [ donor_id, normal_sample_id, tumor_sample_id, baf_filtered_file ]
        }

      baf_filtered_files = filter_baf_files.filter{ it[3] =~ /.baf.filtered.txt/ }
      baf_filtered_segments_bedpe_files = filter_baf_files.filter{ it[3] =~ /.baf.segments.bedpe/ }
      baf_filtered_segments_txt_files = filter_baf_files.filter{ it[3] =~ /.baf.segments.txt/ }
      baf_binned_100kb_files = filter_baf_files.filter{ it[3] =~ /.baf.binned.100kb.txt/ }
      baf_binned_1mb_files = filter_baf_files.filter{ it[3] =~ /.baf.binned.1mb.txt/ }

    }

    if ( params.optional.svs.gripss_filtered_files_dir ) {
      gripss_filtered_files = extractGripssFilteredFiles( params.optional.svs.gripss_filtered_files_dir )
    } else {
      input_gripss_filter_files = gripss_somatic_filtered_vcfs
        .combine( baf_filtered_files, by: [0,1,2] )
        .combine( cobalt_filtered_readcounts_1kb_files, by: [0,1,2] )

      FilterGripss( input_gripss_filter_files )
      filter_gripss_files = FilterGripss.out
        .transpose()
        .map{
          donor_id, normal_sample_id, tumor_sample_id, gripss_filtered_file ->
          file_name = gripss_filtered_file.getName()
          gripss_filtered_file.copyTo("${params.out_dir}/intermediate/svs/Breakends/${donor_id}/${normal_sample_id}/${file_name}")
          [ donor_id, normal_sample_id, tumor_sample_id, gripss_filtered_file ]
        }
      gripss_filtered_files = filter_gripss_files.filter{ it[3] =~ /.svs.postfilter.vcf/ }

    }
    if ( params.optional.svs.integrated_sv_files_dir ) {
      integrated_sv_files = extractIntegratedSvFiles( params.optional.svs.integrated_sv_files_dir )
    } else {
      input_integrate_sv_files = baf_filtered_files
        .combine( cobalt_filtered_readcounts_1kb_files, by: [0,1,2] )
        .combine( baf_filtered_segments_bedpe_files, by: [0,1,2] )
        .combine( cobalt_filtered_readcounts_segments_bedpe_files, by: [0,1,2] )
        .combine( gripss_filtered_files, by: [0,1,2] )

      IntegrateSvFiles( input_integrate_sv_files )
      integrated_sv_files = IntegrateSvFiles.out
        .transpose()
        .map{
          donor_id, normal_sample_id, tumor_sample_id, gripss_filtered_file ->
          file_name = gripss_filtered_file.getName()
          gripss_filtered_file.copyTo("${params.out_dir}/intermediate/svs/Integration/${donor_id}/${normal_sample_id}/${file_name}")
          [ donor_id, normal_sample_id, tumor_sample_id, gripss_filtered_file ]
        }
      integrated_cnv_files = integrated_sv_files.filter{ it[3] =~ /.integrated.cnvs.txt/ }
      integrated_sv_filtered_files = integrated_sv_files.filter{ it[3] =~ /.integrated.svs.filtered.vcf/ }
    }

    input_create_sv_plots = cobalt_filtered_readcounts_100kb_files
      .combine( cobalt_filtered_readcounts_1mb_files, by: [0,1,2] )
      .combine( cobalt_filtered_readcounts_segments_txt_files, by: [0,1,2] )
      .combine( baf_binned_100kb_files, by: [0,1,2] )
      .combine( baf_filtered_segments_txt_files, by: [0,1,2] )
      .combine( integrated_cnv_files, by: [0,1,2] )

    CreateSvPlots( input_create_sv_plots )
    sv_plots = CreateSvPlots.out
    .transpose()
    .map{
      donor_id, normal_sample_id, tumor_sample_id, sv_plot ->
      file_name = sv_plot.getName()
      sv_plot.copyTo("${params.out_dir}/intermediate/svs/Plots/${donor_id}/${normal_sample_id}/${file_name}")
      [ donor_id, normal_sample_id, tumor_sample_id, sv_plot ]
    }

    input_create_circos_config = integrated_cnv_files
      .combine( cobalt_filtered_readcounts_1mb_files, by: [0,1,2] )
      .combine( baf_binned_1mb_files, by: [0,1,2] )
      .combine( cobalt_filtered_readcounts_segments_txt_files, by: [0,1,2] )
      .combine( baf_filtered_segments_txt_files, by: [0,1,2] )
      .combine( integrated_sv_filtered_files, by: [0,1,2] )

    CreateCircosConfig( input_create_circos_config )
    circos_files = CreateCircosConfig.out
    .transpose()
    .map{
      donor_id, normal_sample_id, tumor_sample_id, config_file ->
      file_name = config_file.getName()
      config_file.copyTo("${params.out_dir}/intermediate/svs/Circos/configs/${donor_id}/${normal_sample_id}/${file_name}")
      [ donor_id, normal_sample_id, tumor_sample_id, config_file ]
    }

    circos_config_files = circos_files.filter{ it[3] =~ /.circos.*.conf/ }
    circos_txt_files = circos_files.filter{ it[3] =~ /.circos.*.txt/ }
    circos_input_files = circos_config_files
      .combine(circos_txt_files, by: [0,1,2] )
      .groupTuple( by: [0,1,2,3 ])

    Circos( circos_input_files )
    circos_plots = Circos.out
    .transpose()
    .map{
      donor_id, normal_sample_id, tumor_sample_id, circos_plot ->
      file_name = circos_plot.getName()
      circos_plot.copyTo("${params.out_dir}/intermediate/svs/Circos/plots/${donor_id}/${normal_sample_id}/${file_name}")
      [ donor_id, normal_sample_id, tumor_sample_id, circos_plot ]
    }
}

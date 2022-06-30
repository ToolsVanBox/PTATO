// def check_params( input_params ) {
//   // run_template_file = file('/hpc/pmc_vanboxtel/tools/ToolsVanBox_DEV/PTATO_DEV/configs/run-template.config')
//
//   // println run_template_file.text
//   def var = 'run'
//   println input_params[var].keySet()
//
//   for var in input_params.keySet() {
//     println var
//   }
// }

// def params_expected = [
//   'run': [ 'snvs', 'QC', 'svs', 'indels' ]
//   'train': 'version'
// ]
//
  // input_vcfs_dir = ''
  // bams_dir = ''
  // out_dir = ''
  // bulk_names = [
  //   ['donor_id', 'sample_id'],
  // ]
  //
  // snvs {
  //   rf_rds = ''
  // }
  //
  // indels {
  //   rf_rds = ''
  //   pon = ''
  // }
  // optional {
  //
  //   germline_vcfs_dir = ''
  //
  //   short_variants {
  //     somatic_vcfs_dir = ''
  //     walker_vcfs_dir = ''
  //     phased_vcfs_dir = ''
  //     ab_tables_dir = ''
  //     context_beds_dir = ''
  //     features_beds_dir = ''
  //   }
  //
  //   snvs {
  //     rf_tables_dir = ''
  //     ptato_vcfs_dir = ''
  //   }
  //
  //   indels {
  //     rf_tables_dir = ''
  //     ptato_vcfs_dir = ''
  //   }
  //
  //   qc {
  //     wgs_metrics_dir = ''
  //     alignment_summary_metrics_dir = ''
//     }
//
//     svs {
//       gridss_driver_vcfs_dir = ''
//       gridss_unfiltered_vcfs_dir = ''
//       gripss_somatic_filtered_vcfs_dir = ''
//       cobalt_ratio_tsv_dir = ''
//
//       cobalt_filtered_readcounts_dir = ''
//       baf_filtered_files_dir = ''
//       gripss_filtered_files_dir = ''
//       integrated_sv_files_dir = ''
//     }
//   }
//
//
// }

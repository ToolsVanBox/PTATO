include { train_indel_rf } from '../../NextflowModules/Utils/randomForest.nf' params(params)

workflow get_indel_rf_files {
  take:
    rf_tables
    label_info
  main:
    input_train_indel_rf = rf_tables
        .combine( label_info, by: [0,1] )
        .map{
          donor_id, sample_id, rf_table, label ->
          [ label, rf_table ]
        }
        .groupTuple( by: [0] )
        .collect()

    train_indel_rf( input_train_indel_rf )

    randomforest_files = train_indel_rf.out
      .map{
        confusion_file, importance_file, rdata_file, rds_file ->
        confusion_name = confusion_file.getName()
        importance_name = importance_file.getName()
        rdata_name = rdata_file.getName()
        rds_name = rds_file.getName()
        confusion_file = confusion_file.copyTo("${params.out_dir}/indels/randomforest/${confusion_name}")
        importance_file = importance_file.copyTo("${params.out_dir}/indels/randomforest/${importance_name}")
        rdata_file = rdata_file.copyTo("${params.out_dir}/indels/randomforest/${rdata_name}")
        rds_file = rds_file.copyTo("${params.out_dir}/indels/randomforest/${rds_name}")
        [ confusion_file, importance_file, rdata_file, rds_file ]
      }

  emit:
      randomforest_files
}
//
// workflow test {
//   take:
//     ab_tables
//     features_beds
//     somatic_vcfs
//
//   main:
//     if ( params.optional.rf_tables_dir ) {
//       rf_tables = extractRfTableRdsFromDir( params.optional.rf_tables_dir )
//     } else {
//       createRFtable( ab_tables.combine( features_beds, by: [0,1] ) )
//       rf_tables = createRFtable.out
//         .map{
//           donor_id, sample_id, rf_table_rds ->
//           rds_name = rf_table_rds.getName()
//           rf_table_rds = rf_table_rds.copyTo("${params.out_dir}/snvs/intermediate/rf/${donor_id}/${rds_name}")
//           [ donor_id, sample_id, rf_table_rds ]
//         }
//     }
//
//     input_testRF = somatic_vcfs.combine( rf_tables, by: [0,1] )
//     testRF( input_testRF )
//
//     bgzip( testRF.out )
//     tabix( bgzip.out )
//     pap_vcfs = tabix.out
//       .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
//         vcf_name = vcf_gz.getName()
//         tbi_name = vcf_tbi.getName()
//         vcf_gz = vcf_gz.copyTo("${params.out_dir}/snvs/${donor_id}/${vcf_name}")
//         vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/snvs/${donor_id}/${tbi_name}")
//         [ donor_id, sample_id, vcf_gz, vcf_tbi ]
//       }
//   emit:
//       pap_vcfs
// }

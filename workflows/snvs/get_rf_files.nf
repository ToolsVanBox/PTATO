include { train_snv_rf } from '../../NextflowModules/Utils/randomForest.nf' params(params)

workflow get_snv_rf_files {
  take:
    rf_tables
    label_info
  main:
    input_train_snv_rf = rf_tables
        .combine( label_info, by: [0,1] )
        .map{
          donor_id, sample_id, rf_table, label ->
          [ label, rf_table ]
        }
        .groupTuple( by: [0] )
        .collect()

    train_snv_rf( input_train_snv_rf )

    randomforest_files = train_snv_rf.out
      .map{
        confusion_file, importance_file, rdata_file, rds_file ->
        confusion_name = confusion_file.getName()
        importance_name = importance_file.getName()
        rdata_name = rdata_file.getName()
        rds_name = rds_file.getName()
        confusion_file = confusion_file.copyTo("${params.out_dir}/snvs/randomforest/${confusion_name}")
        importance_file = importance_file.copyTo("${params.out_dir}/snvs/randomforest/${importance_name}")
        rdata_file = rdata_file.copyTo("${params.out_dir}/snvs/randomforest/${rdata_name}")
        rds_file = rds_file.copyTo("${params.out_dir}/snvs/randomforest/${rds_name}")
        [ confusion_file, importance_file, rdata_file, rds_file ]
      }

  emit:
      randomforest_files
}

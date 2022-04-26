include { closest } from '../../NextflowModules/bedtools/2.30.0/closest.nf' params(params)
include {
  intersect;
  intersectAll
} from '../../NextflowModules/bedtools/2.30.0/intersect.nf' params(params)
include {
  merge;
  mergeAll
} from '../../NextflowModules/bedtools/2.30.0/merge.nf' params(params)


workflow closest_feature {
  take:
    input_sample_beds
    input_feature_beds
  main:
    input_files = input_sample_beds.combine( input_feature_beds )
    closest( input_files )
    merge( closest.out )
    features_beds = merge.out
      .map{ donor_id, sample_id, bed ->
        bed_name = bed.getName()
        bed = bed.copyTo("${params.out_dir}/intermediate/features/${donor_id}/${sample_id}/${bed_name}")
        [ donor_id, sample_id, bed ]
      }

  emit:
    features_beds
}

workflow intersect_feature {
  take:
    input_sample_beds
    input_feature_beds
  main:
    input_files = input_sample_beds.combine(input_feature_beds)
    intersect( input_files )
    merge( intersect.out )
    features_beds = merge.out
      .map{ donor_id, sample_id, bed ->
        bed_name = bed.getName()
        bed = bed.copyTo("${params.out_dir}/intermediate/short_variants/features/${donor_id}/${sample_id}/${bed_name}")
        [ donor_id, sample_id, bed ]
      }
  emit:
    features_beds
}

workflow merge_features {
  take:
    input_sample_beds
    input_feature_merged_beds
  main:
    input_files = input_sample_beds.join(input_feature_merged_beds, by: [0,1] )
    intersectAll( input_files )
    mergeAll( intersectAll.out )
    features_beds = mergeAll.out
      .map{ donor_id, sample_id, bed ->
        bed_name = bed.getName()
        bed = bed.copyTo("${params.out_dir}/intermediate/short_variants/features/${donor_id}/${bed_name}")
        [ donor_id, sample_id, bed ]
      }
  emit:
    features_beds
}

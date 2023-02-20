include { closest } from '../../NextflowModules/bedtools/2.30.0/closest.nf' params(params)
include {
  intersect;
  intersectAll
} from '../../NextflowModules/bedtools/2.30.0/intersect.nf' params(params)
include {
  groupby;
  groupbyAll
} from '../../NextflowModules/bedtools/2.30.0/groupby.nf' params(params)


workflow closest_feature {
  take:
    input_sample_beds
    input_feature_beds
  main:
    input_files = input_sample_beds.combine( input_feature_beds )
    closest( input_files )
    groupby( closest.out )
    features_beds = groupby.out
      .map{ donor_id, sample_id, bed ->
        bed_name = bed.getName()
        bed = bed.copyTo("${params.out_dir}/intermediate/short_variants/features/${donor_id}/${sample_id}/${bed_name}")
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
    groupby( intersect.out )
    features_beds = groupby.out
      .map{ donor_id, sample_id, bed ->
        bed_name = bed.getName()
        bed = bed.copyTo("${params.out_dir}/intermediate/short_variants/features/${donor_id}/${sample_id}/${bed_name}")
        [ donor_id, sample_id, bed ]
      }
  emit:
    features_beds
}

workflow groupby_features {
  take:
    input_sample_beds
    input_feature_groupby_beds
  main:
    input_files = input_sample_beds.join(input_feature_groupby_beds, by: [0,1] )
    intersectAll( input_files )
    groupbyAll( intersectAll.out )
    features_beds = groupbyAll.out
      .map{ donor_id, sample_id, bed ->
        bed_name = bed.getName()
        bed = bed.copyTo("${params.out_dir}/intermediate/short_variants/features/${donor_id}/${bed_name}")
        [ donor_id, sample_id, bed ]
      }
  emit:
    features_beds
}

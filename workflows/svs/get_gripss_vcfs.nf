include { GripssApplicationKt } from '../../NextflowModules/gridss-purple-linx/1.3.2/GripssApplicationKt.nf' params(params)
include { GripssHardFilterApplicationKt } from '../../NextflowModules/gridss-purple-linx/1.3.2/GripssHardFilterApplicationKt.nf' params(params)

workflow get_gripss_vcfs {
  take:
    gridss_unfiltered_vcfs
  main:
    GripssApplicationKt( gridss_unfiltered_vcfs )
    gripss_somatic_vcfs = GripssApplicationKt.out
    .map{
      donor_id, normal_sample_id, tumor_sample_id, gripss_somatic_vcf, gripss_somatic_tbi->
      vcf_name = gripss_somatic_vcf.getName()
      tbi_name = gripss_somatic_tbi.getName()
      gridss_driver_vcf = gripss_somatic_vcf.copyTo("${params.out_dir}/svs/intermediate/gripss/${donor_id}/${normal_sample_id}/${vcf_name}")
      gripss_somatic_tbi = gripss_somatic_tbi.copyTo("${params.out_dir}/svs/intermediate/gripss/${donor_id}/${normal_sample_id}/${tbi_name}")
      [ donor_id, normal_sample_id, tumor_sample_id, gripss_somatic_vcf, gripss_somatic_tbi ]
    }
    GripssHardFilterApplicationKt( gripss_somatic_vcfs )
    gripss_somatic_filtered_vcfs = GripssHardFilterApplicationKt.out
    .map{
      donor_id, normal_sample_id, tumor_sample_id, gripss_somatic_filtered_vcf, gripss_somatic_filtered_tbi->
      vcf_name = gripss_somatic_filtered_vcf.getName()
      tbi_name = gripss_somatic_filtered_tbi.getName()
      gripss_somatic_filtered_vcf = gripss_somatic_filtered_vcf.copyTo("${params.out_dir}/svs/intermediate/gripss/${donor_id}/${normal_sample_id}/${vcf_name}")
      gripss_somatic_filtered_tbi = gripss_somatic_filtered_tbi.copyTo("${params.out_dir}/svs/intermediate/gripss/${donor_id}/${normal_sample_id}/${tbi_name}")
      [ donor_id, normal_sample_id, tumor_sample_id, gripss_somatic_filtered_vcf, gripss_somatic_filtered_tbi ]
    }
  emit:
    gripss_somatic_filtered_vcfs
}

def extractNoptaVcfGzFromDir( nopta_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    nopta_dir = nopta_dir.tokenize().collect{"$it/*/*.vcf.gz"}
    Channel
      .fromPath(nopta_dir, type:'file')
      .ifEmpty { error "No .vcf.gz files found in ${nopta_dir}." }
      .map { nopta_vcf_path ->
          nopta_vcf_file = nopta_vcf_path
          nopta_tbi_file = nopta_vcf_path+".tbi"
          nopta_sample_id = nopta_vcf_path.getName().toString().replaceAll(/.vcf.gz$/, '')
          file(nopta_vcf_path.getBaseName()).getBaseName()
          nopta_donor_id = nopta_vcf_path.getParent().getName()
          [nopta_donor_id, nopta_sample_id, nopta_vcf_file, nopta_tbi_file]
      }
}

def extractPtaVcfGzFromDir( pta_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    pta_dir = pta_dir.tokenize().collect{"$it/*/*.vcf.gz"}
    Channel
      .fromPath(pta_dir, type:'file')
      .ifEmpty { error "No .vcf.gz files found in ${pta_dir}." }
      .map { pta_vcf_path ->
          pta_vcf_file = pta_vcf_path
          pta_tbi_file = pta_vcf_path+".tbi"
          pta_sample_id = pta_vcf_path.getName().toString().replaceAll(/.vcf.gz$/, '')
          pta_donor_id =  pta_vcf_path.getParent().getName()
          [pta_donor_id, pta_sample_id, pta_vcf_file, pta_tbi_file]
      }
}

def extractInputVcfGzFromDir( input_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    input_dir = input_dir.tokenize().collect{"$it/*/*.vcf.gz"}
    Channel
      .fromPath(input_dir, type:'file')
      .ifEmpty { error "No .vcf.gz files found in ${input_dir}." }
      .map { input_vcf_path ->
          input_vcf_file = input_vcf_path
          input_tbi_file = input_vcf_path+".tbi"
          input_sample_id = input_vcf_path.getName().toString().replaceAll(/.vcf.gz$/, '')
          input_donor_id =  input_vcf_path.getParent().getName()
          [input_donor_id, input_sample_id, input_vcf_file, input_tbi_file]
      }
}

def extractGermlineVcfGzFromDir( germline_vcfs_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    germline_vcfs_dir = germline_vcfs_dir.tokenize().collect{"$it/*/*vcf.gz"}
    Channel
      .fromPath(germline_vcfs_dir, type:'file')
      .ifEmpty { error "No .vcf.gz files found in ${germline_vcfs_dir}." }
      .map { germline_vcf_path ->
          germline_vcf_file = germline_vcf_path
          germline_tbi_file = germline_vcf_path+".tbi"
          germline_sample_id = germline_vcf_path.getName().toString().replaceAll(/(.germline)*.vcf.gz$/, '')
          germline_donor_id =  germline_vcf_path.getParent().getName()
          [germline_donor_id, germline_sample_id, germline_vcf_file, germline_tbi_file]
      }
}

def extractSomaticVcfGzFromDir( somatic_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    somatic_dir = somatic_dir.tokenize().collect{"$it/*/*.vcf.gz"}
    Channel
      .fromPath(somatic_dir, type:'file')
      .ifEmpty { error "No .vcf.gz files found in ${somatic_dir}." }
      .map { somatic_vcf_path ->
          somatic_vcf_file = somatic_vcf_path
          somatic_tbi_file = somatic_vcf_path+".tbi"
          somatic_sample_id = somatic_vcf_path.getName().toString().replaceAll(/(.SMuRF.filtered)*.vcf.gz$/, '')
          somatic_donor_id =  somatic_vcf_path.getParent().getName()
          [somatic_donor_id, somatic_sample_id, somatic_vcf_file, somatic_tbi_file]
      }
}

def extractBamsFromDir( bams_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    bams_dir = bams_dir.tokenize().collect{"$it/*/*.bam"}
    Channel
      .fromPath(bams_dir, type:'file')
      .ifEmpty { error "No bam files found in ${bams_dir}." }
      .map { bam_path ->
          bam_file = bam_path
          bai_file = bam_path.toString().replace('.bam','.bai')
          bam_sample_id = bam_path.getName().toString().replaceAll(/(_dedup)*.bam$/, '')
          bam_donor_id =  bam_path.getParent().getName()
          [bam_donor_id, bam_sample_id, bam_file, bai_file]
      }
}

def extractWalkerVcfGzFromDir( walker_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  walker_dir = walker_dir.tokenize().collect{"$it/*/*.walker.vcf.gz"}
  Channel
    .fromPath(walker_dir, type:'file')
    .ifEmpty { error "No .walker.vcf.gz files found in ${walker_dir}." }
    .map { walker_vcf_path ->
        walker_vcf_file = walker_vcf_path
        walker_tbi_file = walker_vcf_path+".tbi"
        walker_sample_id = walker_vcf_path.getName().toString().replaceAll(/.walker.vcf.gz$/, '')
        walker_donor_id =  walker_vcf_path.getParent().getName()
        [walker_donor_id, walker_sample_id, walker_vcf_file, walker_tbi_file]
    }
}

def extractAbTableFromDir( ab_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  ab_dir = ab_dir.tokenize().collect{"$it/*/*.abtable.txt"}
  Channel
    .fromPath(ab_dir, type:'file')
    .ifEmpty { error "No .abtable.txt files found in ${ab_dir}." }
    .map { ab_path ->
        ab_file = ab_path
        ab_sample_id = ab_path.getName().toString().replaceAll(/.abtable.txt$/, '')
        ab_donor_id =  ab_path.getParent().getName()
        [ab_donor_id, ab_sample_id, ab_file]
    }
}

def extractContextBedFromDir( context_bed_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  context_bed_dir = context_bed_dir.tokenize().collect{"$it/*/*.bed"}
  Channel
    .fromPath(context_bed_dir, type:'file')
    .ifEmpty { error "No .bed files found in ${context_bed_dir}." }
    .map { context_bed_path ->
        context_bed_file = context_bed_path
        context_sample_id = context_bed_path.getName().toString().replaceAll(/(.context)*(.sorted)*.bed$/, '')
        context_donor_id =  context_bed_path.getParent().getName()
        [context_donor_id, context_sample_id, context_bed_file]
    }
}

def extractFeaturesBedFromDir( features_bed_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  features_bed_dir = features_bed_dir.tokenize().collect{"$it/*/*.bed"}
  Channel
    .fromPath(features_bed_dir, type:'file')
    .ifEmpty { error "No .bed files found in ${features_bed_dir}." }
    .map { features_bed_path ->
        features_bed_file = features_bed_path
        features_sample_id = features_bed_path.getName().toString().replaceAll(/(.features)*(.merged)*.bed$/, '')
        features_donor_id =  features_bed_path.getParent().getName()
        [features_donor_id, features_sample_id, features_bed_file]
    }
}

def extractPhasedVcfGzFromDir( phased_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  phased_dir = phased_dir.tokenize().collect{"$it/*/*.phased.vcf.gz"}
  Channel
    .fromPath(phased_dir, type:'file')
    .ifEmpty { error "No .phased.vcf.gz files found in ${phased_dir}." }
    .map { phased_vcf_path ->
        phased_vcf_file = phased_vcf_path
        phased_tbi_file = phased_vcf_path+".tbi"
        phased_sample_id = phased_vcf_path.getName().toString().replaceAll(/\.(.+?)\.phased.vcf.gz$/, '')
        phased_donor_id =  phased_vcf_path.getParent().getName()
        m = phased_vcf_path =~ /.+\.(.+?)\.phased.vcf.gz/
        phased_chrom = m[0][1]
        [phased_donor_id, phased_sample_id, phased_chrom, phased_vcf_file, phased_tbi_file]
    }
}

def extractSnvRfTableRdsFromDir( snv_rf_tables_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  snv_rf_tables_dir = snv_rf_tables_dir.tokenize().collect{"$it/*/*.rftable.rds"}
  Channel
    .fromPath(snv_rf_tables_dir, type:'file')
    .ifEmpty { error "No .rftable.rds files found in ${snv_rf_tables_dir}." }
    .map { snv_rf_rds_path ->
        snv_rf_rds_file = snv_rf_rds_path
        snv_rf_rds_sample_id = snv_rf_rds_path.getName().toString().replaceAll(/.rftable.rds$/, '')
        snv_rf_rds_donor_id =  snv_rf_rds_path.getParent().getName()
        [snv_rf_rds_donor_id, snv_rf_rds_sample_id, snv_rf_rds_file]
    }
}

def extractIndelRfTableRdsFromDir( indel_rf_tables_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  indel_rf_tables_dir = indel_rf_tables_dir.tokenize().collect{"$it/*/*.rftable.rds"}
  Channel
    .fromPath(indel_rf_tables_dir, type:'file')
    .ifEmpty { error "No .rftable.rds files found in ${indel_rf_tables_dir}." }
    .map { indel_rf_rds_path ->
        indel_rf_rds_file = indel_rf_rds_path
        indel_rf_rds_sample_id = indel_rf_rds_path.getName().toString().replaceAll(/.rftable.rds$/, '')
        indel_rf_rds_donor_id =  indel_rf_rds_path.getParent().getName()
        [indel_rf_rds_donor_id, indel_rf_rds_sample_id, indel_rf_rds_file]
    }
}

def extractPtatoVcfGzFromDir( ptato_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  ptato_vcfs_dir = ptato_vcfs_dir.tokenize().collect{"$it/*/*.ptato.vcf.gz"}
  Channel
    .fromPath(ptato_vcfs_dir, type:'file')
    .ifEmpty { error "No .ptato.vcf.gz files found in ${ptato_vcfs_dir}." }
    .map { ptato_vcf_path ->
        ptato_vcf_file = ptato_vcf_path
        ptato_tbi_file = ptato_vcf_path+".tbi"
        ptato_sample_id = ptato_vcf_path.getName().toString().replaceAll(/.ptato.vcf.gz$/, '')
        ptato_donor_id =  ptato_vcf_path.getParent().getName()
        [ptato_donor_id, ptato_sample_id, ptato_vcf_file, ptato_tbi_file]
    }
}

def extractWGSMetricsFromDir( wgs_metrics_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  wgs_metrics_dir = wgs_metrics_dir.tokenize().collect{"$it/*/*.wgs_metrics*"}
  Channel
    .fromPath(wgs_metrics_dir, type:'file')
    .ifEmpty { error "No .wgs_metrics* files found in ${wgs_metrics_dir}." }
    .map { wgs_metrics_path ->
        wgs_metrics_file = wgs_metrics_path
        wgs_metrics_sample_id = wgs_metrics_path.getName().toString().replaceAll(/.wgs_metrics.*$/, '')
        wgs_metrics_donor_id =  wgs_metrics_path.getParent().getName()
        [wgs_metrics_donor_id, wgs_metrics_sample_id, wgs_metrics_file]
    }
}

def extractAlignmentSummaryMetricsFromDir( alignment_summary_metrics_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  alignment_summary_metrics_dir = alignment_summary_metrics_dir.tokenize().collect{"$it/*/*.alignment_summary_metrics*"}
  Channel
    .fromPath(alignment_summary_metrics_dir, type:'file')
    .ifEmpty { error "No .alignment_summary_metrics* files found in ${alignment_summary_metrics_dir}." }
    .map { alignment_summary_metrics_path ->
        alignment_summary_metrics_file = alignment_summary_metrics_path
        alignment_summary_metrics_sample_id = alignment_summary_metrics_path.getName().toString().replaceAll(/(.multiple_metrics)*.alignment_summary_metrics.*$/, '')
        alignment_summary_metrics_donor_id =  alignment_summary_metrics_path.getParent().getName()
        [alignment_summary_metrics_donor_id, alignment_summary_metrics_sample_id, alignment_summary_metrics_file]
    }
}

def extractGridssDriverVcfGzFromDir( gridss_driver_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  gridss_driver_vcfs_dir = gridss_driver_vcfs_dir.tokenize().collect{"$it/*/*.gridss.driver.vcf.gz"}
  Channel
    .fromPath(gridss_driver_vcfs_dir, type:'file')
    .ifEmpty { error "No .gridss.driver.vcf.gz files found in ${gridss_driver_vcfs_dir}." }
    .map { gridss_driver_vcf_path ->
        gridss_driver_vcf_file = gridss_driver_vcf_path
        gridss_driver_tbi_file = gridss_driver_vcf_path+".tbi"
        gridss_driver_tumor_sample_id = gridss_driver_vcf_path.getName().toString().replaceAll(/.gridss.driver.vcf.gz$/, '')
        gridss_driver_normal_sample_id =  gridss_driver_vcf_path.getParent().getName()
        gridss_driver_donor_id =  gridss_driver_vcf_path.getParent().getParent().getName()
        [gridss_driver_donor_id, gridss_driver_normal_sample_id, gridss_driver_tumor_sample_id, gridss_driver_vcf_file, gridss_driver_tbi_file]
    }
}


def extractGridssUnfilteredVcfGzFromDir( gridss_unfiltered_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  gridss_unfiltered_vcfs_dir = gridss_unfiltered_vcfs_dir.tokenize().collect{"$it/*/*.gridss.unfiltered.vcf.gz"}
  Channel
    .fromPath(gridss_unfiltered_vcfs_dir, type:'file')
    .ifEmpty { error "No .gridss.unfiltered.vcf.gz files found in ${gridss_unfiltered_vcfs_dir}." }
    .map { gridss_unfiltered_vcf_path ->
        gridss_unfiltered_vcf_file = gridss_unfiltered_vcf_path
        gridss_unfiltered_tbi_file = gridss_unfiltered_vcf_path+".tbi"
        gridss_unfiltered_tumor_sample_id = gridss_unfiltered_vcf_path.getName().toString().replaceAll(/.gridss.unfiltered.vcf.gz$/, '')
        gridss_unfiltered_normal_sample_id =  gridss_unfiltered_vcf_path.getParent().getName()
        gridss_unfiltered_donor_id =  gridss_unfiltered_vcf_path.getParent().getParent().getName()
        [gridss_unfiltered_donor_id, gridss_unfiltered_normal_sample_id, gridss_unfiltered_tumor_sample_id, gridss_unfiltered_vcf_file, gridss_unfiltered_tbi_file]
    }
}

def extractGripssSomaticFilteredVcfGzFromDir( gripss_somatic_filtered_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  gripss_somatic_filtered_vcfs_dir = gripss_somatic_filtered_vcfs_dir.tokenize().collect{"$it/*/*.gripss.somatic.filtered.vcf.gz"}
  Channel
    .fromPath(gripss_somatic_filtered_vcfs_dir, type:'file')
    .ifEmpty { error "No .gripss.somatic.filtered.vcf.gz files found in ${gripss_somatic_filtered_vcfs_dir}." }
    .map { gripss_somatic_filtered_vcf_path ->
        gripss_somatic_filtered_vcf_file = gripss_somatic_filtered_vcf_path
        gripss_somatic_filtered_tbi_file = gripss_somatic_filtered_vcf_path+".tbi"
        gripss_somatic_filtered_tumor_sample_id = gripss_somatic_filtered_vcf_path.getName().toString().replaceAll(/.gripss.somatic.filtered.vcf.gz$/, '')
        gripss_somatic_filtered_normal_sample_id =  gripss_somatic_filtered_vcf_path.getParent().getName()
        gripss_somatic_filtered_donor_id =  gripss_somatic_filtered_vcf_path.getParent().getParent().getName()
        [gripss_somatic_filtered_donor_id, gripss_somatic_filtered_normal_sample_id, gripss_somatic_filtered_tumor_sample_id, gripss_somatic_filtered_vcf_file, gripss_somatic_filtered_tbi_file]
    }
}

def extractCobaltRatioTsvFromDir( cobalt_ratio_tsv_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  cobalt_ratio_tsv_dir = cobalt_ratio_tsv_dir.tokenize().collect{"$it/*/*.cobalt.ratio.tsv"}
  Channel
    .fromPath(cobalt_ratio_tsv_dir, type:'file')
    .ifEmpty { error "No .cobalt.ratio.tsv files found in ${cobalt_ratio_tsv_dir}." }
    .map { cobalt_ratio_tsv_path ->
        cobalt_ratio_tsv_file = cobalt_ratio_tsv_path
        cobalt_ratio_tumor_sample_id = cobalt_ratio_tsv_path.getName().toString().replaceAll(/.cobalt.ratio.tsv$/, '')
        cobalt_ratio_normal_sample_id =  cobalt_ratio_tsv_path.getParent().getName()
        cobalt_ratio_donor_id =  cobalt_ratio_tsv_path.getParent().getParent().getName()
        [cobalt_ratio_donor_id, cobalt_ratio_normal_sample_id, cobalt_ratio_tumor_sample_id, cobalt_ratio_tsv_file]
    }
}

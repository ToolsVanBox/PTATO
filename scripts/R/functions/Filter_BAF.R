### TO DO
## Make a excludebaflist, like the copynumber / COBALT excludelist
# Or, overlap the BAF with the COBALT excludelist
## Fine-map the CNVs. First run it at 100k, then finemap the edges by running it at 1-10kb
## Filter SNPs next to centromeres

Filter_BAF <- function(SNV_VCF_File, OUTPUT_dir, SAMPLE, NORMAL,
                       CENTROMERES, CYTOBAND, MIN_DEPTH_SAMPLE = 5, MIN_DEPTH_NORMAL = 10, CHR_SIZES){

  BAF_Output_File <-  paste(OUTPUT_dir, ".baf.filtered.txt", sep = "")

  print(paste("# Reading: ", SNV_VCF_File, sep = ""))
  SNV_VCF <- Read_SNV_VCF(VCF_file = SNV_VCF_File, SAMPLE = SAMPLE, NORMAL = NORMAL)
  print("# Filtering VAFs")
  SNV_VAFs <- Filter_VAFs(SNVs = SNV_VCF, SAMPLE = SAMPLE, NORMAL = NORMAL, CENTROMERES = CENTROMERES,
                          MIN_DEPTH_SAMPLE = MIN_DEPTH_SAMPLE,
                          MIN_DEPTH_NORMAL = MIN_DEPTH_NORMAL)

  #print(head(SNV_VAFs))

  print("# Binning VAFs")
  BAF_binned <- Merge_VAF_Bins(SNV_VAFs,
                               Binsize = 100000, Chr_sizes = CHR_SIZES, Value_column = which(names(SNV_VAFs) == "VAF_normalized"))

  #print("# 1 MB")
  BAF_binned_1mb <- Merge_VAF_Bins(SNV_VAFs,
                               Binsize = 1000000, Chr_sizes = CHR_SIZES, Value_column = which(names(SNV_VAFs) == "VAF_normalized"))

  print("# Segmenting VAF bins")
  BAF_Segments <- Segment_BAFs(BAF = BAF_binned[,c(1,2,4)], Binsize = 100000, Cytoband = CYTOBAND)

  print("# Finemapping VAF segments")
  BAF_Segments_Finemapped <- Finemap_BAF_Segments(SegmentsCoarse = BAF_Segments[,c(2, 4,5,7)],
                                                  SegmentsFine = SNV_VAFs[,c(1,2,2, 7)],
                                                  Cutoff_Low = 0.16,
                                                  Cutoff_High = 0.4,
                                                  Max_Finemapping_Dist = 200000)

  print(paste("# Writing: ", OUTPUT_dir, ".baf.segments.raw.txt", sep = ""))
  write.table(BAF_Segments, file = paste(OUTPUT_dir, ".baf.segments.raw.txt", sep = ""), quote = F, row.names = F, sep = "\t")

  print(paste("# Writing: ", OUTPUT_dir, ".baf.segments.txt", sep = ""))
  write.table(BAF_Segments_Finemapped, file = paste(OUTPUT_dir, ".baf.segments.txt", sep = ""), quote = F, row.names = F, sep = "\t")

  print(paste("# Writing: ", OUTPUT_dir, ".baf.binned.100kb.txt", sep = ""))
  write.table(BAF_binned,
              file = paste(OUTPUT_dir, ".baf.binned.100kb.txt", sep = ""), quote = F, row.names = F, sep = "\t")

  print(paste("# Writing: ", OUTPUT_dir, ".baf.binned.1mb.txt", sep = ""))
  write.table(BAF_binned_1mb,
              file = paste(OUTPUT_dir, ".baf.binned.1mb.txt", sep = ""), quote = F, row.names = F, sep = "\t")

  #print(paste("# Writing: ", OUTPUT_dir, "BAF_segments.bedpe", sep = ""))
  ### bedpe files can be imported by the the StructuralVariantAnnotation package. df_to_bedpe in "Filter_COBALT.R" script
  BAF_bedpe <- .df_to_bedpe(BAF_Segments_Finemapped, score = "mean", output_file = paste(OUTPUT_dir, ".baf.segments.bedpe", sep = ""))

  print(paste("# Writing: ", BAF_Output_File, sep = ""))
  write.table(SNV_VAFs,file = BAF_Output_File, quote = F, row.names = F, sep = "\t")

  return(SNV_VAFs)
}

Read_SNV_VCF <- function(VCF_file, SAMPLE, NORMAL){
  param <- ScanVcfParam(samples = c(SAMPLE, NORMAL), info=NA)

  ### Read and filter the GATK/SMuRF vcf in chunks. Still takes minutes
  SNVs_HET <- NULL
  tab <- VcfFile(VCF_file, yieldSize=100000)
  open(tab)
  while (nrow(vcf_yield <- readVcf(tab, "hg38", param=param))){
    cat("vcf dim:", dim(vcf_yield), "\n")

    SNVs <- vcf_yield[isSNV(vcf_yield),] # for now only select the SNVs, as indel VAF may be inaccurate
    SNVs_yield <- SNVs[which(geno(SNVs)$GT[,NORMAL] == "0/1" | geno(SNVs)$GT[,NORMAL] == "0|1"),]
    # Remove low quality SNPs that are marked by SMuRF:
    # SNVs_yield <- SNVs_yield[rowRanges(SNVs_yield)$FILTER %in% c("KnownVariant", "ControlEvidence")]
    SNVs_yield <- SNVs_yield[rowRanges(SNVs_yield)$FILTER %in% c("PASS")]
    if(length(SNVs_HET) > 0){
      SNVs_HET <- rbind(SNVs_HET, SNVs_yield)
    } else {
      SNVs_HET <- SNVs_yield
    }
  }
  close(tab)
  SNVs_HET <- SNVs_HET[!duplicated(names(SNVs_HET))]
  return(SNVs_HET)
}

Filter_VAFs <- function(SNVs, SAMPLE, NORMAL, CENTROMERES, MIN_DEPTH_SAMPLE = 5, MIN_DEPTH_NORMAL = 10,
                        MIN_VAF_NORMAL = 0.3, MAX_VAF_NORMAL = 0.7){

  ADs <- geno(SNVs)$AD[names(SNVs), SAMPLE] # VAFs are calculated by dividing: ALT / (REF+ALT) = AD[2] / (AD[1] + AD[2])
  ADs_normal <- geno(SNVs)$AD[names(SNVs), NORMAL]
  # Calculate for each SNP how much the VAF is divergent from 0.5 (eg if the VAF is 0.9, it's 0.4 from 0.5)
  VAFs <- data.frame(row.names = names(unlist(lapply(ADs, function(x) x[2]/(x[1]+x[2])))),
                     chromosome = seqnames(rowRanges(SNVs)),
                     position = start(rowRanges(SNVs)),
                     DP_Sample = geno(SNVs)$DP[names(SNVs), SAMPLE],
                     DP_Normal = geno(SNVs)$DP[names(SNVs), NORMAL],
                     VAF = unlist(lapply(ADs, function(x) x[2]/(x[1]+x[2]))),
                     VAF_Normal = unlist(lapply(ADs_normal, function(x) x[2]/(x[1]+x[2]))))
  VAFs$VAF_normalized <- abs(0.5-VAFs$VAF)

  # Remove the SNPs with loo low coverage in either the sample or the normal
  VAFs_filtered <- VAFs[which(VAFs$DP_Sample > MIN_DEPTH_SAMPLE), ]
  VAFs_filtered <- VAFs_filtered[which(VAFs_filtered$DP_Normal > MIN_DEPTH_NORMAL), ]

  # Remove SNPs on decoy chromosomes and SNPs with divergent VAFs in the normal
  VAFs_filtered <- VAFs_filtered[which(VAFs_filtered$VAF_Normal > MIN_VAF_NORMAL & VAFs_filtered$VAF_Normal < MAX_VAF_NORMAL),]
  VAFs_filtered <- VAFs_filtered[which(VAFs_filtered$chromosome %in% c(1:22, "X", "Y")),]

  # Overlap the SNPs with the centromeres
  VAFs_filtered_g <- GRanges(seqnames = VAFs_filtered$chromosome, IRanges(start = VAFs_filtered$position, end = VAFs_filtered$position))
  overlap_SNPs_centromeres <- findOverlaps(VAFs_filtered_g, CENTROMERES)

  if(length(overlap_SNPs_centromeres) > 0){
    print(paste("# Removed ", length(overlap_SNPs_centromeres), " of ", nrow(VAFs_filtered), " SNPs close to centromeres/telomeres", sep = ""))
    VAFs_filtered <- VAFs_filtered[-queryHits(overlap_SNPs_centromeres),]
  }
  VAFs_filtered <- VAFs_filtered[!is.na(VAFs_filtered$VAF_normalized),]
  return(VAFs_filtered)
}


Merge_VAF_Bins <- function(Input, Value_column = 5, Binsize = 100000, Chr_sizes = NULL) {
  df <- data.frame()
  # Use mixedsort to order the chromosomes alphanumerically (allows different ref genomes, but causes dependency on gtools)
  Chromosomes <- gtools::mixedsort(as.character(unique((Input[,1]))))
  #print(chromosomes)
  for (i in c(1:length(Chromosomes))) {
    Chrom <- Chromosomes[i]

    # Select chromosome
    tmp <- subset( Input, chromosome==Chrom)

    # Binning of positions
    grouping <- cut(tmp$position,
                    c(seq(1,max(tmp$position),Binsize), max(tmp$position)),
                    labels=seq(1,max(tmp$position),Binsize), max(tmp$position))

    # Calculate median readdepth
    medians <- data.frame(tapply(tmp[,Value_column], grouping, median))
    rm(tmp)

    colnames(medians) <- "medianBAF"

    # Create clean data farme
    medians$Chromosome <- Chrom
    medians$start.pos <- as.numeric(as.character(rownames(medians)))
    medians$end.pos <- as.numeric(as.character(rownames(medians)))+Binsize-1
    # Limit end of bins to end of chromosome
    if(!is.null(Chr_sizes)){
      medians$end.pos[medians$end.pos >  as.numeric(Chr_sizes[Chrom])] <- as.numeric(Chr_sizes[Chrom])
    }
    medians <- medians[-which(is.na(medians$medianBAF)),]

    # Merge data
    df <- rbind(df, medians)
  }
  df <- df[,c("Chromosome","start.pos", "end.pos", "medianBAF")]
  return(df)
}

Segment_BAFs <- function(BAF, Binsize, Cytoband = NULL, Winsorize = TRUE){
  # ReadCounts need to be Chromosome, start.pos, value
  if(Winsorize == TRUE){
    BAF.wins <- winsorize(data=BAF,verbose=FALSE) # remove outlier bins
  } else {
    BAF.wins <- BAF
  }

  BAF.wins_g <- GRanges(seqnames = BAF.wins[,1], IRanges(start = BAF.wins[,2], end = BAF.wins[,2]))

  if(!is.null(Cytoband)){
    if(class(Cytoband) == "GRanges"){
      olap_BAF_Cytoband <- nearest(BAF.wins_g, Cytoband)
      ChromosomeArms_Binned <- Cytoband$arm[olap_BAF_Cytoband]
      BAF.seg <- pcf(data=BAF.wins, gamma=100,verbose=FALSE, arms = ChromosomeArms_Binned)
    } else {
      print(paste("!!! Cytoband coordinates not provided", sep = ""))
      print("Using hg19 coordinates to determine chromosome arm positions")
      BAF.seg <- pcf(data=BAF.wins, gamma=100,verbose=FALSE, assembly = "hg19")
    }
  }
  # the CopyNumber package works with only one position per bin (so not start and end, but just start, center or end). Chose for start.pos, and change the end.pos of the segments by subtracting the binsize-1
  BAF.seg$end.pos <- BAF.seg$end.pos + (Binsize-1)

  return(BAF.seg)
}


Finemap_BAF_Segments <- function(SegmentsCoarse, SegmentsFine, Cutoff_Low, Cutoff_High, Max_Finemapping_Dist = 200000){
  # Segments should be chrom - start - end - value
  SegmentsCoarse_g <- GRanges(seqnames = SegmentsCoarse[,1], IRanges(start = SegmentsCoarse[,2], end = SegmentsCoarse[,3]))
  SegmentsFine_g <- GRanges(seqnames = SegmentsFine[,1], IRanges(start = SegmentsFine[,2], end = SegmentsFine[,3]), Value = SegmentsFine[,4])
  Segments_Finemapped <- SegmentsCoarse

  for(i in 1:nrow(Segments_Finemapped)){
    #print(i)
    if(Segments_Finemapped[i, 4] < Cutoff_Low | Segments_Finemapped[i,4 ] > Cutoff_High){
      Event <- GRanges(seqnames = Segments_Finemapped[i,1], IRanges(start = Segments_Finemapped[i,2], end = Segments_Finemapped[i,3]))
      if(Segments_Finemapped[i,4] < Cutoff_Low ){
        Small_segments <- SegmentsFine_g[SegmentsFine_g$Value < Cutoff_Low]
      } else if (Segments_Finemapped[i,4] > Cutoff_High) {
        Small_segments <- SegmentsFine_g[SegmentsFine_g$Value > Cutoff_High]
      }
      olap <- findOverlaps(Event, Small_segments)
      if(length(olap) > 0){
        # Adjust the start pos of the segment (if the overlapping small segments are within MAX_FINEMAPPING_DIST)
        if(abs(Segments_Finemapped[i,2] - min(start(Small_segments[subjectHits(olap)]))) < Max_Finemapping_Dist){
          Segments_Finemapped[i,2] <- min(start(Small_segments[subjectHits(olap)]))
        }
        # Adjust the end pos of the segment (if the overlapping small segments are within MAX_FINEMAPPING_DIST)
        if(abs(Segments_Finemapped[i,3] - max(end(Small_segments[subjectHits(olap)]))) < Max_Finemapping_Dist){
          Segments_Finemapped[i,3] <- max(end(Small_segments[subjectHits(olap)]))
        }
        # Adjust the start pos of the subsequent segment if the segment and the subsequent segment overlap after finemapping (only if it's not the last segment)
        if(i + 1 <= nrow(Segments_Finemapped)){
          if(Segments_Finemapped[i+1,2] <= Segments_Finemapped[i,3] & Segments_Finemapped[i,1] == Segments_Finemapped[i+1,1]){
            Segments_Finemapped[i+1,2] <- Segments_Finemapped[i,3]+1
          }
        }
        # Adjust the end pos of the previous segment if the segment and the previous segment overlap after finemapping (only if it's not the first segment)
        if(i-1 > 0){
          if(Segments_Finemapped[i-1,3] >= Segments_Finemapped[i,2] & Segments_Finemapped[i,1] == Segments_Finemapped[i-1,1]){
            Segments_Finemapped[i-1,3] <- Segments_Finemapped[i,2]-1
          }
        }
      }
    }
  }

  return(Segments_Finemapped)
}

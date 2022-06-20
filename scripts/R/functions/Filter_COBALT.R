### The functions in this script are used to normalize PTA COBALT readcounts against a consensus / PON COBALT readcounts list
### In addition, the normalized read counts (in 1kb bins) are binned in larger bins
### ReadCounts are also segmented (eg bins with similar readcounts are grouped together)

.NormalizeCobalt <- function(COBALT_sample, COBALT_PON, Chromosomes, sample_colum = 4, PON_mean_column = 15){
  
  Autosomes <- as.numeric(names(Chromosomes))[!is.na( as.numeric(names(Chromosomes)))]
  
  # First normalize to the same total number of counts as the pon file (autosomes only)
  normalization_factor <- sum(COBALT_sample[which(COBALT_sample$chromosome %in% Autosomes),sample_colum]) / sum(COBALT_PON[which(COBALT_PON$chromosome %in% Autosomes),PON_mean_column])
  
  COBALT_normalized <- COBALT_sample[,c(1,2, sample_colum)] # assume column 1 and 2 are chromosome and position

  COBALT_normalized[,"normReadCount"] <- COBALT_normalized[,3] / normalization_factor

  # Next the normReadCount are divided by the PON_meanReadCount to normalize the ReadCounts for recurrent PTA coverage deviations
  COBALT_output <- merge(COBALT_normalized, COBALT_PON[,c(1,2, PON_mean_column)], by = c(1,2))
  names(COBALT_output)[ncol(COBALT_output)] <- "PON_meanReadCount"
  COBALT_output$filteredReadCount <- COBALT_output[,"normReadCount"] / COBALT_output[,ncol(COBALT_output)]
  
  # Sex chromosomes
  # The ReadCount PON only contains XY samples. 
  COBALT_output$filteredReadCount[which(COBALT_output$chromosome %in% c("X", "Y"))] <- COBALT_output$filteredReadCount[which(COBALT_output$chromosome %in% c("X", "Y"))]/2
  
  # Chromosomes are ordered based on the order in the chromosome/CHR_SIZES input
  COBALT_output$chromosome <- factor(COBALT_output$chromosome, levels = names(Chromosomes))
  COBALT_output <- COBALT_output[order(COBALT_output$chromosome, COBALT_output$position),]
  
  return(COBALT_output)
}



.RemoveCobaltOutliers <- function(COBALT_input, CENTROMERES = NULL, COBALT_MIN_COV = 0.01, COBALT_MAX_COV = 0.99){
  # Filter the bins with extreme high or low PON_meanReadCount
  COBALT_input$FILTER <- ""
  COBALT_input$FILTER[COBALT_input$PON_meanReadCount < quantile(COBALT_input$PON_meanReadCount, COBALT_MIN_COV)] <- "LOW_MEAN_COUNT"
  COBALT_input$FILTER[COBALT_input$PON_meanReadCount > quantile(COBALT_input$PON_meanReadCount, COBALT_MAX_COV)] <- "HIGH_MEAN_COUNT"

  if(!is.null(CENTROMERES)){
    if(class(CENTROMERES) == "GRanges"){
      # Overlap the COBALT bins with the centromeres
      COBALT_g <- GRanges(seqnames = COBALT_input$chromosome, IRanges(start = COBALT_input$position,
                                                                      end = COBALT_input$position + 999))
      olap_COBALT_CENT <- findOverlaps(COBALT_g, CENTROMERES)

      COBALT_input$FILTER[queryHits(olap_COBALT_CENT)][COBALT_input$FILTER[queryHits(olap_COBALT_CENT)] != ""] <- paste(COBALT_input$FILTER[queryHits(olap_COBALT_CENT)][COBALT_input$FILTER[queryHits(olap_COBALT_CENT)] != ""], "CENTROMERIC", sep = ";")
      COBALT_input$FILTER[queryHits(olap_COBALT_CENT)][COBALT_input$FILTER[queryHits(olap_COBALT_CENT)] == ""] <- "CENTROMERIC"
    }
  }

  COBALT_input$FILTER[COBALT_input$FILTER == ""] <- "PASS"
  return(COBALT_input)
}


.MergeBins <- function(Input, Chr_sizes = NULL, Value_column = 3, Binsize = 100000, Ploidy = 2) {
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

    colnames(medians) <- "medianReadCount"

    # Create clean data farme
    medians$Chromosome <- Chrom
    medians$start.pos <- as.numeric(as.character(rownames(medians)))
    medians$end.pos <- as.numeric(as.character(rownames(medians)))+Binsize-1
    # Limit end of bins to end of chromosome
    if(!is.null(Chr_sizes)){
      medians$end.pos[medians$end.pos >  as.numeric(Chr_sizes[Chrom])] <-  as.numeric(Chr_sizes[Chrom])
    }
    medians$CopyNumber <- medians$medianReadCount*Ploidy
    medians <- medians[-which(is.na(medians$medianReadCount)),]

    # Merge data
    df <- rbind(df, medians)
  }
  df <- df[,c("Chromosome","start.pos", "end.pos", "medianReadCount", "CopyNumber")]
  return(df)
}

.Segment_ReadCounts <- function(ReadCounts, Binsize, Chr_Sizes = NULL, Cytoband = NULL, Winsorize = TRUE){
  # ReadCounts need to be Chromosome, start.pos, value

  if(Winsorize == TRUE){
    ReadCounts.wins <- winsorize(data=ReadCounts,verbose=FALSE) # remove outlier bins
  } else {
    ReadCounts.wins <- ReadCounts
  }

  ReadCounts.wins_g <- GRanges(seqnames = ReadCounts.wins[,1], IRanges(start = ReadCounts.wins[,2], end = ReadCounts.wins[,2]))

  if(!is.null(Cytoband)){
    if(class(Cytoband) == "GRanges"){
      olap_ReadCounts_Cytoband <- nearest(ReadCounts.wins_g, Cytoband)
      ChromosomeArms_Binned <- Cytoband$arm[olap_ReadCounts_Cytoband]
      ReadCounts.seg <- pcf(data=ReadCounts.wins, gamma=100,verbose=FALSE, arms = ChromosomeArms_Binned)
    } else {
      print(paste("!!! Cytoband coordinates not provided", sep = ""))
      print("Using hg19 coordinates to determine chromosome arm positions")
      ReadCounts.seg <- pcf(data=ReadCounts.wins, gamma=100,verbose=FALSE, assembly = "hg19")
    }
  }
  # the CopyNumber package works with only one position per bin (so not start and end, but just start, center or end). Chose for start.pos, and change the end.pos of the segments by subtracting the binsize-1
  ReadCounts.seg$end.pos <- ReadCounts.seg$end.pos + (Binsize-1)

  # Limit end of segments to end of chromosomes
  if(!is.null(Chr_Sizes)){
    for(i in 1:nrow(ReadCounts.seg)){
      if(ReadCounts.seg[i,"end.pos"] > as.numeric(Chr_Sizes[names(Chr_Sizes) == ReadCounts.seg[i,"chrom"]])){
        ReadCounts.seg[i,"end.pos"] <- as.numeric(Chr_Sizes[names(Chr_Sizes) == ReadCounts.seg[i,"chrom"]])
      }
    }
  }

  return(ReadCounts.seg)
}

.FinemapSegments <- function(SegmentsCoarse, SegmentsFine, Cutoff_Low, Cutoff_High, Max_Finemapping_Dist = 200000){
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

FilterCobalt <- function(COBALT_input, COBALT_PON, OUTPUT_dir, CENTROMERES, CYTOBAND, COBALT_MIN_COV, COBALT_MAX_COV, CHR_SIZES){
  print("# Normalizing COBALT ReadCounts")
  NormalizedReadCount <- .NormalizeCobalt(COBALT_sample = COBALT_input, COBALT_PON = COBALT_PON, Chromosomes = CHR_SIZES)
  print("# Filtering COBALT ReadCounts")
  FilteredReadCount <- .RemoveCobaltOutliers(COBALT_input = NormalizedReadCount,
                                            COBALT_MIN_COV = COBALT_MIN_COV,
                                            COBALT_MAX_COV = COBALT_MAX_COV,
                                            CENTROMERES = CENTROMERES_g)
  #print(head(FilteredReadCount))
  print("# Filtered bins:")
  print(summary(factor(FilteredReadCount$FILTER)))

  print("# Binning ReadCounts")
  binnedReadCount_100kb <- .MergeBins(FilteredReadCount[FilteredReadCount$FILTER == "PASS",c("chromosome","position", "filteredReadCount")],
                               Binsize = 100000, Chr_sizes = CHR_SIZES)
  binnedReadCount_1mb <- .MergeBins(FilteredReadCount[FilteredReadCount$FILTER == "PASS",c("chromosome","position", "filteredReadCount")],
                               Binsize = 1000000, Chr_sizes = CHR_SIZES) # Used for Circos
  #print(head(binnedReadCount))

  print("# Segmenting filtered ReadCounts")
  ReadCount_Segments100kb <- .Segment_ReadCounts(ReadCounts = binnedReadCount_100kb[,c(1,2, 5)], Binsize = 100000, Cytoband = CYTOBAND, Chr_Sizes = CHR_SIZES)
  ReadCount_Segments1kb <- .Segment_ReadCounts(ReadCounts = binnedReadCount_100kb[,c(1,2, 5)], Binsize = 1000, Cytoband = CYTOBAND,  Chr_Sizes = CHR_SIZES)
  #print(head(ReadCount_Segments100kb))

  print("# Finemapping ReadCounts segments")
  #SegmentsCoarse, SegmentsFine, Cutoff_Low, Cutoff_High, Max_Finemapping_Dist = 200000
  ReadCount_Segments <- .FinemapSegments(SegmentsCoarse = ReadCount_Segments100kb[,c(2,4,5,7)],
                  SegmentsFine = ReadCount_Segments1kb[,c(2,4,5,7)],
                  Cutoff_Low = 1.5, Cutoff_High = 2.5, Max_Finemapping_Dist = 200000)

  # Export segments
  print(paste("# Writing: ", OUTPUT_dir, "readcounts.segments.raw.txt", sep = "."))
  write.table(ReadCount_Segments100kb, file = paste(OUTPUT_dir, "readcounts.segments.raw.txt", sep = "."), quote = F, row.names = F, sep = "\t")

  print(paste("# Writing: ", OUTPUT_dir, "readcounts.segments.txt", sep = "."))
  write.table(ReadCount_Segments, file = paste(OUTPUT_dir, "readcounts.segments.txt", sep = "."), quote = F, row.names = F, sep = "\t")
  # Export Readcounts
  print(paste("# Writing: ", OUTPUT_dir, "readcounts.filtered.1kb.txt", sep = "."))
  write.table(FilteredReadCount,
              file = paste(OUTPUT_dir, "readcounts.filtered.1kb.txt", sep = "."), quote = F, row.names = F, sep = "\t")

  print(paste("# Writing: ", OUTPUT_dir, "readcounts.filtered.100kb.txt", sep = "."))
  write.table(binnedReadCount_100kb,
              file = paste(OUTPUT_dir, "readcounts.filtered.100kb.txt", sep = "."), quote = F, row.names = F, sep = "\t")

  print(paste("# Writing: ", OUTPUT_dir, "readcounts.filtered.1mb.txt", sep = "."))
  write.table(binnedReadCount_1mb,
              file = paste(OUTPUT_dir, "readcounts.filtered.1mb.txt", sep = "."), quote = F, row.names = F, sep = "\t")

  ### bedpe files can be imported by the the StructuralVariantAnnotation package
  COBALT_bedpe <- .df_to_bedpe(ReadCount_Segments100kb, score = "mean", output_file = paste(OUTPUT_dir, "readcounts.segments.raw.bedpe", sep = "."))
  COBALT_finemapped_bedpe <- .df_to_bedpe(ReadCount_Segments, score = "mean", output_file = paste(OUTPUT_dir, "readcounts.segments.bedpe", sep = "."))
  print("### Finished ReadCount filtering")

  return(FilteredReadCount)
}

### Functions to merge BAF and CopyNumber data into a CNV table

Integrate_RD_BAF <- function(BAF, RD,
                             BAF_segments, RD_segments,
                             Type = "Loss"){

  # Select the segments of the right type (loss or gain)
  if(Type == "Loss" | Type == "LOH"){
    RD_seg <- RD_segments[which(RD_segments$QUAL < 1.5),]
    BAF_seg <- BAF_segments[which(BAF_segments$QUAL >= 0.3),]
  } else if (Type == "Gain") {
    RD_seg <- RD_segments[which(RD_segments$QUAL > 2.5),]
    BAF_seg <- BAF_segments[which(BAF_segments$QUAL > 0.16 & BAF_segments$QUAL < 0.3),]
  }

  # Overlap the ReadDepth and BAF segments of the same type (loss or gain)
  overlap_RD_BAF <- findBreakpointOverlaps(RD_seg, BAF_seg)

  # Select the regions that have both BAF and CN support for losses
  if(length(overlap_RD_BAF) > 0){
    if(Type != "LOH"){
      merged_CNVs <- intersect(RD_seg[queryHits(overlap_RD_BAF)],
                               BAF_seg[subjectHits(overlap_RD_BAF)])
    } else {
      # LOH regions only have "loss" at the BAF-level, but not at the ReadDepth (eg no overlap between LOH and loss regions)
      merged_CNVs <- BAF_seg[-subjectHits(overlap_RD_BAF),"QUAL"]
    }
    
    # If merged_CNVs is empty, then return an empty GRanges object.
    if (!length(merged_CNVs)){
      merged_CNVs2 <- GRanges()
      return(merged_CNVs2)
    }

    # Merge the losses if there's a gap of max 1e6 bp between two events
    merged_CNVs2 <- reduce(merged_CNVs, min.gapwidth = 1e6)
    merged_CNVs2$CopyNumber <- Type
    merged_CNVs2$RD <- NA
    merged_CNVs2$BAF <- NA

    # The coordinates of the CNVs can be different from the coordinates of the ReadDepth and/or BAF segments. Therefore recalculate the mean ReadDepth and BAF for each segment using the newly defined coordinates
    for(i in 1:length(merged_CNVs2)){
      #print(i)
      olap_CNVs_RD <- findOverlaps(merged_CNVs2[i], RD)
      merged_CNVs2[i]$RD <- mean(RD$RD[subjectHits(olap_CNVs_RD)], na.rm = T)

      olap_CNVs_BAF <- findOverlaps(merged_CNVs2[i], BAF)
      merged_CNVs2[i]$BAF <- mean(BAF$BAF[subjectHits(olap_CNVs_BAF)], na.rm = T)
    }
    if(Type == "LOH"){
      # Some (usually smaller) segments with LOH actually also have a low read depth, but they were not classified as Loss in the readdephth filtering. Correct for this:
      merged_CNVs2$CopyNumber[which(merged_CNVs2$RD < 1.5 & merged_CNVs2$BAF > 0.3)] <- "Loss"
      merged_CNVs2 <- merged_CNVs2[which(merged_CNVs2$BAF > 0.3),]
    }
  } else {
    merged_CNVs2 <- GRanges()
  }
  #print(head(merged_CNVs2))
  return(merged_CNVs2)
}

# Main function
Merge_CNVs <- function(BAF, RD,
                       BAF_segments, RD_segments, SVs, Output_dir, SAMPLE){
  print("# Integrating BAF and ReadCounts segments")
  losses <- Integrate_RD_BAF(BAF = BAF,
                             RD = RD,
                             BAF_segments = BAF_segments,
                             RD_segments = RD_segments, Type = "Loss")

  gains <- Integrate_RD_BAF(BAF = BAF,
                            RD = RD,
                            BAF_segments = BAF_segments,
                            RD_segments = RD_segments, Type = "Gain")

  loh <- Integrate_RD_BAF(BAF = BAF,
                          RD = RD,
                          BAF_segments = BAF_segments,
                          RD_segments = RD_segments,
                          Type = "LOH")
  CNVs <- c(losses, gains, loh)
  CNVs <- CNVs[order(seqnames(CNVs), start(CNVs)),]
  CNV_output <- as.data.frame(CNVs)
  names(CNV_output)[1] <- "chrom"

  # Overlap the CNVs with translocations to determine which translocations are true unbalanced translocations
  SVs_CTX <- Filter_Unbalanced_Translocations(CNVs = CNVs, SVs = SVs)
  SVs_CTX_Filtered <- SVs_CTX[which(rowRanges(SVs_CTX)$FILTER == "PASS"),]
  SVs_CTX_Filtered <- SVs_CTX_Filtered[as.vector(seqnames(rowRanges(SVs_CTX_Filtered))) %in% c(1:22, "X", "Y")]
  
  # Filter breakends of large CNVs if they do not overlap with the readcounts/BAF determined CNVs 
  SVs_Filtered <- Filter_Breakends_LargeCNVs(SVs = SVs_CTX_Filtered, CNVs = CNV_output)
  SVs_Filtered <- SVs_Filtered[which(rowRanges(SVs_Filtered)$FILTER == "PASS"),]
  
  print(summary(factor( VariantAnnotation::fixed(SVs_Filtered)$FILTER)))
  print(summary(factor(info(SVs_Filtered)$SVTYPE)))
  
  if(Output_dir != ""){
    print("# Writing output files")
    writeVcf(SVs_CTX_Filtered, filename = paste(Output_dir, ".integrated.svs.filtered.vcf", sep = ""))
    writeVcf(SVs_CTX, filename = paste(Output_dir, ".integrated.svs.annotated.vcf", sep = ""))
  }

  CNV_output_file <- paste(Output_dir, ".integrated.cnvs.txt", sep = "")
  print(paste("### Writing: ", CNV_output_file, sep = ""))
  write.table(CNV_output, file = CNV_output_file, quote = F, row.names = F, col.names = T, sep = "\t")

  return(CNV_output)
}

# Determine which CNVs have a breakend
# Filter unbalanced translocations

Filter_Unbalanced_Translocations <- function(CNVs, SVs){

  Translocations <- rowRanges(SVs[info(SVs)$SVTYPE == "CTX"])

  NearestNeighbour <- distanceToNearest(Translocations)
  Translocations$Dist <- NA
  Translocations$Dist[queryHits(NearestNeighbour)] <- mcols(NearestNeighbour)$distance
  # Assume that co-breakends of balanced translocations are within 1000000bp of eachother
  Unbalanced_Translocations <- SVs[names(Translocations)[which(Translocations$Dist > 1000000)], ]

  Breakends_CTX <- rowRanges(Unbalanced_Translocations)
  CNV_Starts <- GRanges(seqnames = seqnames(CNVs), IRanges(start(CNVs), start(CNVs)))
  CNV_Ends <- GRanges(seqnames = seqnames(CNVs), IRanges(end(CNVs), end(CNVs)))

  CNV_Breaks <- c(CNV_Starts, CNV_Ends)

  olap_Breakends_CNV <- findOverlaps(CNV_Breaks, Breakends_CTX, maxgap = 1e5)

  # Filter the seemingly unbalanced translocations not flanking a CNV (so not truly a translocation)
  if(length(olap_Breakends_CNV) > 0 ){
    Filtered_Translocations <- names(Breakends_CTX[-subjectHits(olap_Breakends_CNV),])
  } else {
    # Filter all CTX's if none overlaps with the CNVs
    Filtered_Translocations <- names(Breakends_CTX)
  }
  Filtered_Translocations2 <- unique(info(SVs)[Filtered_Translocations, "EVENT"])

  VariantAnnotation::fixed(SVs)$FILTER[which(unlist(info(SVs)$EVENT) %in% Filtered_Translocations2 & rowRanges(SVs)$FILTER != "" & rowRanges(SVs)$FILTER != "PASS")] <- paste(VariantAnnotation::fixed(SVs)$FILTER[which(unlist(info(SVs)$EVENT) %in% Filtered_Translocations2 & rowRanges(SVs)$FILTER != "" & rowRanges(SVs)$FILTER != "PASS")], "CTX_NO_MATE", sep = ";")
  VariantAnnotation::fixed(SVs)$FILTER[which(unlist(info(SVs)$EVENT) %in% Filtered_Translocations2 & rowRanges(SVs)$FILTER == "PASS")] <- "CTX_NO_MATE"
  VariantAnnotation::fixed(SVs)$FILTER[which(unlist(info(SVs)$EVENT) %in% Filtered_Translocations2 & rowRanges(SVs)$FILTER == "")] <- "CTX_NO_MATE"
  #print(summary(factor( VariantAnnotation::fixed(SVs)$FILTER)))
  return(SVs)
}

Filter_Breakends_LargeCNVs <- function(CNVs, SVs){
  
  
  CNV_Breakends <- SVs[info(SVs)$SVTYPE %in% c("DEL", "DUP")]
  CNV_Breakends_10mb <- CNV_Breakends[which(info(CNV_Breakends)$SVLEN > 10e6)]
  CNVs_g <- GRanges(seqnames = CNVs[,1], IRanges(start = CNVs[,2], end = CNVs[,3]))
  if(length(CNV_Breakends_10mb) > 0){
    olap_CNVs <- findOverlaps(rowRanges(CNV_Breakends_10mb), CNVs_g)
    if(length(olap_CNVs) > 0){
      Breakends_overlapping_CNVs <- names(CNV_Breakends_10mb[queryHits(olap_CNVs)])
      Filtered_Breakends <- CNV_Breakends_10mb[-Breakends_overlapping_CNVs,]
    } 
    else {
      Filtered_Breakends <- CNV_Breakends_10mb
    }
    Filtered_Events <- unique(info(SVs)[names(Filtered_Breakends), "EVENT"])
    
    VariantAnnotation::fixed(SVs)$FILTER[which(unlist(info(SVs)$EVENT) %in% Filtered_Events & rowRanges(SVs)$FILTER != "" & rowRanges(SVs)$FILTER != "PASS")] <- paste(VariantAnnotation::fixed(SVs)$FILTER[which(unlist(info(SVs)$EVENT) %in% Filtered_Events & rowRanges(SVs)$FILTER != "" & rowRanges(SVs)$FILTER != "PASS")], "NO_RD_SUPPORT", sep = ";")
    VariantAnnotation::fixed(SVs)$FILTER[which(unlist(info(SVs)$EVENT) %in% Filtered_Events & rowRanges(SVs)$FILTER == "PASS")] <- "NO_RD_SUPPORT"
    VariantAnnotation::fixed(SVs)$FILTER[which(unlist(info(SVs)$EVENT) %in% Filtered_Events & rowRanges(SVs)$FILTER == "")] <- "NO_RD_SUPPORT"
    
  } 

  return(SVs)
}
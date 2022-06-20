### This script can be used to filter a GRIPSS VCF against a blacklist, determine the BAF in the SV regions and determine the relative read depth in SV regions.

### Requirements
# StructuralVariantAnnotation R-package (can be installed in guix)

### Steps:
# 1. Filter SV VCF against blacklist
# 2. Overlap SV candidates with SNVs to determine BAF
# 3. Overlap SV candidates with Cobalt read depth data to determine coverage

### To do
# Add translocation filtering (eg balanced vs unbalanced)

### Limitations
# currently all INV < 1000bp are removed, most seem to be PTA artefacts

Filter_Breakends <- function(GRIPSS_PON,
                             SV_VCF_file,
                             BAFs,
                             ReadCounts,
                             Output_dir = "",
                             max_dist = 100,
                             max_SV_size_annotation = 1e6,
                             MaxBreakendCov = 100){

  print(paste("# Reading SV VCF: ", SV_VCF_file, sep = ""))
  SV_VCF <- readVcf(SV_VCF_file, "hg38")
  
  
  VariantAnnotation::fixed(SV_VCF)$FILTER[VariantAnnotation::fixed(SV_VCF)$FILTER == "PASS"] <- ""
  SV_VCF <- SV_VCF[as.vector(seqnames(rowRanges(SV_VCF))) %in% c(1:22, "X", "Y")]
  
  SVGR_raw <- breakpointRanges(SV_VCF) # Some SVs are not included here, because they don't have a partner breakend

  # Remove the breakends without a partner breakend
  #SV_VCF <- SV_VCF[which(names(SV_VCF) %in% names(SVGR_raw)),]
  SV_VCF <- SV_VCF[names(SVGR_raw),]
  
  # Rescue translocations with one filtered breakend from the unfiltered GRIPSS vcf
  GRIPSS_raw_file <- gsub(x = SV_VCF_file, pattern = ".filtered", replacement = "")
  ## Doesnt work at the moment in nextflow, needs to be added
  #SV_VCF_CTX <- .Rescue_Translocation_Breakends(SV_VCF, Unfiltered_Datafile = GRIPSS_raw_file)
  
  SV_VCF_CTX <- SV_VCF
  
  SVGR <- breakpointRanges(SV_VCF_CTX)
  
  # Add the SVTYPE to the VCF
  info(SV_VCF_CTX)$SVTYPE <- simpleEventType(SVGR[names(SV_VCF_CTX)])

  # Add the SVLEN to the VCF
  # header line needs to be added to the VCF:
  header_lines_SVLEN <- data.frame(Number = 1, Type = "Float", Description = "Size of the rearrangement (bp)")
  row.names(header_lines_SVLEN) <- "SVLEN"
  info(header(SV_VCF_CTX)) <- rbind(info(header(SV_VCF_CTX)), header_lines_SVLEN)
  
  info(SV_VCF_CTX)$SVLEN <- SVGR[names(SV_VCF_CTX)]$svLen

  SV_VCF_PON <- .Filter_PON(SV_VCF = SV_VCF_CTX, GRIPSS_PON = GRIPSS_PON, max_dist = max_dist)

  ### 2. Calculate the B-allele frequency of SNPs overlapping the candidate SVs
  simplegr <- SVGR[simpleEventType(SVGR) %in% c("INS", "INV", "DEL", "DUP")]
  simplebed <- data.frame(
    chrom=seqnames(simplegr),
    # call the centre of the homology/inexact interval
    start=as.integer((start(simplegr) + end(simplegr)) / 2),
    end=as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
    name=simpleEventType(simplegr),
    score=simplegr$QUAL,
    strand="."
  )
  # Only select each SV event once
  row.names(simplebed) <- names(simplegr)
  simplebed <- simplebed[simplebed$start < simplebed$end,]
  simplebed_g <- GRanges(simplebed)

  # Select only the SVs with a size < max_SV_size_annotation (large SVs take a long time to run)
  simplebed_g_small <- simplebed_g[width(simplebed_g) < max_SV_size_annotation,]
  simplebed_g_small <- simplebed_g_small[seqnames(simplebed_g_small) %in% c(1:22, "X", "Y"),]

  print("# Overlapping SV candidates with SNPs")
  SV_VCF_BAF <- .BAF_Annotate_Breakends(SVs = simplebed_g_small,
                                       SV_VCF = SV_VCF_PON,
                                       SNVs = BAFs)

  print("# Overlapping SV candidates with ReadCounts")
  SV_VCF_BAF_RD <- .ReadCounts_Annotate_Breakends(SV_VCF = SV_VCF_BAF, SVs = simplebed_g_small,
                                                 ReadCounts = ReadCounts[,c(1,2,2,6)])

  SV_VCF_PostFilter <- .Postfilter_Breakends(SV_VCF = SV_VCF_BAF_RD, SV_BED = simplebed_g, MaxBreakendCov = MaxBreakendCov)
  
  #
  VariantAnnotation::fixed(SV_VCF_PostFilter)$FILTER[VariantAnnotation::fixed(SV_VCF_PostFilter)$FILTER == ""] <- "PASS"
  SV_VCF_Filtered <- SV_VCF_PostFilter[which(rowRanges(SV_VCF_PostFilter)$FILTER %in% c("PASS", "RESCUE")),]
  SV_VCF_Filtered <- SV_VCF_Filtered[as.vector(seqnames(rowRanges(SV_VCF_Filtered))) %in% c(1:22, "X", "Y")]

  print(summary(factor(rowRanges(SV_VCF_PostFilter)$FILTER)))
  print(summary(factor(info(SV_VCF_Filtered)$SVTYPE)))

  # info(SV_VCF_Filtered)$SVLEN
  # info(SV_VCF_Filtered)$SVTYPE

  if(Output_dir != ""){
    print("# Writing output files")
    writeVcf(SV_VCF_Filtered, filename = paste(Output_dir, ".svs.filtered.vcf", sep = ""))
    writeVcf(SV_VCF_BAF_RD, filename = paste(Output_dir, ".svs.annotated.vcf", sep = ""))
    writeVcf(SV_VCF_PostFilter, filename = paste(Output_dir, ".svs.postfilter.vcf", sep = ""))
  }

  print("# Done")

  return(SV_VCF_PostFilter)
}


.Rescue_Translocation_Breakends <- function(Breakends, Unfiltered_Datafile, SearchWindow = 1e5){
  
  Breakends_output <- Breakends
  
  # Select the translocations (CTX)
  BreakendsGR <- breakpointRanges(Breakends) # Some SVs are not included here, because they don't have a partner breakend
  info(Breakends)$SVTYPE <- simpleEventType(BreakendsGR[names(Breakends)])
  Breakends_CTX <- Breakends[info(Breakends)$SVTYPE == "CTX"]
  if(length(Breakends_CTX) > 0){
    Breakends_CTX <- Breakends_CTX[rowRanges(Breakends_CTX)$FILTER != "PON"]
    
    # Read part of the VCF at locations surrounding the translocations:
    Translocs <- GRanges(seqnames = seqnames(rowRanges(Breakends_CTX)), IRanges(start = start(rowRanges(Breakends_CTX)), end = end(rowRanges(Breakends_CTX))))
    start(Translocs) <- start(Translocs) - SearchWindow
    end(Translocs) <- end(Translocs) + SearchWindow
    GRIPSS_raw_data <- readVcf(Unfiltered_Datafile, "hg38", param=Translocs)
    GRIPSS_GR <- breakpointRanges(GRIPSS_raw_data)
    GRIPSS_GR$svtype <- simpleEventType(GRIPSS_GR)
    GRIPSS_CTX <- GRIPSS_GR[GRIPSS_GR$svtype == "CTX",]
    
    # Select the translocation breakends not present in the GRIPSS vcf:
    GRIPSS_CTX2 <- GRIPSS_CTX[-which(names(GRIPSS_CTX) %in% names(Breakends_CTX))]
    GRIPSS_CTX2 <- GRIPSS_CTX2[seqnames(GRIPSS_CTX2) %in% c(1:22, "X", "Y")]
    
    if(length(GRIPSS_CTX2) > 0){
      print(paste("Rescued ", length(GRIPSS_CTX2)/2, " translocation breakpoint junctions", sep = ""))
      # print(GRIPSS_CTX2)
      rescued <- GRIPSS_raw_data[names(GRIPSS_CTX2),]
      VariantAnnotation::fixed(GRIPSS_raw_data[names(GRIPSS_CTX2),])$FILTER <- "RESCUE"
      
      Breakends_output <- rbind(Breakends_output, GRIPSS_raw_data[names(GRIPSS_CTX2),])
    }
    
  }
  print(paste("# ", length(Breakends_output) - length(Breakends), " rescued translocation breakends",sep = ""))
  #Breakends_output <- Breakends_output[order(seqnames(rowRanges(Breakends_output)), start(rowRanges(Breakends_output))),]
  return(Breakends_output)
}

.Filter_PON <- function(SV_VCF, GRIPSS_PON, max_dist = 100){

  PON_g <- GRanges(seqnames = GRIPSS_PON[,1], IRanges(start = GRIPSS_PON[,2], end = GRIPSS_PON[,3]))
  #print(head(PON_g))

  print(paste("# Number of raw breakends: ", length(SV_VCF), sep = ""))
  olap_Breakends_PON <- findOverlaps(SV_VCF, PON_g, maxgap = max_dist)

  print(paste("# Number of variants overlapping with blacklist: ", length(olap_Breakends_PON), sep = ""))

  Filtered_Breakends <- names(SV_VCF)[queryHits(olap_Breakends_PON)]
  Filtered_Events <- info(SV_VCF)[Filtered_Breakends, "EVENT"]

  VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER != "")] <- paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER != "")], "BLACKLIST_PARTNER", sep = ";")
  VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER == "")] <- "BLACKLIST_PARTNER"

  VariantAnnotation::fixed(SV_VCF)$FILTER[which(names(SV_VCF) %in% Filtered_Breakends & rowRanges(SV_VCF)$FILTER != "")] <- paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(names(SV_VCF) %in% Filtered_Breakends & rowRanges(SV_VCF)$FILTER != "")], "BLACKLIST", sep = ";")
  VariantAnnotation::fixed(SV_VCF)$FILTER[which(names(SV_VCF) %in% Filtered_Breakends & rowRanges(SV_VCF)$FILTER == "")] <- "BLACKLIST"

  print(summary(factor(rowRanges(SV_VCF)$FILTER)))
  return(SV_VCF)
}


.BAF_Annotate_Breakends <- function(SVs, SV_VCF, SNVs){

  SNVs_g <- GRanges(seqnames = SNVs[,1], IRanges(start = SNVs[,2], end = SNVs[,2]), VAF_normalized = SNVs$VAF_normalized)
  Filtered_Events <- c()

  # These lines need to be added to the output VCF header
  header_lines_BAF <- data.frame(Number = c(rep(1, 4)), Type = rep("Float", 4), Description = c("Number of SNPs heterozygous in control overlapping the SV",
                                                                                                "Mean VAF of the SNPs overlapping the SV",
                                                                                                "Median VAF of the SNPs overlapping the SV",
                                                                                                "Standard deviation of the mean VAF of  SNPs overlapping the SV"))
  row.names(header_lines_BAF) <- c("het_SNPs", "VAF_SNP_mean", "VAF_SNP_median", "VAF_SNP_sd")
  info(header(SV_VCF)) <- rbind(info(header(SV_VCF)), header_lines_BAF)

  info(SV_VCF)$het_SNPs <- 0
  info(SV_VCF)$VAF_SNP_mean <- NA
  info(SV_VCF)$VAF_SNP_median <- NA
  info(SV_VCF)$VAF_SNP_sd <- NA

  for(i in 1:length(SVs)){
    #print(i)
    SV <- SVs[i,]
    event <- info(SV_VCF)[names(SV), "EVENT"]
    if(SV$name %in% c("DEL", "DUP", "INV")){
      #print(SV)
      VAFs <- c()
      olap_SV_SNV <- findOverlaps(SV, SNVs_g)
      if(length(olap_SV_SNV) > 0){
        VAFs <- SNVs_g$VAF_normalized[subjectHits(olap_SV_SNV)]
        #
        # snvs_SV <- snv_vcf_het[subjectHits(olap_SV_SNV),]
        #
        # DPs <- geno(snvs_SV)$DP[names(snvs_SV), sample_name]
        # ADs <- geno(snvs_SV)$AD[names(snvs_SV), sample_name]
        # ADs2 <- ADs[which(DPs > min_SNP_depth)]
        # # Calculate for each SNP how much the VAF is divergent from 0.5 (eg if the VAF is 0.9, it's 0.4 from 0.5)
        # VAFs <- unlist(lapply(ADs2, function(x) abs(0.5-x[2]/(x[1]+x[2]))))
        #
        if(length(VAFs) > 0){
          info(SV_VCF)[info(SV_VCF)$EVENT == event, "het_SNPs"] <- length(VAFs)
          info(SV_VCF)[info(SV_VCF)$EVENT == event, "VAF_SNP_mean"] <- mean(VAFs[!is.na(VAFs)])
          info(SV_VCF)[info(SV_VCF)$EVENT == event, "VAF_SNP_median"] <- median(VAFs[!is.na(VAFs)])
          info(SV_VCF)[info(SV_VCF)$EVENT == event, "VAF_SNP_sd"] <- sd(VAFs[!is.na(VAFs)])
        }
        if(SV$name == "DUP" & mean(VAFs[!is.na(VAFs)]) < 0.18 & length(VAFs) > 1){
          Filtered_Events <- c(Filtered_Events, event)
        }
        if(SV$name == "DEL" & mean(VAFs[!is.na(VAFs)]) < 0.4 & length(VAFs) > 1){
          Filtered_Events <- c(Filtered_Events, event)
        }
      }
    }
  }
  print(paste("# Breakends without BAF support: ", length(Filtered_Events), sep = ""))
  if(length(Filtered_Events) > 0 ){
    VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER != "")] <- paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER != "")], "NO_BAF_SUPPORT", sep = ";")
    VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER == "")] <- "NO_BAF_SUPPORT"
  }
  return(SV_VCF)
}


.ReadCounts_Annotate_Breakends <- function(SVs, SV_VCF, ReadCounts){

  ReadCounts_g <- GRanges(seqnames = ReadCounts[,1], IRanges(start = ReadCounts[,2], end = ReadCounts[,3]+999), RD = ReadCounts[,4])
  Filtered_Events <- c()

  # header line needs to be added to the VCF:
  header_lines_readdepth <- data.frame(Number = 1, Type = "Float", Description = "Read depth in sample relative to control samples")
  row.names(header_lines_readdepth) <- "Rel_read_depth"
  info(header(SV_VCF)) <- rbind(info(header(SV_VCF)), header_lines_readdepth)

  info(SV_VCF)$Rel_read_depth <- NA

  for(i in 1:length(SVs)){
    #print(i)
    SV <- SVs[i,]
    if(SV$name %in% c("DEL", "DUP", "INV")){
      #print(SV)
      olap_SV_ReadCounts <- findOverlaps(SV, ReadCounts_g, minoverlap = 500)
      if(length(olap_SV_ReadCounts) > 0){

        EVENT <- info(SV_VCF)[names(SV), "EVENT"]
        regions <- ReadCounts[subjectHits(olap_SV_ReadCounts),]

        ### Add read depth info to the VCF entry of both breakends
        info(SV_VCF)[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"]),"Rel_read_depth"] <- mean(regions[,4])

        if(SV$name == "DEL" & unique(info(SV_VCF)[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"]),"Rel_read_depth"]) > 1.5){
          Filtered_Events <- c(Filtered_Events, EVENT)
        } else if (SV$name == "DUP" & unique(info(SV_VCF)[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"]),"Rel_read_depth"]) < 2.5){
          Filtered_Events <- c(Filtered_Events, EVENT)
        }
      }
    }
  }
  print(paste("# Breakends without ReadDepth support: ", length(Filtered_Events), sep = ""))
  if(length(Filtered_Events) > 0 ){
    VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER != "")] <- paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER != "")], "NO_RD_SUPPORT", sep = ";")
    VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Events & rowRanges(SV_VCF)$FILTER == "")] <- "NO_RD_SUPPORT"
  }
  return(SV_VCF)
}

#
#   print(SV)
#   print(info(SV_VCF)[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"]),"Rel_read_depth"])
#
#   info(SV_VCF)[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"]),"Rel_read_depth"]
#
#   if(VariantAnnotation::fixed(SV_VCF)$FILTER[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"])] == ""){
#     VariantAnnotation::fixed(SV_VCF)$FILTER[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"])] <- "NO_RD_SUPPORT"
#   } else {
#     VariantAnnotation::fixed(SV_VCF)$FILTER[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"])] <-
#       paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"])], "NO_RD_SUPPORT", sep = ";")
#   }
# }

# if(SV$name == "DUP" &  info(SV_VCF)[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"]),"Rel_read_depth"] < 2.5){
#   print(SV)
#   print(info(SV_VCF)[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"]),"Rel_read_depth"])
#
#   info(SV_VCF)[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"]),"Rel_read_depth"]
#
#   if(VariantAnnotation::fixed(SV_VCF)$FILTER[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"])] == ""){
#     VariantAnnotation::fixed(SV_VCF)$FILTER[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"])] <- "NO_RD_SUPPORT"
#   } else {
#     VariantAnnotation::fixed(SV_VCF)$FILTER[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"])] <-
#       paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(info(SV_VCF)$EVENT == info(SV_VCF)[names(SV), "EVENT"])], "NO_RD_SUPPORT", sep = ";")
#   }
#
# }

.Postfilter_Breakends <- function(SV_VCF, SV_BED, MaxBreakendCov = 100){
  # Flag the Breakends with excessive coverage
  Excess_COV_Events <- unique(info(SV_VCF)[which(info(SV_VCF)$REF > MaxBreakendCov), "EVENT"])

  VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Excess_COV_Events & rowRanges(SV_VCF)$FILTER != "")] <- paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Excess_COV_Events & rowRanges(SV_VCF)$FILTER != "")], "EXCESS_COV", sep = ";")
  VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Excess_COV_Events & rowRanges(SV_VCF)$FILTER == "")] <- "EXCESS_COV"

  # Many seeming inversions are called in PTA data. True inversions should have both ++ and -- reads. However, in the PTA data these candidates usually often only have ++ or -- and not both. Possible small circular amplifications that are formed during the PTA amplification.
  Inversions <- SV_BED[SV_BED$name == "INV",]
  if(length(Inversions) > 0){
    Inversions$Dist <- 1e6
    NearestNeighbour <- distanceToNearest(Inversions)
    Inversions$Dist[queryHits(NearestNeighbour)] <- mcols(NearestNeighbour)$distance
    # Assume that co-breakends are within 100bp of eachother
    Filtered_Inversions <- info(SV_VCF)[names(Inversions)[which(Inversions$Dist > 100)], "EVENT"]
    
    VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Inversions & rowRanges(SV_VCF)$FILTER != "")] <- paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Inversions & rowRanges(SV_VCF)$FILTER != "")], "INV_NO_MATE", sep = ";")
    VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Inversions & rowRanges(SV_VCF)$FILTER == "")] <- "INV_NO_MATE"
  } else {
    
  }
  

  # Filter all INV smaller than 1kb (most appear to be PTA artefacts)
  Filtered_Inversions_Size <- info(SV_VCF)[which(info(SV_VCF)$SVLEN < 1000 & info(SV_VCF)$SVTYPE == "INV"), "EVENT"]
  VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Inversions_Size & rowRanges(SV_VCF)$FILTER != "")] <- paste(VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Inversions_Size & rowRanges(SV_VCF)$FILTER != "")], "INV_SIZE", sep = ";")
  VariantAnnotation::fixed(SV_VCF)$FILTER[which(unlist(info(SV_VCF)$EVENT) %in% Filtered_Inversions_Size & rowRanges(SV_VCF)$FILTER == "")] <- "INV_SIZE"

  # Similar, CTXs balanced translocations should have breakends (breakpoint junctions) on both derivative chromosomes. Unbalanced translocations can have just one breakpoint junction, but there should be a gain or loss next to them. Check if this true

  #print(summary( factor(VariantAnnotation::fixed(SV_VCF)$FILTER)))
  return(SV_VCF)
}

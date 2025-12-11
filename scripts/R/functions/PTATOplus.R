# This script contains functions to analyze PTATO SNV/Indel data

library(dplyr)
library(ggplot2)
library(reshape2)
# This function selects the PTAprobCutoffs from the header of a FILTERED PTATO VCF
getPTAprobCutoff <- function(VCF_file){
  VCF_header <- scanVcfHeader(VCF_file)
  PTAprobCutoff_raw <- meta(VCF_header)$PTAprobCutoff[1,1]
  PTAprobCutoff_A  <- gsub(pattern = "\"", replacement = "", x = PTAprobCutoff_raw, fixed=TRUE)
  PTAprobCutoff_B <- gsub(pattern = "\\(|\\)|- ", replacement = "", x = PTAprobCutoff_A)
  PTAprobCutoffs <- unlist(strsplit(PTAprobCutoff_B, split = " "))
  names(PTAprobCutoffs) <- c("PTAprobCutoff", "Walker", "Cossim")
  return(PTAprobCutoffs)
}


# This function can be used to transform a VariantAnnotation vcf to a MutationalPatterns grangeslist (human only)
VCF_to_GR <- function(VCF, chromosomes = c(1:22, "X", "Y"), add_chr = TRUE){
  GR <- rowRanges(VCF)
  GR <- GR[which(seqnames(GR) %in% chromosomes),]
  seqlevels(GR) <- seqlevels(GR)[1:24]
  # Mutationalpatterns usually uses chr
  if(add_chr == TRUE){
    seqlevels(GR) <- paste("chr", seqlevels(GR), sep = "")
  } 
  return(GR)
}

# Variants with a VAF<VAF_threshold will be flagged with FAIL_VAF in the FILTER field.
# Variants with PTAprob>PTAprobCutoff will be flagged with FAIL in the FILTER field.
readPTATOvcf <- function(vcf_unfiltered, 
                         vcf_filtered = "",
                         cutoff = 1, # 1 = mean, 2 = Lira, 3 = Cossim
                         type = "snv",
                         VAF_threshold = 0.25,
                         chromosomes = c(1:22, "X", "Y"),
                         genome = "hg38"){
  
  vcf <- readVcf(vcf_unfiltered, genome = genome)
  vcf <- vcf[as.vector(seqnames(vcf)) %in% chromosomes,]
  GR <- VCF_to_GR(vcf, chromosomes = chromosomes)
  
  if(type != "indel"){
    PTAprobCutoff <- getPTAprobCutoff(vcf_filtered)
  } else {
    # indels that are removed have a PTAprob of 1
    PTAprobCutoff <- 0.5
  }
  
  GR$VAF <- as.numeric(geno(vcf)$VAF)
  # For each variant, PTATO reports 3 probabilities in the geno field. The 1st is calculated using all random forest features. For some variants, not all features can be determined. For these variants, a random forest without XX or XX is run, for which the probabilities are in the 2nd and 3rd column respectively.
  PTAprobs <- as.data.frame(geno(vcf)$PTAprob)
  # Replace the probs that could not be determine by the second RF by the probs of the 3rd RF
  PTAprobs[which(is.na(PTAprobs[,2])),2] <- PTAprobs[which(is.na(PTAprobs[,2])),3]
  # Replace the probs that could not be determine by the first RF by the probs of the 2nd RF
  PTAprobs[which(is.na(PTAprobs[,1])),1] <- PTAprobs[which(is.na(PTAprobs[,1])),2]
  
  GR$PTAprob <- as.numeric(PTAprobs[,1])

  GR$PTAprobCutoff <- PTAprobCutoff[1]
  GR$PTAprobCutoff_Walker <- PTAprobCutoff[2]
  GR$PTAprobCutoff_Cossim <- PTAprobCutoff[3]
  
  GR$FILTER[GR$VAF < VAF_threshold] <- "FAIL_VAF"
  GR$FILTER[GR$PTAprob > PTAprobCutoff[cutoff]] <- "FAIL"
  
  return(GR)
}

# Read all SNV VCFs from one PTATO run
readPTATO_SNVs <- function(PTATO_dir, MinimalVAF = 0.2){
  
  SBSs_raw <- list()
  
  vcfs_unfiltered <- list.files(PTATO_dir,
                               pattern = ".snvs.ptato.vcf.gz$", full.names = T, recursive = T)
  #print(vcfs_unfiltered)
  Samples <- gsub(x = gsub(pattern = ".*_", replacement = "", vcfs_unfiltered), pattern = ".snvs.ptato.vcf.gz", replacement = "")
  
  for(Sample in Samples){
    print(Sample)
    vcf_unfiltered <- list.files(PTATO_dir,
                                 pattern = paste(Sample, ".snvs.ptato.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(vcf_unfiltered)
    vcf_filtered <- list.files(PTATO_dir,
                               pattern = paste(Sample, ".snvs.ptato.filtered.vcf.gz$", sep = ""), full.names = T, recursive = T)
    print(vcf_filtered)
    SBSs_raw[[Sample]] <- readPTATOvcf(vcf_unfiltered = vcf_unfiltered, 
                                       vcf_filtered = vcf_filtered, 
                                       VAF_threshold = MinimalVAF)
  }
  
  return(SBSs_raw)
  
}

# This plots the number of variants on the Y-axis and the PTAprobs on the X-axis
plot_PTAprobs_cumulative <- function(GR, Labels = "", Cutoff_label = TRUE, Colors = NULL){
  Probs <- data.frame()
  for(sample in names(GR)){
    print(sample)
    if(length(Labels) == length(names(GR))){
      Label <- Labels[which(sample == names(GR))]

    } else {
      Label <- sample
    }
    probs_sample <- data.frame(Sample = sample, 
                               Label = Label,
                               PTAprob = GR[[sample]]$PTAprob, 
                               PTAprobCutoff = GR[[sample]]$PTAprobCutoff, 
                               TotalVariants = length(GR[[sample]]))
    Probs <- rbind(Probs, probs_sample)
  }
  
  if(Cutoff_label == TRUE){
    Probs$Label <- paste(Probs$Label, " (", Probs$PTAprobCutoff, ")", sep = "")
  }
  
  if(!is.null(Colors)){
    names(Colors) <- unique(Probs$Label)
  }

  Probs$Label <- factor(Probs$Label, levels = unique(Probs$Label[order(Probs$TotalVariants , decreasing = T)]))
  Probs$Sample <- factor(Probs$Sample, levels = names(GR))

  plot <- ggplot(Probs, aes(x = PTAprob, col = Label)) +
    geom_step(aes(len=TotalVariants,y=..y.. * len), stat="ecdf", linewidth = 0.5) + 
    geom_vline(aes(xintercept = as.numeric(PTAprobCutoff), col = Label), linetype = 2, alpha = 0.7) +
    scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + 
    theme_classic(base_size = 8) + theme(panel.grid.major = element_line()) +
    labs(y = "Somatic variants (#)")
  if(!is.null(Colors)){
    plot <- plot + scale_color_manual(values = Colors)
  }
  return(plot)
}

# This plots the distribution of all PTAprobs in each sample
plot_PTAprobs_distribution <- function(GRL, Labels = "", Cutoff_label = TRUE, Colors = NULL){
  Probs <- data.frame()
  for(sample in names(GRL)){
    print(sample)
    if(length(Labels) == length(names(GRL))){
      Label <- Labels[which(sample == names(GRL))]
    } else {
      Label <- sample
    }
    probs_sample <- data.frame(Sample = sample, 
                               Label = Label,
                               PTAprob = GRL[[sample]]$PTAprob, 
                               PTAprobCutoff = GRL[[sample]]$PTAprobCutoff, 
                               TotalVariants = length(GRL[[sample]]))
    Probs <- rbind(Probs, probs_sample)
  }
  Probs$Sample <- factor(Probs$Sample, levels = names(GRL))
  
  if(Cutoff_label == TRUE){
    Probs$Label <- paste(Probs$Label, " (", Probs$PTAprobCutoff, ")", sep = "")
  }
  if(!is.null(Colors)){
    names(Colors) <- unique(Probs$Label)
  }
  Probs$Label <- factor(Probs$Label, levels = unique(Probs$Label[order(Probs$TotalVariants , decreasing = T)]))
  
  plot <- ggplot(Probs, aes(x = PTAprob, color = Label)) + 
    geom_density(linewidth = 0.5) + scale_y_continuous(expand = c(0, 0))  + 
    geom_vline(aes(xintercept = as.numeric(PTAprobCutoff), col = Label), linetype = 2) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic(base_size = 8) + theme(panel.grid.major = element_line())
  
  if(!is.null(Colors)){
    plot <- plot + scale_color_manual(values = Colors)
  }
  
  return(plot)
}


# Plot the number of variants passing or failing the FILTER
plot_PTAprobs_filter <- function(GRL, mode = "absolute", Labels = ""){
  Probs <- data.frame()
  for(sample in names(GRL)){
    print(sample)
    if(length(Labels) == length(names(GRL))){
      Label <- Labels[which(sample == names(GRL))]
    } else {
      Label <- sample
    }
    probs_sample <- data.frame(Sample = sample, 
                               Label = Label,
                               PTAprob = GRL[[sample]]$PTAprob, 
                               PTAprobCutoff = GRL[[sample]]$PTAprobCutoff, 
                               TotalVariants = length(GRL[[sample]]),
                               Filter = GRL[[sample]]$FILTER)
    Probs <- rbind(Probs, probs_sample)
  }
  
  Probs_freq <- Probs %>% group_by(Label, Filter) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))

  #print(Probs_freq)
  
  #Probs_freq$Sample <- factor(Probs_freq$Sample, levels = names(GRL))
  
  if(mode == "absolute"){
    Probs_freq$Text <- Probs_freq$n
    Probs_freq$Text[Probs_freq$n < max(Probs_freq$n)*0.05] <- ""
    
    plot <- ggplot(Probs_freq, aes(fill = Filter, x = Label, y = n)) + 
      geom_bar(stat = "identity", width = 0.8) + 
      geom_text(data = Probs_freq, aes(label = Text, y=n),
                position = position_stack(vjust = 0.5), size = 2.5) +
      scale_fill_manual(values = c("FAIL" = "#D55E00", "FAIL_VAF" = "#669bbc", "PASS"= "#009e73")) +
      #scale_fill_brewer(palette = "Set1") + 
      scale_y_continuous(expand = c(0, 0)) + 
      labs(y = "Somatic variants", x = "Sample") + 
      theme_classic(base_size = 8) + 
      theme(axis.text.x = element_text(angle = 60,  hjust=1))
  } else if (mode == "relative"){
    Probs_freq$Text <- paste(round(Probs_freq$freq*100,0), "%", sep ="")
    Probs_freq$Text[Probs_freq$freq < max(Probs_freq$freq)*0.1] <- ""
    
    plot <- ggplot(Probs_freq, aes(fill = Filter, x = Label, y = freq)) + 
      geom_bar(stat = "identity", width = 0.8) + 
      geom_text(data = Probs_freq, aes(label = Text, y=freq),
                position = position_stack(vjust = 0.5), size = 2.5) +
      theme_classic(base_size = 8) + 
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = c("FAIL" = "#D55E00", "FAIL_VAF" = "#669bbc", "PASS"= "#009e73")) +
      #scale_fill_brewer(palette = "Set1") + 
      theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
      labs(y = "Somatic variants (% of total)", x = "Sample")
  }
  return(plot)
}

# Plot the amount of CallableLoci basepairs for each sample in Samples
plotCallableLoci <- function(CallableLoci_Stats_Dir, Samples, plot = TRUE){
  stats_files <- list.files(CallableLoci_Stats_Dir, pattern = ".txt$", full.names = T)
  inputs <- list()
  stats <- data.frame()
  for(sample in Samples){
    inputs[[sample]] <- stats_files[grep(sample, stats_files)]
  }
  
  for(sample in names(inputs)){
    if(length(inputs[[sample]]) > 0){
      stats_sample <- read.delim(inputs[[sample]])
      stats_sample$Sample <- sample
      stats <- rbind(stats, stats_sample)
    } else {
      print(paste("! No CallableLoci.stats file found for sample: ",sample, sep = ""))
    }

  }
  stats <- stats[stats$State != "REF_N",]
  stats$State <- factor(stats$State, levels = c("EXCESSIVE_COVERAGE","LOW_COVERAGE","NO_COVERAGE","POOR_MAPPING_QUALITY", "CALLABLE"))
  colors <- c("EXCESSIVE_COVERAGE" = "#A65628","LOW_COVERAGE" = "#FF7F00","NO_COVERAGE" = "#E41A1C","POOR_MAPPING_QUALITY" = "#FFFF33", "CALLABLE" = "#377EB8")
  stats$Label <- paste(round(stats$Freq*100,1),  "%", sep ="")
  stats$Label[stats$Freq < 0.1] <- ""
  
  if(plot == TRUE){
    stats_plot <- ggplot(stats, aes(x = Sample, y = Freq, fill = State)) + geom_bar(stat = "identity") +
      geom_text(data = stats, aes(label = Label, y=Freq),
                position = position_stack(vjust = 0.5), size = 2.5) +
      scale_fill_manual(values = colors) + 
      theme_classic(base_size = 8) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 60,  hjust=1))
    return(stats_plot)
  } else {
    return(stats)
  }
 
}

# Plot the precision-recall curves based on the LiRA 
# It will search for the ".ptatotable.txt" file within the inputdir, so best to make the inputdir a bit specific (eg "inputdir = [PROJECT_DIR]/[individual]/PTATO/snvs/" )
plot_precision_recall <- function(inputdir, sample, PTAprobsCutoff = 0){
  ptatotable_files <- list.files(path = inputdir, pattern = ".ptatotable.txt", recursive = T, full.names = T)
  #ptatotable_file <- ptatotable_files[grep(sample, ptatotable_files)]
  ptatotable_file <- ptatotable_files[grep(paste0(sample, ".ptatotable.txt"), ptatotable_files)]
  print(ptatotable_file)
  if(file.exists(ptatotable_file) == T & length(ptatotable_file) == 1){
    print(ptatotable_file)
    ptatotable <- read.table(ptatotable_file, header = T)
    precision <- melt(ptatotable[,c("PTAprob_cutoff", "precision", "recall")], id.vars = "PTAprob_cutoff")
    
    plot <- ggplot(precision, aes(x = PTAprob_cutoff, y = value, col = variable)) + 
      geom_line() + 
      geom_vline(xintercept = PTAprobsCutoff[i], linetype= 2) +
      theme_classic(base_size = 8) + 
      scale_y_continuous(expand = c(0, 0.01)) + 
      scale_x_continuous(expand = c(0, 0)) + 
      ggtitle(sample) +
      theme(panel.grid.major = element_line(color = "lightgray",
                                            linetype = 2), plot.title = element_text(hjust = 0.5)) +
      labs(y = "Value", col = "Variable")
    
    return(plot)
  } else {
    print(paste("PTATO table not found: ", sample, sep =""))
  }
}

# Plot the precision-recall curves based on the LiRA for multiple samples (within one PTATO folder)
# It will search for the ".ptatotable.txt" file within the inputdir, so best to make the inputdir a bit specific (eg "inputdir = [PROJECT_DIR]/[individual]/PTATO/snvs/" )
plot_precision_recall_multisample <- function(ptatotable_files, samples, PTAprobsCutoffs = 0){
  #ptatotable_files <- list.files(path = paste(inputdir, "snvs/", sep = ""), pattern = ".ptatotable.txt", recursive = T, full.names = T)
  print(ptatotable_files)
  precision_merged <- data.frame()
  i <- 1
  for(sample in samples){
    print(sample)
    #ptatotable_file <- ptatotable_files[grep(sample, ptatotable_files)]
    ptatotable_file <- ptatotable_files[grep(paste0(sample, ".ptatotable.txt"), ptatotable_files)]
    print(ptatotable_file)
    if(file.exists(ptatotable_file) == T & length(ptatotable_file) == 1){
      #print(ptatotable_file)
      ptatotable <- read.table(ptatotable_file, header = T)
      precision <- melt(ptatotable[,c("PTAprob_cutoff", "precision", "recall")], id.vars = "PTAprob_cutoff")
      precision$Sample <- sample
      # Recalculate walker cutoff if they are not supplied in arguments
      if(PTAprobsCutoffs == 0){
        ptatotable2 <- ptatotable[which( abs(ptatotable$recall-ptatotable$precision) > 0 ),]
        ptatotable2$prec_recall <- abs(ptatotable2$precision-ptatotable2$recall)
        lowest <- ptatotable2[which(ptatotable2$prec_recall == min(ptatotable2$prec_recall, na.rm = T)),]
        if(nrow(lowest) > 1){
          lowest <- lowest[1,]
        }
        precision$Cutoff <- lowest$PTAprob_cutoff
      } else {
        precision$Cutoff <- PTAprobsCutoffs[i]
        
      }
      #print(head(precision))
      i <- i+1
      precision_merged <- rbind(precision_merged, precision)
    } else {
      print(paste("PTATO table not found: ", sample, sep =""))
    }
  }
  #print(head(precision_merged))
  plot <- ggplot(precision_merged, aes(x = PTAprob_cutoff, y = value, col = variable)) + 
    geom_line(linewidth = 1) + 
    geom_vline(aes(xintercept = Cutoff), linetype= 2) +
    facet_wrap(~Sample, ncol = 3) +
    theme_classic(base_size = 8) + 
    scale_y_continuous(expand = c(0, 0.01)) + 
    scale_x_continuous(expand = c(0, 0)) + 
    theme(panel.grid.major = element_line(color = "grey85",
                                          linetype = 2, size = 0.2), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, size = 0.75)) +
    labs(y = "Value", col = "Variable")
  
  return(plot)
}


# Plot the number of mutations that pass or fail the PTAprobsCutoff at different Cutoffs
plot_mutations_cutoff <- function(mut_mat_cosim){
  Variants_PASS <- data.frame(Cutoff = as.numeric(colnames(mut_mat_cosim[["PASS"]])),
                              Variants  = colSums(mut_mat_cosim[["PASS"]]), 
                              Type = "PASS")
  Variants_FAIL <- data.frame(Cutoff = as.numeric(colnames(mut_mat_cosim[["FAIL"]])),
                              Variants  = colSums(mut_mat_cosim[["FAIL"]]), 
                              Type = "FAIL")
  Variants <- rbind(Variants_PASS, Variants_FAIL)
  
  plot <- ggplot(Variants, aes(x = Cutoff, y = Variants, col = Type)) + geom_line() +
    geom_vline(xintercept = mut_mat_cosim[["PTAprobsCutoff"]]) +
    geom_vline(xintercept = mut_mat_cosim[["CosSim_Cutoff"]], linetype = 2)
  
  return(plot)
}


### Walker / Linked read functions
# The plots can be generated separately, or in one report file: makeWalkerReport(walker_vcf_directory, individual, samples, output_dir)

# ! Sizes of plots are not yet adjusted for the number of input samples

# Standard set of colors used 
walker_colors <- c("#D55E00", "#669bbc", "#009e73", "gray")
names(walker_colors) <- c("Low", "Intermediate", "High", NA)

# This function reads the Walker output vcf.gz files in the "walker_vcf_directory". "samples" should be in the vcf.gz filenames.
# Can only run for snv or indel seperately
# Output is a dataframe with all the variants with their walker scores
# Thresholds are now determine arbitarily
RetrieveWalkerScores <- function(walker_vcf_files, samples, type = "snv",
                                 high_score_threshold = 1000, 
                                 low_score_threshold = 1){
  # walker_vcf_directory
  #walker_vcf_files <- list.files(walker_vcf_directory, pattern = ".vcf.gz$", full.names = T)
  walker_scores <- data.frame()
  for(sample in samples){
    print(sample)
    walker_vcf_file <- walker_vcf_files[grep(paste0(sample, ".walker.vcf.gz"), walker_vcf_files)]
    print(walker_vcf_file)
    # Read the Walker VCF
    walker_vcf <- readVcf(walker_vcf_file, genome = "hg38")
    # Select only the snv or the indel
    if(type == "snv"){
      walker_variants <- walker_vcf[isSNV(walker_vcf),]
    } else if (type == "indel"){
      walker_variants <- walker_vcf[isIndel(walker_vcf),]
    }
    # Select autosomes only
    walker_variants <- walker_variants[as.vector(seqnames(walker_variants)) %in% c(1:22),]
    
    # Collect all the walker scores
    walker_scores_sample <- data.frame(Variant = row.names(geno(walker_variants)$WS), 
                                       WS = as.numeric(geno(walker_variants)$WS), check.names = F)
    walker_scores_sample$WS_Count <- walker_scores_sample$WS
    walker_scores_sample$WS_Count[which(as.numeric(walker_scores_sample$WS) < low_score_threshold)] <- "Low"
    walker_scores_sample$WS_Count[which(as.numeric(walker_scores_sample$WS) > high_score_threshold)] <- "High"
    walker_scores_sample$WS_Count[which(as.numeric(walker_scores_sample$WS) >= low_score_threshold  &as.numeric(walker_scores_sample$WS) <= high_score_threshold)] <- "Intermediate"
    walker_scores_sample$WS_Count[is.na(walker_scores_sample$WS_Count)] <- NA
    walker_scores_sample$WS_Count <- factor(walker_scores_sample$WS_Count, levels = c(NA, "Low", "Intermediate", "High"))
    walker_scores_sample$sample <- sample
    walker_scores <- rbind(walker_scores, walker_scores_sample)
  }
  return(walker_scores)
}

# This calculaes the number of variants for each category (Low, Intermediate, High, NA) per sample
CalculateWalkerCounts <- function(walker_scores){
  walker_counts <- walker_scores %>% group_by(sample, WS_Count) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  return(walker_counts)
}


# The Walker scores are extrapolated by multiplying the amount of undetermined variants (WS_Count=NA) by the frequency of the called variants (WS_Count!=NA)
Extrapolate_Walker <- function(walker_counts, walker_scores){
  walker_counts_extrapolated <- data.frame()
  walker_counts_called <- walker_scores[!is.na(walker_scores$WS_Count),] %>% group_by(sample, WS_Count) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  ## Extrapolate walker
  for(sample in unique(walker_counts$sample)){
    #print(sample)
    walker_counts_uncalled_sample <- walker_counts[walker_counts$sample == sample & is.na(walker_counts$WS_Count),]
    walker_counts_called_sample <- walker_counts_called[walker_counts_called$sample == sample,]
    walker_counts_extrapolated_sample <- data.frame(sample = sample, 
                                                    WS_Count = walker_counts_called_sample$WS_Count, 
                                                    n = walker_counts_uncalled_sample$n * walker_counts_called_sample$freq)
    walker_counts_extrapolated_sample$Extrapolated <- walker_counts_extrapolated_sample$n + walker_counts_called_sample$n
    #print(walker_counts_extrapolated_sample)
    walker_counts_extrapolated <- rbind(walker_counts_extrapolated, walker_counts_extrapolated_sample)
  }
  return(walker_counts_extrapolated)
}


plot_walker_counts <- function(walker_counts, mode = "absolute"){
  if(mode == "absolute"){
    plot <- ggplot(walker_counts, aes(x = sample, fill = WS_Count, y = n)) + 
      geom_bar(stat = "identity") + 
      theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 60,  hjust=1)) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_manual(values = walker_colors) + 
      labs(y = "Somatic variants (Autosomal)", fill = "Walker Score", x = "Sample")
  } else {
    plot <- ggplot(walker_counts, aes(x = sample, fill = WS_Count, y = freq)) + geom_bar(stat = "identity") +
      geom_text(aes(label = paste(round(freq*100,0), "%", sep =""), y=freq), position = position_stack(vjust = 0.5), size= 2) +
      theme_classic(base_size = 8) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_manual(values = walker_colors) +
      theme(axis.text.x = element_text(angle = 60,  hjust=1)) + 
      labs(y = "Somatic variants (% of total)", fill = "Walker Score", x = "Sample")
  }
  return(plot)
}

plot_Walker_Extrapolated <- function(Walker_Extrapolated){
  plot <- ggplot(Walker_Extrapolated, aes(x = sample, fill = WS_Count, y = Extrapolated)) + 
    geom_bar(stat = "identity", col = "black") +
    geom_text(aes(label = round(Extrapolated,0), y = Extrapolated), position = position_stack(vjust = 0.5), size= 2) +
    theme_classic(base_size = 8) + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_fill_manual(values = walker_colors) +
    theme(axis.text.x = element_text(angle = 60,  hjust=1)) + 
    labs(y = "Somatic variants (Extrapolated autosomes)", fill = "Walker Score", x = "Sample")
  return(plot)
}

plot_Walker_distribution <- function(walker_scores){
  plot <- ggplot(walker_scores, aes(x = sample, y = WS, fill = sample)) +
    geom_hline(yintercept = c(100, 1000), linetype = 2, col = "grey") +
    geom_boxplot() +  
    coord_cartesian(ylim = c(0,3000)) + 
    theme_classic(base_size = 8) + 
    theme(axis.text.x = element_text(angle = 60,  hjust=1), legend.position = "none") +
    labs(y = "Walker Score", x = "Sample")
  return(plot)
}

# This function performs all other functions and generates one page with all plots for one individual
makeWalkerReport <- function(walker_vcf_directory, individual, samples, output_dir){
  for(type in c("snv", "indel")){
    print(type)
    Walker_Scores <- RetrieveWalkerScores(walker_vcf_directory = walker_vcf_directory, samples = samples, type = type)
    Walker_Counts <- CalculateWalkerCounts(Walker_Scores)
    Walker_Extrapolated <- Extrapolate_Walker(walker_counts = Walker_Counts, walker_scores = Walker_Scores)
    
    p1 <- plot_walker_counts(walker_counts = Walker_Counts, mode = "absolute")
    p2 <- plot_walker_counts(walker_counts = Walker_Counts, mode = "relative")
    p3 <- plot_Walker_distribution(walker_scores = Walker_Scores)
    p4 <- plot_Walker_Extrapolated(Walker_Extrapolated)
    
    pdf(paste(output_dir, individual, "_Walker_Overview_",type,".pdf", sep = ""), width = 8, height = 8, pointsize = 8)
    plot_grid(p1, p2,  p3,  p4, ncol = 2)
    dev.off()
  }
}



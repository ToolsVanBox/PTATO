
CNV_cols <- c("steelblue4","red3","darkorange", "grey60")
names(CNV_cols) <- c("Gain","Loss","LOH", "Normal")

# Main function
Generate_Plots <- function(ReadCounts_100kb_file,
                           ReadCounts_1mb_file,
                           ReadCounts_segments_file,
                           BAF_binned_file,
                           BAF_segments_file,
                           CNV_file,
                           OUTPUT_dir, SAMPLE = "", max_level_to_plot = 4){

  # Read the data
  # ReadDepth
  ReadCounts_100kb <- read.delim(ReadCounts_100kb_file)
  ReadCounts_100kb$type <- "ReadDepth"
  ReadCounts_1mb <- read.delim(ReadCounts_1mb_file)
  ReadCounts_1mb$type <- "ReadDepth"

  # Temporary fix for the chrX
  ReadCounts_100kb$CopyNumber[ReadCounts_100kb$Chromosome == "X"] <- ReadCounts_100kb$CopyNumber[ReadCounts_100kb$Chromosome == "X"]/2
  ReadCounts_1mb$CopyNumber[ReadCounts_1mb$Chromosome == "X"] <- ReadCounts_1mb$CopyNumber[ReadCounts_1mb$Chromosome == "X"] / 2

  ReadCounts_segments <- read.delim(ReadCounts_segments_file)
  ReadCounts_segments$type <- "ReadDepth"
  # Temp fix for chrX
  ReadCounts_segments[ReadCounts_segments$chrom == "X",4] <-   ReadCounts_segments[ReadCounts_segments$chrom == "X",4] / 2

   # BAF
  BAF_binned <- read.delim(BAF_binned_file)
  BAF_binned$type <- "BAF"
  BAF_segments <- read.delim(BAF_segments_file)
  BAF_segments$type <- "BAF"

  # CNVs
  CNVs <- read.delim(CNV_file)

  ### Colors
  # CNV_cols <- c("darkblue","red3","darkorange")
  # names(CNV_cols) <- c("Gain","Loss","LOH")

  ReadCounts_plot <- ReadCounts_100kb[,c(1,2,3,5, 6)]
  names(ReadCounts_plot) <- c("chrom", "start.pos","end.pos","value","type")
  BAF_plot <- BAF_binned
  names(BAF_plot) <- c("chrom", "start.pos","end.pos","value","type")
  if(nrow(CNVs) > 0){
    CNVs_plot <- CNVs[,c(1,2,3,6)]
    names(CNVs_plot) <- c("chrom", "start.pos","end.pos","value")
  } else {
    CNVs_plot <- NULL
  }

  plot_karyogram(CNV_data = CNVs, ReadCounts = ReadCounts_100kb, SAMPLE = SAMPLE)
  ggsave(paste(OUTPUT_dir, ".karyogram.100kb.pdf", sep = ""), width = 10, height = 2.5, dpi = 300)
  ggsave(paste(OUTPUT_dir, ".karyogram.100kb.png", sep = ""), width = 10, height = 2.5, dpi = 300)

  plot_karyogram(CNV_data = CNVs, ReadCounts = ReadCounts_1mb, SAMPLE = SAMPLE, point_size = 0.7)
  ggsave(paste(OUTPUT_dir, ".karyogram.1mb.pdf", sep = ""), width = 10, height = 2.5, dpi = 300)
  ggsave(paste(OUTPUT_dir, ".karyogram.1mb.png", sep = ""), width = 10, height = 2.5, dpi = 300)



  plot_BAF(segment_data = ReadCounts_segments, segment_data2 = BAF_segments,
           point_data = ReadCounts_plot, point_data2 = BAF_plot, CNV_data = CNVs_plot)
  ggsave(paste(OUTPUT_dir, ".copynumber.baf.pdf", sep = ""), width = 11, height = 5, dpi = 300)
  ggsave(paste(OUTPUT_dir, ".copynumber.baf.png", sep = ""), width = 11, height = 5, dpi = 300)


  plot_CN(segment_data = ReadCounts_segments,
          point_data = ReadCounts_plot, title = "Read depth",
          CNV_data = CNVs_plot)

  ggsave(paste(OUTPUT_dir, ".copynumber.pdf", sep = ""), width = 10, height = 9, dpi = 300)
  ggsave(paste(OUTPUT_dir, ".copynumber.png", sep = ""), width = 10, height = 9, dpi = 300)


  plot_CN(segment_data = BAF_segments,
          point_data = BAF_plot, title = "B-allele frequency",
          CNV_data = CNVs_plot)

  ggsave(paste(OUTPUT_dir, ".baf.pdf", sep = ""), width = 10, height = 9, dpi = 300)
  ggsave(paste(OUTPUT_dir, ".baf.png", sep = ""), width = 10, height = 9, dpi = 300)
}


plot_karyogram <- function(CNV_data = NULL, ReadCounts, SAMPLE = "",
                           max_level_to_plot = 4, SVtype_column = 6,
                           point_size = 0.3){
  print(SAMPLE)
  ReadCounts$Type <- "Normal"

  if(!is.null(CNV_data)){
    CNV_g <- GRanges(seqnames = CNV_data[,1], IRanges(start = CNV_data[,2], end = CNV_data[,3]), Type = CNV_data[,SVtype_column])
    ReadCounts_g <- GRanges(seqnames = ReadCounts[,1], IRanges(start = ReadCounts[,2], end = ReadCounts[,3]))
    Olap_ReadCounts_CNV <- findOverlaps(ReadCounts_g, CNV_g)

    ReadCounts$Type[queryHits(Olap_ReadCounts_CNV)] <- CNV_g$Type[subjectHits(Olap_ReadCounts_CNV)]
  }

  ReadCounts$CopyNumber[which(ReadCounts$CopyNumber > max_level_to_plot)] <- max_level_to_plot
  ReadCounts$Chromosome <- factor(ReadCounts$Chromosome, levels = c(1:22, "X", "Y"))
  plot <- ggplot(data = ReadCounts) +

    # points for the 2n (normal) regions:
    geom_point(data = ReadCounts[which(ReadCounts$Type == "Normal"),], aes(x = (start.pos+end.pos) / 2, y = CopyNumber, color = Type), alpha = 0.7, size = point_size) +
    # points for the CNVs:
    geom_point(data = ReadCounts[which(ReadCounts$Type != "Normal"),], aes(x = (start.pos+end.pos) / 2, y = CopyNumber, color = Type), alpha = 0.9, size = point_size+0.1) +

    facet_grid(.~Chromosome, space = "free_x", scales = "free", switch="y") +
    scale_color_manual(values = CNV_cols) +

    labs(x = "Genomic Position (bp)", y = "Copy Number") + ggtitle(SAMPLE) +
    theme_grey(base_size = 10) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(
      axis.line.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),

      panel.grid.minor=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(colour="grey90", linetype=2, size=0.2),
      panel.spacing = unit(0.2, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position="bottom",
      panel.background=element_blank(),
      panel.border=element_blank(),
      strip.placement = "outside")
  return(plot)
}

plot_BAF <- function(segment_data, segment_data2 = NULL,
                     point_data, point_data2 = NULL,
                     CNV_data = NULL,
                     chrom = "", start = "", end = "", max_level_to_plot = 4){
  names(segment_data)[ncol(segment_data)-1] <- "value"
  if(!is.null(segment_data2)){
    names(segment_data2)[ncol(segment_data2)-1] <- "value"
    merged_segments <- rbind(segment_data[,c("chrom", "start.pos", "end.pos", "value", "type")],
                             segment_data2[,c("chrom", "start.pos", "end.pos", "value", "type")])
    merged_segments$type <- factor(merged_segments$type, levels = c(unique(segment_data$type), unique(segment_data2$type)))
  } else {
    merged_segments <- segment_data[,c("chrom", "start.pos", "end.pos", "value", "type")]
  }

  # Determine the vertical parts of the segments (y to yend)
  segment_plot_data <- data.frame()
  for(type in unique(merged_segments$type)){
    print(type)
    segment_data_type <- merged_segments[merged_segments$type == type,]

    for(chromosome in unique(segment_data_type$chrom)){
      chrom_data <- segment_data_type[segment_data_type$chrom == chromosome,]
      for(i in 1:nrow(chrom_data)){
        #print(i)
        if(i < nrow(chrom_data)){
          chrom_data$ymin[i] <- chrom_data$value[i+1]
        } else {
          chrom_data$ymin[i] <- chrom_data$value[i]
        }
      }
      segment_plot_data <- rbind(segment_plot_data, chrom_data)
    }
  }
  segment_plot_data$chrom <- factor(segment_plot_data$chrom, levels = c(1:22,"X","Y"))
  segment_plot_data <- segment_plot_data[order(segment_plot_data$type, segment_plot_data$chrom),]

  names(point_data)[ncol(point_data) - 1] <- "value"
  if(!is.null(point_data2)){
    names(point_data2)[ncol(point_data2) - 1] <- "value"
    merged_data <- rbind(point_data, point_data2)
    merged_data$type <- factor(merged_data$type, levels = c(unique(point_data$type), unique(point_data2$type)))

  } else {
    merged_data <- point_data
  }

  merged_data$value[merged_data$value > max_level_to_plot] <- max_level_to_plot
  merged_data$chrom <- factor(merged_data$chrom, levels = c(1:22, "X", "Y"))
  merged_data <- merged_data[order(merged_data$type, merged_data$chrom),]

  if(!is.null(CNV_data)){
    CNV_data$chrom <- factor(CNV_data$chrom, levels = c(1:22, "X", "Y"))
    CNV_data <- CNV_data[order(CNV_data$chrom),]
  }

  plot <- ggplot(data = segment_plot_data) +
    geom_rect(data = CNV_data, mapping = aes(xmin = start.pos, xmax = end.pos, ymin = -Inf, ymax = Inf, fill = value), alpha = 0.2) +

    geom_point(data = merged_data, aes(x = start.pos, y = value), color = "darkgray", alpha = 0.8, size = 0.3) +

    geom_segment(aes(x = start.pos, xend = end.pos, y = value, yend = value), col = "red") +
    geom_segment(aes(x = end.pos, xend = end.pos, y = ymin, yend = value), col = "red") +

    #geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) +
    facet_grid(type~chrom, space = "free_x", scales = "free", switch="y") +
    scale_fill_manual(values = CNV_cols) +
    labs(x = "Genomic Position (bp)", y = "") +
    theme_grey(base_size = 8) +
    theme(
      axis.line.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),

      panel.grid.minor=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(colour="grey90", linetype=2, size=0.2),
      panel.spacing = unit(0.2, "lines"),

      #legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      strip.placement = "outside")
  return(plot)
}




### test
### could make this one function. main difference if facet_wrap vs facet_grid
plot_CN <- function(segment_data, segment_data2 = NULL,
                     point_data, point_data2 = NULL,
                     CNV_data = NULL,
                     chrom = "", start = "", end = "", max_level_to_plot = 4, nrow = 5, title = ""){
  names(segment_data)[ncol(segment_data)-1] <- "value"
  if(!is.null(segment_data2)){
    names(segment_data2)[ncol(segment_data2)-1] <- "value"
    merged_segments <- rbind(segment_data[,c("chrom", "start.pos", "end.pos", "value", "type")], segment_data2[,c("chrom", "start.pos", "end.pos", "value", "type")])
    merged_segments$type <- factor(merged_segments$type, levels = c(unique(segment_data$type), unique(segment_data2$type)))
  } else {
    merged_segments <- segment_data[,c("chrom", "start.pos", "end.pos", "value", "type")]
  }

  # Determine the vertical parts of the segments (y to yend)
  segment_plot_data <- data.frame()
  for(type in unique(merged_segments$type)){
    print(type)
    segment_data_type <- merged_segments[merged_segments$type == type,]

    for(chromosome in unique(segment_data_type$chrom)){
      chrom_data <- segment_data_type[segment_data_type$chrom == chromosome,]
      for(i in 1:nrow(chrom_data)){
        #print(i)
        if(i < nrow(chrom_data)){
          chrom_data$ymin[i] <- chrom_data$value[i+1]
        } else {
          chrom_data$ymin[i] <- chrom_data$value[i]
        }
      }
      segment_plot_data <- rbind(segment_plot_data, chrom_data)
    }
  }
  segment_plot_data$chrom <- factor(segment_plot_data$chrom, levels = c(1:22,"X","Y"))
  segment_plot_data <- segment_plot_data[order(segment_plot_data$type, segment_plot_data$chrom),]

  #print(head(segment_plot_data))
  names(point_data)[ncol(point_data) - 1] <- "value"
  if(!is.null(point_data2)){
    names(point_data2)[ncol(point_data2) - 1] <- "value"
    merged_data <- rbind(point_data, point_data2)
    merged_data$type <- factor(merged_data$type, levels = c(unique(point_data$type), unique(point_data2$type)))

  } else {
    merged_data <- point_data
  }

  merged_data$value[merged_data$value > max_level_to_plot] <- max_level_to_plot
  merged_data$chrom <- factor(merged_data$chrom, levels = c(1:22, "X", "Y"))
  merged_data <- merged_data[order(merged_data$type, merged_data$chrom),]

  #print(head(merged_data))

  if(!is.null(CNV_data)){
    CNV_data$chrom <- factor(CNV_data$chrom, levels = c(1:22, "X", "Y"))
    CNV_data <- CNV_data[order(CNV_data$chrom),]
  }

  plot <- ggplot(data = segment_plot_data) +
    geom_rect(data = CNV_data, mapping = aes(xmin = start.pos, xmax = end.pos, ymin = -Inf, ymax = Inf, fill = value), alpha = 0.2) +

    geom_point(data = merged_data, aes(x = start.pos, y = value), color = "darkgray", alpha = 0.8, size = 0.3) +

    geom_segment(aes(x = start.pos, xend = end.pos, y = value, yend = value), col = "red") +
    geom_segment(aes(x = end.pos, xend = end.pos, y = ymin, yend = value), col = "red") +

    #geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) +
    facet_wrap(~chrom,  scales = "free_x", ncol = nrow) +
    scale_fill_manual(values = CNV_cols) +
    labs(x = "Genomic Position (bp)", y = "") +
    ggtitle(title) +
    theme_grey(base_size = 8) +
    theme(
      axis.line.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),

      panel.grid.minor=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(colour="grey90", linetype=2, size=0.2),
      panel.spacing = unit(0.2, "lines"),

      #legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      strip.placement = "outside", plot.title = element_text(hjust = 0.5))
  return(plot)
}

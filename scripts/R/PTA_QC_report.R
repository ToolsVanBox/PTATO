# input variables:
# at least 1: the output directory,
# optionaly: 2 is name for report (default is the dir name), 3 is output directory

# One major drawback of single cell WGS is allelic and locus dropout, leading to regions without coverage. Of course variants can never be detected in regions that are not covered by the WGS. Therefore we need to know how much of the genome is amplified.
# In this script we'll look at the WGS metrics data produced by Picard.
# ------------------------------

# libraries
library(cowplot)
library(ggplot2)
library(scales)

# get input arguments
args <- commandArgs(trailingOnly = TRUE)
samples <- strsplit(args[1],",")[[1]]
alignment_files <- strsplit(args[2],",")[[1]]
wgsmetrics_files <- strsplit(args[3],",")[[1]]

order_i = order(samples)
samples = samples[order_i]
alignment_files = alignment_files[order_i]
wgsmetrics_files = wgsmetrics_files[order_i]
n_wgsm_files = length(wgsmetrics_files)
outfile = args[4]

# setwd("~/hpc/pmc_vanboxtel/projects/PTA_manuscript/2_Code/work/01/b9b7fd18a81106cd5d8553e1ba5378/")

# set labels size based on the number of samples
if (n_wgsm_files <= 10){
    label_size = 3.5
    axis_text_size = 8.8
} else if (n_wgsm_files <= 20){
    label_size = 2.5
    axis_text_size = 7
} else{
    label_size = 1.5
    axis_text_size = 4.5
}

# Create dataframe matching samples to a number
sample_nrs = seq_len(n_wgsm_files)
sample_nr_df = data.frame("sample" = samples, "sample_nr" = factor(sample_nrs, levels = sample_nrs))

# Set colours
set.seed(42)
palette = scales::hue_pal(c(0, 360) + 30, l = 50, c = 75)(n_wgsm_files)
palette = palette[sample.int(n_wgsm_files, n_wgsm_files)]


# Log input info
# cat("current folder: ", getwd(), fill = T)
# cat("input = 1: ", in_dir, "| 2: ", outfile, "| 3: ", out_dir, fill = T)

# import the wgsmetrics files
wgsmetrics_merged_l <- lapply(seq(length(wgsmetrics_files)), function(n_samp) {
    wgsmetrics <- read.delim(wgsmetrics_files[n_samp], skip = 11, header = F)
    wgsmetrics$sample <- samples[n_samp]
    wgsmetrics$cumsum <- (sum(wgsmetrics[,2]) - cumsum(as.numeric(wgsmetrics[,2]))) / sum(wgsmetrics[,2]) * 100
    #wgsmetrics$percent = 100 * wgsmetrics$V2 / sum(wgsmetrics$V2)
    return(wgsmetrics)
})
wgsmetrics_merged = do.call(rbind, wgsmetrics_merged_l)
wgsmetrics_merged$sample <- factor(wgsmetrics_merged$sample, levels = samples)
wgsmetrics_merged = merge(wgsmetrics_merged, sample_nr_df, by = "sample", sort = FALSE)

# Next look at how many mapped reads are discarded for example because they're PCR duplicates
wgsmetrics_merged2_l <- lapply(seq(length(wgsmetrics_files)), function(n_samp) {
    wgsmetrics2 <- read.delim(wgsmetrics_files[n_samp], skip = 6, header = T, nrows = 1)
    wgsmetrics2$sample <- samples[n_samp]
    wgsmetrics2
})
wgsmetrics_merged2 = do.call(rbind, wgsmetrics_merged2_l)
wgsmetrics_merged2$sample <- factor(wgsmetrics_merged2$sample, levels = samples)
wgsmetrics_merged2$PCT_1X_2 = wgsmetrics_merged2$PCT_1X * 100
wgsmetrics_merged2$PCT_5X_2 = wgsmetrics_merged2$PCT_5X * 100
wgsmetrics_merged2 = merge(wgsmetrics_merged2, sample_nr_df, by = "sample", sort = FALSE)

# Check how many sequenced and mapped bases are discarded due to quality control
select_columns = c("sample","PCT_EXC_MAPQ", "PCT_EXC_DUPE","PCT_EXC_TOTAL")
wgsmetrics_merged3 <- wgsmetrics_merged2[,select_columns]
wgsmetrics_merged3$PCT_EXC_OTHER <- wgsmetrics_merged3$PCT_EXC_TOTAL - wgsmetrics_merged3$PCT_EXC_DUPE - wgsmetrics_merged3$PCT_EXC_MAPQ
rownames(wgsmetrics_merged3) = wgsmetrics_merged3$sample
wgsmetrics_merged3 = wgsmetrics_merged3[,-c(1, 4), drop = FALSE]
wgsmetrics_merged3_m = reshape(wgsmetrics_merged3, direction = "long", idvar = "sample",
        timevar = "variable", varying = list(names(wgsmetrics_merged3)), times = names(wgsmetrics_merged3), ids = row.names(wgsmetrics_merged3))
colnames(wgsmetrics_merged3_m)[2] = "value"
wgsmetrics_merged3_m$sample <- factor(wgsmetrics_merged3_m$sample, levels = rev(samples))
wgsmetrics_merged3_m = merge(wgsmetrics_merged3_m, sample_nr_df, by = "sample", sort = FALSE)
wgsmetrics_merged3_m$sample_nr = factor(wgsmetrics_merged3_m$sample_nr, levels = rev(levels(wgsmetrics_merged3_m$sample_nr)))

# Read alignment file for nr of reads and mismatch rate
alignment_metrics_merged_l <- lapply(seq(length(alignment_files)), function(n_samp) {
    alignment_metrics <- read.delim(alignment_files[n_samp], skip = 6, header = T, nrows = 1)
    alignment_metrics$sample <- samples[n_samp]
    alignment_metrics
})
alignment_metrics_merged = do.call(rbind, alignment_metrics_merged_l)

# Check how many (unmapped/filtered) reads there are
alignment_metrics_merged$Filtered = alignment_metrics_merged$TOTAL_READS - alignment_metrics_merged$PF_READS
alignment_metrics_merged$PF_Unaligned = alignment_metrics_merged$PF_READS - alignment_metrics_merged$PF_READS_ALIGNED
read_counts_df = alignment_metrics_merged[,c("Filtered", "PF_Unaligned", "PF_READS_ALIGNED"), drop = FALSE]
rownames(read_counts_df) = alignment_metrics_merged$sample
colnames(read_counts_df)[3] = "PF_Aligned"
read_counts_df_long = reshape(read_counts_df, direction = "long", idvar = "sample",
        timevar = "Reads", varying = list(names(read_counts_df)), times = names(read_counts_df), ids = row.names(read_counts_df))
colnames(read_counts_df_long)[2] = "Count"
read_counts_df_long$Reads = factor(read_counts_df_long$Reads, levels = c("Filtered", "PF_Unaligned", "PF_Aligned"))
read_counts_df_long = merge(read_counts_df_long, sample_nr_df, by = "sample", sort = FALSE)


# Check how many (indel) mismatches and HQ errors there are
error_rate_df = alignment_metrics_merged[,c("PF_MISMATCH_RATE", "PF_HQ_ERROR_RATE", "PF_INDEL_RATE"), drop = FALSE]
rownames(error_rate_df) = alignment_metrics_merged$sample
colnames(error_rate_df) = c("PF_Mismatch", "PF_HQ_ERROR", "PF_INDEL")
error_rate_df_long = reshape(error_rate_df, direction = "long", idvar = "sample",
                              timevar = "Error", varying = list(names(error_rate_df)), times = names(error_rate_df), ids = row.names(error_rate_df))
colnames(error_rate_df_long)[2] = "Rate"
error_rate_df_long$Error = factor(error_rate_df_long$Error, levels = c("PF_Mismatch", "PF_HQ_ERROR", "PF_INDEL"))
error_rate_df_long = merge(error_rate_df_long, sample_nr_df, by = "sample", sort = FALSE)


# Combine all qc data to write out
all_merged_metrics = wgsmetrics_merged2[,c("sample", "PCT_1X", "PCT_5X", "MEAN_COVERAGE", "HET_SNP_Q", "PCT_EXC_DUPE", "PCT_EXC_MAPQ")]
all_merged_metrics$PCT_EXC_OTHER <- wgsmetrics_merged2$PCT_EXC_TOTAL - wgsmetrics_merged2$PCT_EXC_DUPE - wgsmetrics_merged2$PCT_EXC_MAPQ
read_counts_df$sample = rownames(read_counts_df)
all_merged_metrics = merge(all_merged_metrics, read_counts_df, by = "sample", sort = FALSE)
error_rate_df$sample = rownames(error_rate_df)
all_merged_metrics = merge(all_merged_metrics, error_rate_df, by = "sample", sort = FALSE)
all_merged_metrics = merge(all_merged_metrics, sample_nr_df, by = "sample", sort = FALSE)

outfile_txt = gsub(".pdf", ".txt", outfile)
write.table(all_merged_metrics, outfile_txt, quote = FALSE, sep = "\t", row.names = FALSE)


## Plot coverage
p1 = ggplot(wgsmetrics_merged, aes(x = V1, y = cumsum, col = sample_nr)) +
    geom_line(lwd = 0.8) +
    geom_vline(xintercept = 15) +
    theme_bw() +
    scale_color_manual(values = palette) +
    scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
    scale_x_continuous(expand = c(0,0)) +
    coord_cartesian(xlim= c(0, 50)) +
    labs(y = "Percentage of bases", x = "Coverage (X)") +
    theme(legend.position = 'none')

# ggplot(wgsmetrics_merged[wgsmetrics_merged$V1 != 0,], aes(x = V1, y = percent, color = sample_nr)) +
#     geom_line(lwd = 0.8) +
#     theme_bw() +
#     geom_vline(xintercept = 15) +
#     theme_bw() +
#     scale_color_manual(values = palette) +
#     #scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
#     scale_x_continuous(expand = c(0,0)) +
#     coord_cartesian(xlim= c(1, 50)) +
#     labs(y = "Percentage of bases", x = "Coverage (X)") +
#     theme(legend.position = 'none')



# Next quantify how much of the genome has more than 0 or 5 coverage
p2 = ggplot(wgsmetrics_merged2, aes(x = sample_nr, y = PCT_1X_2, fill = sample)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(PCT_1X_2, 1)), vjust=-0.25, size = label_size) +
    scale_y_continuous(expand = c(0,0), limits = c(0,105)) +
    scale_fill_manual(values = palette) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = axis_text_size)) +
    labs(y = "% of the genome > 0X coverage", x = "")

p3 = ggplot(wgsmetrics_merged2, aes(x = sample_nr, y = PCT_5X_2, fill = sample)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_text(aes(label = round(PCT_5X_2, 1)), vjust=-0.25, size = label_size) +
    scale_y_continuous(expand = c(0,0), limits = c(0,105)) +
    scale_fill_manual(values = palette) +
    theme_classic() +
    labs(y = "% of the genome >= 5X coverage", x = "") +
    theme(legend.position = 'none',
          axis.text.x = element_text(size = axis_text_size))

# plot the mean and SD coverage per sample
p4 = ggplot(wgsmetrics_merged2, aes(x = sample_nr, y = MEAN_COVERAGE, fill = sample)) +
    geom_bar(stat= "identity") +
    geom_errorbar(aes(ymin = MEAN_COVERAGE-SD_COVERAGE, ymax = MEAN_COVERAGE+SD_COVERAGE), width=.2)+
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = palette) +
    labs(y = 'mean coverage', x = '') +
    theme(legend.position = 'none',
          axis.text.x = element_text(size = axis_text_size))

# Picard also estimates the sensitivity to detect heterozygous SNVs, based on the coverage
p5 = ggplot(wgsmetrics_merged2, aes(x = sample_nr, y = HET_SNP_SENSITIVITY, fill = sample)) +
    geom_text(aes(label = round(HET_SNP_SENSITIVITY, 2)), vjust=-0.25, size = label_size) +
    geom_bar(stat= "identity") +
    theme_classic() +
    scale_y_continuous(expand = c(0,0), limits = c(0,1.05)) +
    scale_fill_manual(values = palette) +
    labs(x = '', y = "heterozygous SNP sensitivity") +
    theme(legend.position = 'none',
          axis.text.x = element_text(size = axis_text_size))

# Settings used for the next three plots
qc_colors = c("#E41A1C", "#377EB8", "#4DAF4A")
small_legend = theme(legend.key.size = unit(0.3, "cm"),
                     legend.spacing.x = unit(0.05, 'cm'),
                     legend.spacing.y = unit(0.05, "cm"),
                     legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
                     legend.text = element_text(size = 7),
                     legend.title = element_text(size = 9))

# the number of reads discarded due to quality
p6 = ggplot(wgsmetrics_merged3_m, aes(x = sample_nr, fill = variable, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Sample", y = "Filtered bases" , fill = "Exclusion reason") +
    theme_classic() +
    scale_y_continuous(expand = c(0,0), labels = scales::percent) +
    scale_fill_manual(values = qc_colors) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(axis.text.y = element_text(size = axis_text_size)) +
    small_legend

# Plot (unmapped/filtered) read counts
p7 = ggplot(read_counts_df_long, aes(fill = Reads, y = Count, x = sample_nr)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(x = "Sample", y = "Read count") +
    scale_fill_manual(values = qc_colors) +
    theme_classic() +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(axis.text.x = element_text(size = axis_text_size)) +
    small_legend


# Plot mismatches and error rate
p8 = ggplot(error_rate_df_long, aes(x = sample_nr, fill = Error, y = Rate)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_y_continuous(expand = c(0,0), labels = scales::percent) +
    labs(x = "Sample", y = "Error rate") +
    scale_fill_manual(values = qc_colors) +
    theme_classic() +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(axis.text.x = element_text(size = axis_text_size)) +
    small_legend

# Create legend with numbering
wgs_metrics_nr_name = wgsmetrics_merged2
wgs_metrics_nr_name$sample = paste0(seq_len(n_wgsm_files), ": ", wgs_metrics_nr_name$sample)
wgs_metrics_nr_name$sample = factor(wgs_metrics_nr_name$sample, levels = unique(wgs_metrics_nr_name$sample))
p9 = ggplot(wgs_metrics_nr_name, aes(x = sample, y = PCT_5X, fill = sample)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = palette) +
    theme_classic() +
    labs(fill = "Sample")


## plot and save
legend = get_legend(p9)

all_plots = plot_grid(
    legend, p1, p2, p3, p4, p5, p6, p7,
    ncol = 2,
    rel_heights = c(1.4, 1, 1, 1)
)

all_plots = plot_grid(
    all_plots, p8,
    nrow = 2,
    rel_heights = c(4.2, 1)
)

# ggsave(file.path(out_dir, paste0(outfile, ".pdf")), plot = all_plots, width = 8.3, height = 11.7)
ggsave(outfile, plot = all_plots, width = 8.3, height = 16)

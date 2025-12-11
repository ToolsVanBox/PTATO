.getSourceDir <- function() {
  cmdArgs <- commandArgs(trailingOnly=FALSE)
  fileArg <- "--file="
  match <- grep(fileArg, cmdArgs)
  if (length(match) > 0) {
    sourcedir <- dirname(gsub(fileArg, "", cmdArgs[match]))
  }
  return( sourcedir )
}
sourcedir <- .getSourceDir()

source(paste(sourcedir,"/functions/PTATOplus.R",sep=""), chdir = T)

### libraries
library(cowplot)
library(ggplot2)
library(scales)
library(VariantAnnotation)

# get input arguments
args <- commandArgs(trailingOnly = TRUE)
samples <- strsplit(args[1],",")[[1]]
snv_vcfs <- strsplit(args[2],",")[[1]]
snv_filtered_vcfs <- strsplit(args[3],",")[[1]]
walker_vcfs <- strsplit(args[4],",")[[1]]
ptato_table_files <- strsplit(args[5],",")[[1]]
autosomal_files <- strsplit(args[6],",")[[1]]
out_prefix = args[7]

# Get output variables
out_file <- paste0(out_prefix, ".pdf", sep = "")
out_txt <- paste0(out_prefix, ".txt", sep = "")
out <- data.frame(sample = samples)
order_i = order(samples)
samples = samples[order_i]
n_samples <- length(samples)

# set labels size based on the number of samples
if (n_samples <= 10){
  label_size = 3.5
  axis_text_size = 8.8
} else if (n_samples <= 20){
  label_size = 2.5
  axis_text_size = 7
} else{
  label_size = 1.5
  axis_text_size = 4.5
}

# Create dataframe matching samples to a number
sample_nrs = seq_len(n_samples)
sample_nr_df = data.frame("sample" = samples, "sample_nr" = factor(sample_nrs, levels = sample_nrs))

# Set colours
set.seed(42)
palette = scales::hue_pal(c(0, 360) + 30, l = 50, c = 75)(n_samples)
palette = palette[sample.int(n_samples, n_samples)]
names(palette) <- samples

# Read in all autosomal callable loci files 
CHROMOSOMES <- c(1:22, "X", "Y")
AUTOSOMES <- c(1:22)
callableloci_overview <- data.frame()
for (callability_file in autosomal_files){
  callability <- read.delim(callability_file, header = T)
  sample = gsub(x = callability_file, pattern = ".callableloci.autosomal.txt", replacement = "")
  callableloci_overview <- rbind(callableloci_overview, data.frame(sample = sample, CallableLoci = callability$Freq[callability$State == "CALLABLE"]))
}
# Create callable loci plot 
callableloci_overview <- merge(callableloci_overview, sample_nr_df, by = "sample")
callable_plot <- ggplot(callableloci_overview, aes(x = sample_nr, y = CallableLoci)) +
  geom_bar(stat = "identity", fill = palette, width = 0.8) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1), position = "right") +
  coord_flip() +
  labs(x = "Sample") +
  theme_bw(base_size = 6) 

# Create a legend plot 
callableloci_overview$Label <- paste(callableloci_overview$sample_nr, callableloci_overview$sample, sep = ": ")
callableloci_overview$Label <- factor(callableloci_overview$Label, levels = callableloci_overview$Label[order(callableloci_overview$sample_nr)])
palette2 <- palette
names(palette2) <- callableloci_overview$Label
legend_plot <- ggplot(callableloci_overview, aes(x = Label, y = CallableLoci, fill = Label)) +
  geom_bar(stat = "identity", width = 0.8) + scale_fill_manual(values = palette2)
out2 <- merge(out, callableloci_overview[,c("sample", "CallableLoci")], by = "sample")
out2$CallableLoci <- round(out2$CallableLoci,3)

# Retrieve the SNVs from both the filtered as the unfiltered vcf file 
MinimalVAF = 0.2 
SBSs_raw = list()
for(sample in samples){
  SBSs_raw[[sample]] <- readPTATOvcf(vcf_unfiltered = snv_vcfs[grep(paste0(sample, ".snvs.ptato.vcf.gz"), snv_vcfs)], 
                                          vcf_filtered = snv_filtered_vcfs[grep(paste0(sample, ".snvs.ptato.filtered.vcf.gz"), snv_filtered_vcfs)], 
                                          VAF_threshold = MinimalVAF)
}

# Create the plots fro the PTA probability 
p1 <- plot_PTAprobs_cumulative(SBSs_raw, Labels = sample_nr_df$sample_nr, Colors = palette)
p2 <- plot_PTAprobs_distribution(SBSs_raw, Labels = sample_nr_df$sample_nr, Colors = palette)
p3 <- plot_PTAprobs_filter(SBSs_raw, mode = "absolute", Labels = sample_nr_df$sample_nr)
p4 <- plot_PTAprobs_filter(SBSs_raw, mode = "relative", Labels = sample_nr_df$sample_nr)


### Get all the PTA probabilities
Probs <- data.frame()
for(sample in names(SBSs_raw)){
  print(sample)
  probs_sample <- data.frame(sample = sample, 
                             PTAprob = SBSs_raw[[sample]]$PTAprob, 
                             Chr = as.character(seqnames(SBSs_raw[[sample]])),
                             PTAprobCutoff_LiRA = SBSs_raw[[sample]]$PTAprobCutoff_Walker, 
                             PTAprobCutoff_Cossim = SBSs_raw[[sample]]$PTAprobCutoff_Cossim, 
                             PTAprobCutoff = SBSs_raw[[sample]]$PTAprobCutoff, 
                             TotalVariants = length(SBSs_raw[[sample]]),
                             Filter = SBSs_raw[[sample]]$FILTER)
  Probs <- rbind(Probs, probs_sample)
}

# Filter the probabilities
Probs_freq <- Probs[Probs$Chr %in% paste("chr", 1:22, sep = ""),] %>% group_by(sample, Filter, PTAprobCutoff_LiRA, PTAprobCutoff_Cossim, PTAprobCutoff) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
PASS <- Probs_freq[Probs_freq$Filter == "PASS",]
colnames(PASS)[colnames(PASS) == "n"] <- "SNVs_PASS"
FAIL <- Probs_freq[Probs_freq$Filter == "FAIL",]
colnames(FAIL)[colnames(FAIL) == "n"] <- "SNVs_FAIL"
out3 <- merge(out2, PASS[,c("sample", "PTAprobCutoff_LiRA","PTAprobCutoff_Cossim","PTAprobCutoff","SNVs_PASS")], by = "sample")
out3 <- merge(out3, FAIL[,c("sample", "SNVs_FAIL")], by = "sample")
out3$SNV_burden <- round(out3$SNVs_PASS / out3$CallableLoci, 0)
prec_recall_plot <- plot_precision_recall_multisample(ptatotable_files = ptato_table_files, samples = samples)

### Walker / linked read plots
Walker_Scores <- RetrieveWalkerScores(walker_vcf_files =walker_vcfs, samples = samples, type = "snv")
Walker_Counts <- CalculateWalkerCounts(Walker_Scores)
Walker_Extrapolated <- Extrapolate_Walker(walker_counts = Walker_Counts, walker_scores = Walker_Scores)

Walker_PASS <- Walker_Extrapolated[Walker_Extrapolated$WS_Count %in% c("Intermediate", "High"),] %>% group_by(sample) %>%
  dplyr::summarise(SNVs_Walker = round(sum(Extrapolated), 0))
out4 <- merge(out3, Walker_PASS, by = "sample")

p5 <- plot_walker_counts(walker_counts = Walker_Counts, mode = "absolute")
p6 <- plot_walker_counts(walker_counts = Walker_Counts, mode = "relative")
p7 <- plot_Walker_distribution(walker_scores = Walker_Scores)
p8 <- plot_Walker_Extrapolated(Walker_Extrapolated)

## plot and save
legend = get_legend(legend_plot)
# Page1 
plots_page1 = plot_grid(
  legend, callable_plot, p1, p2, p3, p4,
  ncol = 2
)
# Page 3 
plots_page3 = plot_grid( p5, p6, p7, p8,
                         ncol = 2
)

# Final plot 
pdf(out_file, width = 8.3, height = 9)
plots_page1
prec_recall_plot
plots_page3
dev.off()
# Final report 
write.table(x = out4, file = out_txt, sep = "\t", quote = F, row.names = F)


library(VariantAnnotation)
# library(ggplot2)
# library(gridExtra)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
somatic_vcf_fname <- args[1]
germline_vcf_fname <- args[2]
phased_germline_vcf_fname <- args[3]
chrom <- args[4]
bulkname <- args[5]
output_file <- args[6]
flank <- 200000

somatic_vcf <- readVcf(somatic_vcf_fname)
somatic_vcf_chr <- somatic_vcf[seqnames(somatic_vcf) == chrom,]

AB_somatic <- function( somatic_mut ) {
  rngs <- GRanges(as.character(seqnames(somatic_mut)), IRanges(start(somatic_mut)-flank,start(somatic_mut)+flank))
  param <- ScanVcfParam( which = rngs )

  germline_vcf <- readVcf(germline_vcf_fname, "hg38", param)
  phased_germline_vcf <- readVcf( phased_germline_vcf_fname, "hg38", param)

  germline_vcf_phased_pos <- germline_vcf[which(start(germline_vcf) %in% start(phased_germline_vcf)),]
  phased_germline_vcf <- phased_germline_vcf[which(start(phased_germline_vcf) %in% start(germline_vcf_phased_pos)),]

  # Add phasing info
  geno(germline_vcf_phased_pos)$PGT <- geno(phased_germline_vcf)$GT

  # Filter for variants with a rs number that are heterozygous in the bulk
  # germline_vcf_phased_pos_rshet <- germline_vcf_phased_pos[which((grepl("rs",names(germline_vcf_phased_pos)) & geno(germline_vcf_phased_pos)$PGT[,bulkname] != "0|0" & geno(germline_vcf_phased_pos)$PGT[,bulkname] != "1|1"& geno(germline_vcf_phased_pos)$PGT[,bulkname] != ".|.")),]
  germline_vcf_phased_pos_rshet <- germline_vcf_phased_pos
  if (length(germline_vcf_phased_pos_rshet) == 0 | length(phased_germline_vcf) == 0 | length(somatic_vcf_chr) == 0) {
    message("Not enough mutations")
    return(0)
  }

  # Determine samples where variant is called
  gt <- geno(somatic_mut)$GT
  samples_present <- colnames(gt)[gt == "0|1" | gt == "0/1" | gt == "1|1" | gt == "1/1"]
  ab_tables_l <- lapply(samples_present, function(samplename) {plot_AB_sample(germline_vcf_phased_pos_rshet, somatic_mut, samplename)})
  ab_table <- do.call(rbind, ab_tables_l)

  message("Create out table")
  return( na.omit( ab_table ) )
}

plot_AB_sample = function(vcf4, vcf_somatic, samplename) {
  var_start = start(vcf_somatic)
  var_end = end(vcf_somatic)
  varchr = as.character(seqnames(vcf_somatic))
  # filter out empty and non heterozygous variants
  vcf5 <- vcf4[start(vcf4) == var_start | (geno(vcf4)$PGT[,samplename] != ".|." & geno(vcf4)$PGT[,samplename] != "1|1" & geno(vcf4)$PGT[,samplename] != "0|0" & geno(vcf4)$PGT[,samplename] != ".")]

  if (length(vcf5) == 0) {
    message("Not enough mutations after filtering out empty and non heterozygous variants")
    return(0)
  }
  # Determine allele of variants
  allele <- do.call(rbind, strsplit(geno(vcf5)$PGT[,samplename],"|"))[,1]

  # Switch ad for variants phased to allele 1.
  ad <- do.call(rbind, geno(vcf5)$AD[,samplename])
  tmp_ad <- ad[allele == "1",2]
  ad[allele == "1",2] = ad[allele == "1",1]
  ad[allele == "1",1] = tmp_ad
  # tmp_ad <- ad[allele == "1",][,2]
  # ad[allele == "1",][,2] = ad[allele == "1",][,1]
  # ad[allele == "1",][,1] = tmp_ad

  # Calculate vaf
  vaf <- ad[,2]/rowSums(ad)

  # Create df for all variants
  mydf <- data.frame( "a0" = ad[,1], "a1" = ad[,2], "total_ad" = rowSums(ad), "allele" = allele, "vaf" = vaf, "start" = start(vcf5), "end" = end(vcf5))
  mydf <- mydf[which(mydf$vaf < 1 & mydf$vaf > 0),]

  ad_somatic = geno(vcf_somatic)$AD[,samplename][[1]]
  vardf <- data.frame("a0" = ad_somatic[1], "a1" = ad_somatic[2], "total_ad" = sum(ad_somatic), "allele" = NA, "vaf" = ad_somatic[2]/sum(ad_somatic), "start" = var_start, "end" = var_end)
  if (vardf$vaf == 0 | nrow(mydf) <= 10) {
    message("Mutation not in this sample, or not enough surrounding variants")
    return(0)
  }

  # Fit model
  mylo <- loess(vaf~start, weights = total_ad, data=mydf, degree = 2)

  # Add model predictions and 95%CI
  pred <- predict(mylo, mydf$start, se = TRUE)
  mydf$predicted_vaf <- pred$fit
  mydf$predicted_vaf_se <- pred$se.fit
  mydf$lower <- mydf$predicted_vaf - 1.96*mydf$predicted_vaf_se
  mydf$lower <- mydf$predicted_vaf + 1.96*mydf$predicted_vaf_se

  mydf$dist <- abs(mydf$vaf-mydf$predicted_vaf)

  # The allele of the somatic variant is unknown.
  # Switch vaf allele if that makes it closer to the predicted vaf.
  vardf$predicted_vaf <- predict(mylo, vardf$start)
  if ( is.na( vardf$predicted_vaf ) ) {
    message("Cannot predict vaf")
    return(0)
  }
  if ( vardf$predicted_vaf < 0 ) {
    vardf$predicted_vaf <- 0
  } else if ( vardf$predicted_vaf > 1 ) {
    vardf$predicted_vaf <- 1
  }
  if ((vardf$predicted_vaf > 0.5 & vardf$vaf < 0.5) | (vardf$predicted_vaf < 0.5 & vardf$vaf > 0.5)){
    vardf$vaf = 1-vardf$vaf
    a0 = vardf$a1
    vardf$a1 = vardf$a0
    vardf$a0 = a0
  }
  vardf$dist <- abs(vardf$vaf-vardf$predicted_vaf)

  message(paste0("Created model for sample: ", samplename))

  # Create figure of fit
  # vaf_fig <- ggplot(data = mydf, aes(y = vaf, x = pos)) +
  #   geom_point(data = vardf, aes(x = pos, y = vaf, size = total_ad), colour = "blue") +
  #   geom_point(aes(size = total_ad), colour = "darkgrey") +
  #   geom_line(aes(y = predicted_vaf)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  #   guides(colour = "none") +
  #   labs(size = "Allele depth") +
  #   theme_classic()

  # Determine if mutation is different based on predicted vaf
  pval_binom <- binom.test(vardf$a1, n = vardf$total_ad, p = vardf$predicted_vaf, alternative = "two.sided")$p.value

  # Create boxplot of vaf
  # box_fig <- ggplot(data=mydf,aes(y=dist)) +
  #   geom_boxplot() +
  #   geom_point(data=vardf,aes(y=dist),x=0,col="blue",size=3) +
  #   ylim(0,1) +
  #   theme_classic() +
  #   theme_bw()

  # message(paste0("Created figure for sample: ", samplename))

  # Save figure output
  #pdf(file.path(outdir, paste0(samplename, "_", varchr, "_", varpos, ".pdf")))
  #grid.arrange(vaf_fig, box_fig, ncol = 2, nrow = 1, widths = c(4, 1))
  #dev.off()

  # message(paste0("Wrote out figure for sample: ", samplename))

  # Create output table
  mytable <- data.frame("samplename" = samplename,
                        "chrom" = varchr,
                        "start" = var_start-1,
                        "end" = var_end,
                        "mean_dist" = mean(mydf$dist),
                        "median_dist" = median(mydf$dist),
                        "sd_dist" = sd(mydf$dist),
                        "vardist" = vardf$dist,
                        "n_farther" = length(which(mydf$dist > vardf$dist)),
                        "p_farther" = length(which(mydf$dist > vardf$dist))/nrow(mydf),
                        "p_binom" = pval_binom)
  return(mytable)
}

create_empty_ab_table <- function() {
  ab_table <- data.frame("samplename" = NA,
                       "chrom" = NA,
                       "start" = NA,
                       "end" = NA,
                       "mean_dist" = NA,
                       "median_dist" = NA,
                       "sd_dist" = NA,
                       "vardist" = NA,
                       "n_farther" = NA,
                       "p_farther" = NA,
                       "p_binom" = NA)[-1,]
  return( ab_table )
}
if ( length(somatic_vcf_chr) == 0 ) {
  ab_table <- create_empty_ab_table()
} else {
  ab_table <- list()
  for ( i in c(1:length(somatic_vcf_chr))) {
    ab_table_tmp <- AB_somatic( somatic_vcf_chr[i,] )
    if ( length(ab_table_tmp) == 0 ) {
      ab_table_tmp <- create_empty_ab_table()
    } else if ( ab_table_tmp == 0 ) {
      ab_table_tmp <- create_empty_ab_table()
    } else if ( ncol(ab_table_tmp) <= 1 ) {
      ab_table_tmp <- create_empty_ab_table()
    }
    ab_table[[i]] <- ab_table_tmp
  }
  ab_table <- do.call(rbind, ab_table)
  if ( length(ab_table) == 0 ) {
    ab_table <- create_empty_ab_table()
  }
}
write.table(ab_table, output_file, quote = FALSE, col.names = TRUE, row.names = FALSE)

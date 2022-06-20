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

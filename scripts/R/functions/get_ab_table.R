get_ab_table <- function( fname ) {
  df <- read.table(fname,header=T)
  if ( nrow(df) == 0 ) {
    ab_table <- data.frame("CHROM" = NA, "START" = NA, "END" = NA, "DONOR_ID" = NA, "p_binom" = NA)[-1,]
    return( ab_table )
  }
  colnames(df)[c(2,3,4)] <- c("CHROM","START","END")

  df$DONOR_ID <- donor_id
  df$REGION <- paste(df$CHROM,paste(df$START,df$END,sep="_"),sep="_")
  # ab.data <- df %>% group_by(CHROM, END) %>% slice_sample(n = 1) %>% as.data.frame()
  ab.data <- df[unlist(lapply(by(df,df[,c("REGION")],rownames),function(x){ if (length(x) >0 ) { sample(x,1)} })),]
  ab.data <- ab.data[,c("CHROM","START","END","DONOR_ID","p_binom")]

  # Set data types correctly
  if ("CHROM" %in% colnames(ab.data)) {
    ab.data$CHROM <- as.character(ab.data$CHROM)
  }

  return( ab.data )
}

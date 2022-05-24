.df_to_bedpe <- function(df, output_file = "", score, chrom = "chrom", start.pos = "start.pos", end.pos = "end.pos"){
  bedpe <- df[,c(chrom, start.pos, end.pos, chrom, start.pos, end.pos)]
  bedpe$name <- paste(bedpe$chrom, bedpe$start.pos, bedpe$end.pos, sep = "_")
  bedpe$score <- df[,score]
  bedpe$strand1 <- "."
  bedpe$strand2 <- "."
  if(output_file != ""){
    print(paste("# Writing: ", output_file, sep = ""))
    write.table(bedpe, file = output_file, quote = F, row.names = F, col.names = F, sep = "\t")
  }
  return(bedpe)
}

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

library(randomForest)
source(paste(sourcedir,"/functions/get_indel_feature_table.R",sep=""), chdir = T)
source(paste(sourcedir,"/functions/get_ab_table.R",sep=""), chdir = T)

args <- commandArgs(trailingOnly = TRUE)
ab <- args[1]
bed <- args[2]
donor_id <- args[3]
output_file <- args[4]
muttype_levels <- c(
  "C_deletion", "T_deletion","C_insertion","T_insertion",
  "2bp_deletion","3bp_deletion","4bp_deletion","5+bp_deletion",
  "2bp_insertion","3bp_insertion","4bp_insertion","5+bp_insertion",
  "2bp_deletion_with_microhomology","3bp_deletion_with_microhomology","4bp_deletion_with_microhomology","5+bp_deletion_with_microhomology"
)

ft <- get_feature_table( bed )
at <- get_ab_table( ab )
rf_table <- merge(ft, at, by = c("DONOR_ID", "CHROM", "START", "END"), all.x = TRUE, sort = FALSE)

if ("p_binom" %in% colnames(rf_table)) {
  rf_table$p_binom <- log(as.numeric(as.vector(rf_table$p_binom)))
}

rf_table <- rf_table[,-which(colnames(rf_table) == "TSB")]


# write.table(rf_table, file=output_file, quote=F, row.names=F, col.names=T)
saveRDS(rf_table, file=output_file)

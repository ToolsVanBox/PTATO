library(randomForest)

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

get_feature_table <- function( fname ) {
  df <- read.table(fname)

  if (!file.size(fname)){
      df <- data.frame("CHROM" = NA, "START" = NA, "END" = NA, "CONTEXT" = NA, "SCORE" = NA, "STRAND" = NA, "DONOR_ID" = NA)[-1,]
      return( df )
  }

  # Select only indels
  df <- df[!grepl(">",df$V4),]
  if ( nrow(df) == 0 ) {
      df <- data.frame("CHROM" = NA, "START" = NA, "END" = NA, "CONTEXT" = NA, "SCORE" = NA, "STRAND" = NA, "DONOR_ID" = NA)[-1,]
      return( df )
  }

  # Create test data matrix
  feature.data <- matrix(unlist(strsplit(as.character(df$V8),";")),ncol=4,byrow = T)
  colnames(feature.data) <- strsplit(as.character(df$V7),";")[[1]]
  feature.data <- cbind(df[,c(1:6)],feature.data)
  colnames(feature.data)[1:6] <- c("CHROM","START","END","CONTEXT","SCORE","STRAND")

  # Set transcription strand bias
  if ("TSB" %in% colnames(feature.data)) {
    TSB2 <- rep("Unknown",nrow(feature.data))
    TSB2[feature.data$STRAND == "+" & feature.data$TSB == "+"] <- "Untranscribed"
    TSB2[feature.data$STRAND == "+" & feature.data$TSB == "-"] <- "Transcribed"
    TSB2[feature.data$STRAND == "-" & feature.data$TSB == "+"] <- "Transcribed"
    TSB2[feature.data$STRAND == "-" & feature.data$TSB == "-"] <- "Untranscribed"
    feature.data$TSB <- TSB2
  }

  # Get context
  if ("CONTEXT" %in% colnames(feature.data)) {
    feature.data <- feature.data[!grepl(",",feature.data$CONTEXT),]
    c1 <- gsub("(.+)\\(.+","\\1",feature.data$CONTEXT)
    c2 <- gsub(".+\\(.+\\)(.+)","\\1",feature.data$CONTEXT)
    muttype <- gsub(".+\\((.+_.+)-(.+)\\).+","\\1",feature.data$CONTEXT)
    muttype_sub <- gsub(".+\\((.+_.+)-(.+)\\).+","\\2",feature.data$CONTEXT)
    muttype_sub <- gsub("\\+","",muttype_sub)
    df1.context <- data.frame(do.call('rbind', strsplit(as.character(c1),'',fixed=TRUE)))
    df2.context <- data.frame(do.call('rbind', strsplit(as.character(c2),'',fixed=TRUE)))
    colnames(df1.context) <- c(paste("POSm",c(10:1),sep=""))
    colnames(df2.context) <- c(paste("POSp",c(1:10),sep=""))
    df.context <- cbind(df1.context, df2.context)
    df.context <- cbind(df.context, "MUTTYPE" = muttype)
    df.context <- cbind(df.context, "MUTTYPE_SUB" = muttype_sub)
    feature.data <- cbind(feature.data, df.context)
  }

  # Set data types correctly
  if ("CHROM" %in% colnames(feature.data)) {
    feature.data$CHROM <- as.character(feature.data$CHROM)
  }
  if ("GENEBODY" %in% colnames(feature.data)) {
    feature.data$GENEBODY <- as.numeric(as.vector(feature.data$GENEBODY))
  }
  if ("REPLISEQ" %in% colnames(feature.data)) {
    feature.data$REPLISEQ <- as.numeric(as.vector(feature.data$REPLISEQ))
  }
  if ("SIMPLEREPEAT" %in% colnames(feature.data)) {
    feature.data$SIMPLEREPEAT <- as.numeric(as.vector(feature.data$SIMPLEREPEAT))
  }
  if ("TSB" %in% colnames(feature.data)) {
    feature.data$TSB <- factor(feature.data$TSB, levels=c("Transcribed","Unknown","Untranscribed"))
  }
  if ("MUTTYPE" %in% colnames(feature.data)) {
    feature.data$MUTTYPE <- factor(feature.data$MUTTYPE, levels=muttype_levels)
  }
  if ("MUTTYPE_SUB" %in% colnames(feature.data)) {
    feature.data$MUTTYPE_SUB <- as.numeric(as.vector(feature.data$MUTTYPE_SUB))
  }
  if ("POSm1" %in% colnames(feature.data)) {
    context_cols_i = grepl("POS[mp]", colnames(feature.data))
    feature.data[,context_cols_i] = lapply(feature.data[,context_cols_i], factor, levels = c("A","C","G","T"))
  }
  feature.data$DONOR_ID <- donor_id

  return( feature.data )
}

ft <- get_feature_table( bed )
at <- get_ab_table( ab )
rf_table <- merge(ft, at, by = c("DONOR_ID", "CHROM", "START", "END"), all.x = TRUE)

if ("p_binom" %in% colnames(rf_table)) {
  rf_table$p_binom <- log(as.numeric(as.vector(rf_table$p_binom)))
}

rf_table <- rf_table[,-which(colnames(rf_table) == "TSB")]


# write.table(rf_table, file=output_file, quote=F, row.names=F, col.names=T)
saveRDS(rf_table, file=output_file)

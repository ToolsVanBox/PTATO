library(randomForest)

args <- commandArgs(trailingOnly = TRUE)
label_1 <- args[1]
files_1 <- args[2]
label_2 <- args[3]
files_2 <- args[4]
version <- args[5]

# Function to read the training data
read_train_data = function(fname, label){
  df <- readRDS(fname)
  if (nrow(df)){
    df$val <- label
  }
  return(df)
}

# Read and combine training data
train.data.1.l <-lapply(strsplit(files_1,",")[[1]], read_train_data, label_1)
train.data.2.l <-lapply(strsplit(files_2,",")[[1]], read_train_data, label_2)
train.data <- do.call(rbind, c(train.data.1.l, train.data.2.l))

# Clean up data
train.data <- train.data[!(duplicated(train.data[,c("CHROM","START","END","val")])),]
train.data <- train.data[!(duplicated(train.data[,c("CHROM","START","END")]) | duplicated(train.data[,c("CHROM","START","END")], fromLast = T)), ]

train.data.sub <- train.data[,-c(1:7)]

# Train random forest
train.data.sub.1 <- na.omit(train.data.sub)
train.data.sub.1 <- train.data.sub.1[is.finite(train.data.sub.1$p_binom),]
output_forest <- randomForest(as.factor(val) ~ ., data=train.data.sub.1, mtry=4)

# Train random forest without Allelic Imbalance
train.data.sub.2 <- train.data.sub[,-which(colnames(train.data.sub) == "p_binom")]
train.data.sub.2 <- na.omit(train.data.sub.2)
output_forest_ab <- randomForest(as.factor(val) ~ ., data=train.data.sub.2, mtry=4)

# Train random forest without Allelic Imbalance and without Repliseq
train.data.sub.3 <- train.data.sub[,-which(colnames(train.data.sub) == "REPLISEQ" | colnames(train.data.sub) == "p_binom")]
train.data.sub.3 <- na.omit(train.data.sub.3)
output_forest_repliseq <- randomForest(as.factor(val) ~ ., data=train.data.sub.3, mtry=4)

# Write out confusion and importance tables
write.table(output_forest$confusion,file = paste0(paste0("randomforest", version, "_confusion.txt")), quote = F )
write.table(importance(output_forest),file = paste0(paste0("randomforest", version, "_importance.txt")), quote = F)

write.table(output_forest_ab$confusion,file = paste0(paste0("randomforest_ab", version, "_confusion.txt")), quote = F )
write.table(importance(output_forest_ab),file = paste0(paste0("randomforest_ab", version, "_importance.txt")), quote = F)

write.table(output_forest_repliseq$confusion,file = paste0(paste0("randomforest_repliseq", version, "_confusion.txt")), quote = F )
write.table(importance(output_forest_repliseq),file = paste0(paste0("randomforest_repliseq", version, "_importance.txt")), quote = F)

# Save the random forests.
output_forest_list = list(all = output_forest, ab = output_forest_ab, repliseq = output_forest_repliseq)
save(output_forest_list, file = paste0(paste0("randomforest", version, ".Rdata")))
saveRDS(output_forest_list, file = paste0(paste0("randomforest", version, ".rds")))

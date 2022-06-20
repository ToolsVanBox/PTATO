library(randomForest)

args <- commandArgs(trailingOnly = TRUE)
label_1 <- args[1]
files_1 <- args[2]
label_2 <- args[3]
files_2 <- args[4]
version <- args[5]

for (fname in strsplit(files_1,",")[[1]]) {
  # df <- read.table(fname,header=T)
  df <- readRDS(fname)
  df$val <- label_1
  if ( !exists("train.data") ) {
    train.data <<- df
  } else {
    train.data <- rbind(train.data, df)
  }
}

for (fname in strsplit(files_2,",")[[1]]) {
  # df <- read.table(fname,header=T)
  df <- readRDS(fname)
  df$val <- label_2
  if ( !exists("train.data") ) {
    train.data <<- df
  } else {
    train.data <- rbind(train.data, df)
  }
}

train.data <- train.data[!(duplicated(train.data[,c("CHROM","START","END","val")])),]
train.data <- train.data[!(duplicated(train.data[,c("CHROM","START","END")]) | duplicated(train.data[,c("CHROM","START","END")], fromLast = T)), ]

train.data.sub <- train.data[,-c(1:7)]

train.data.sub.1 <- na.omit(train.data.sub)
train.data.sub.1 <- train.data.sub.1[is.finite(train.data.sub.1$p_binom),]
output_forest <- randomForest(as.factor(val) ~ ., data=train.data.sub.1, mtry=4)

train.data.sub.2 <- train.data.sub[,-which(colnames(train.data.sub) == "p_binom")]
train.data.sub.2 <- na.omit(train.data.sub.2)
output_forest_ab <- randomForest(as.factor(val) ~ ., data=train.data.sub.2, mtry=4)

train.data.sub.3 <- train.data.sub[,-which(colnames(train.data.sub) == "REPLISEQ" | colnames(train.data.sub) == "p_binom")]
train.data.sub.3 <- na.omit(train.data.sub.3)
output_forest_repliseq <- randomForest(as.factor(val) ~ ., data=train.data.sub.3, mtry=4)

write.table(output_forest$confusion,file = paste0(paste0("randomforest", version, "_confusion.txt")), quote = F )
write.table(importance(output_forest),file = paste0(paste0("randomforest", version, "_importance.txt")), quote = F)

write.table(output_forest_ab$confusion,file = paste0(paste0("randomforest_ab", version, "_confusion.txt")), quote = F )
write.table(importance(output_forest_ab),file = paste0(paste0("randomforest_ab", version, "_importance.txt")), quote = F)

write.table(output_forest_repliseq$confusion,file = paste0(paste0("randomforest_repliseq", version, "_confusion.txt")), quote = F )
write.table(importance(output_forest_repliseq),file = paste0(paste0("randomforest_repliseq", version, "_importance.txt")), quote = F)

output_forest_list = list(all = output_forest, ab = output_forest_ab, repliseq = output_forest_repliseq)
save(output_forest_list, file = paste0(paste0("randomforest", version, ".Rdata")))
saveRDS(output_forest_list, file = paste0(paste0("randomforest", version, ".rds")))

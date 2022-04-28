library(randomForest)
library(VariantAnnotation)

args <- commandArgs(trailingOnly = TRUE)

rf <- readRDS( args[1] )
test.data <- readRDS( args[2] )
# colnames(test.data)[colnames(test.data) == "TSB"] <- "TSB2"
# colnames(test.data)[colnames(test.data) == "p_binom"] <- "rand_p_binom"
myvcf <- readVcf( args[3] )
out_file <- args[4]

if (!nrow(test.data) & !length(myvcf)){

    meta(header(myvcf))$PTArandomforest <- DataFrame("Value" = paste0('"', args[1], ' ', '"'), row.names = "PTArandomforest")

    infoheader <- rbind(info(header(myvcf)), "PTAprob" = DataFrame("Number" = "1", "Type" = "Float", "Description" = "PTA probability value"))
    rownames(infoheader) <- c(rownames(info(header(myvcf))),"PTAprob")
    info(header(myvcf)) <- infoheader
    info(myvcf)$PTAprob <- numeric()

    outvcf <- file(out_file, open="a")
    writeVcf(myvcf, outvcf)
} else{

    test.data[which(is.infinite(test.data$p_binom)),"p_binom"] <- NA

    ptable <- predict(rf$all, test.data, type='prob')
    test.data$PTAprob <- ptable[,"PTA"]
    test.data$PTAprob2 <- predict(rf$ab, test.data, type="prob")[,"PTA"]
    test.data$PTAprob3 <- predict(rf$repliseq, test.data, type="prob")[,"PTA"]

    # test.data[is.na(test.data$PTAprob),"PTAprob"] <- predict(rf$ab, test.data[is.na(test.data$PTAprob),], type="prob")[,"PTA"]
    # test.data[is.na(test.data$PTAprob),"PTAprob"] <- predict(rf$repliseq, test.data[is.na(test.data$PTAprob),], type="prob")[,"PTA"]

    gr <- granges(myvcf)
    gr.df <- as.data.frame(gr)
    colnames(gr.df)[c(1,3)] <- c("CHROM","END")
    gr.df.merged <- merge(gr.df, test.data, by = c("CHROM", "END"), all.x = TRUE)

    # gr$PTAprob <- gr.df.merged$PTAprob
    gr$PTAprob <- paste(gr.df.merged$PTAprob,gr.df.merged$PTAprob2,gr.df.merged$PTAprob3,sep=",")

    meta(header(myvcf))$PTArandomforest <- DataFrame("Value" = paste0('"', args[1], ' ', '"'), row.names = "PTArandomforest")

    infoheader <- rbind(info(header(myvcf)), "PTAprob" = DataFrame("Number" = "3", "Type" = "Float", "Description" = "PTA probability value"))
    rownames(infoheader) <- c(rownames(info(header(myvcf))),"PTAprob")
    info(header(myvcf)) <- infoheader
    info(myvcf)$PTAprob <- gr$PTAprob

    outvcf <- file(out_file, open="a")
    writeVcf(myvcf, outvcf)
}

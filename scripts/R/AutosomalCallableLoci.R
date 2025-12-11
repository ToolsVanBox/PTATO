

# get input arguments
args <- commandArgs(trailingOnly = TRUE)
callableloci_file <- strsplit(args[1],",")[[1]]
callability_outputfile <- strsplit(args[2],",")[[1]]

## 1. Callable loci
# In this step, a .txt file is generated showing the amount and frequencies of bases that are callable. 
# In contrast to the standard CallableLoci output, this txt file only contains data for the autosomes, to enable comparisons between genders
# This step takes quite a while, because all callableloci.bed files are read 
CHROMOSOMES <- c(1:22, "X", "Y")
AUTOSOMES <- c(1:22)
callableloci_overview <- data.frame()

callable_bed <- read.delim(callableloci_file, header = F)

# Remove the decoy chromosomes:
callableloci  <- callable_bed[which(callable_bed[,1] %in% CHROMOSOMES),]

# Make a summary file (autosomes-only)
callable_autosomes <- callableloci[which(callableloci[,1] %in% AUTOSOMES),]
callable_autosomes$Width <- callable_autosomes[,3] - callable_autosomes[,2]
callability <- aggregate(callable_autosomes$Width, by=list(State=callable_autosomes[,4]), FUN=sum)
names(callability)[2] <- "nBases"

# The "REF_N" regions are not taken into account for CallableLoci
callability$Freq <- callability[,2] / sum(callability$nBases[callability$State != "REF_N"])

# Write the output file 
print(paste("# Writing: ", callability_outputfile, sep = ""))
write.table(x = callability, file = callability_outputfile, sep = "\t", quote = F, row.names = F)

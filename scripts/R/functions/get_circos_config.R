### This script contains functions to make two Circos plots showing copy number changes and structural variant breakpoint junctions
# The first circos plot uses a scatter plot to show read depth and BAF (in 1mb bins)
# The second plot uses segments (histogram) to show read depth and BAF segments

## To DO
# visual improvements

Make_SV_Circos <- function(CNVs_file, ReadCounts_1mb_file, BAF_1mb_file, ReadCounts_segments_file, BAF_segments_file, SVs_file, OUTPUT_dir, Circos_Path ){

  # Input_Dir should contain the output of the PTATO-SV pipeline
  CNVs <- read.delim(CNVs_file)
  ReadCounts_1mb <- read.delim(ReadCounts_1mb_file)
  BAF_1mb <- read.delim(BAF_1mb_file)

  # ReadCounts (scatter)
  ReadCounts_Circos <- Dataframe_to_Circos(ReadCounts_1mb, Value_column = 5)
  ReadCounts_Circos_a <- Overlap_CNVs(Data = ReadCounts_Circos, CNVs =  CNVs[,c(1,2,3,6)])
  print(paste("# Writing: ", OUTPUT_dir, ".circos.readcounts.1mb.txt", sep = ""))
  write.table(ReadCounts_Circos_a,
              file = paste(OUTPUT_dir, ".circos.readcounts.1mb.txt", sep = ""), quote = F, row.names = F, col.names = F,  sep = "\t")

  # BAF (scatter)
  BAF_Circos <- Dataframe_to_Circos(BAF_1mb, Chrom_column = 1, Start_column = 2,End_column = 3, Value_column = 4)
  BAF_Circos_a <- Overlap_CNVs(Data = BAF_Circos, CNVs = CNVs[,c(1,2,3,6)])
  print(paste("# Writing: ", OUTPUT_dir, ".circos.baf.txt", sep = ""))
  write.table(BAF_Circos_a,
              file = paste(OUTPUT_dir, ".circos.baf.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")

  # Links (scatter + segments)
  SVs <- readVcf( SVs_file )
  Links <- SV_VCF_to_Links(SVs = SVs)
  print(paste("# Writing: ", OUTPUT_dir, ".circos.links.txt", sep = ""))
  write.table(Links,
              file = paste(OUTPUT_dir, ".circos.links.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")

  # ReadCount Segments
  ReadCount_Segments <- read.delim(ReadCounts_segments_file)
  ReadCount_Segments_Circos <- Dataframe_to_Circos(ReadCount_Segments, Chrom_column = 1, Start_column = 2,End_column = 3, Value_column = 4)
  ReadCount_Segments_Circos[,4] <- ReadCount_Segments_Circos[,4] - 2
  ReadCount_Segments_Circos$Param <- "color=dgrey,fill_color=dgrey"
  ReadCount_Segments_Circos$Param[ReadCount_Segments_Circos$mean < -0.5] <- "color=dred,fill_color=dred"
  ReadCount_Segments_Circos$Param[ReadCount_Segments_Circos$mean > 1.5] <- "color=dblue,fill_color=dblue"
  write.table(ReadCount_Segments_Circos,
              file = paste(OUTPUT_dir, ".circos.readcounts.segments.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")

  # BAF segments
  BAF_Segments_Input <- read.delim(BAF_segments_file)
  BAF_Segments <- Dataframe_to_Circos(BAF_Segments_Input, Chrom_column = 1, Start_column = 2,End_column = 3, Value_column = 4)
  BAF_Segments$Param <- "color=dgrey,fill_color=vlgrey"
  BAF_Segments$Param[BAF_Segments[,4] > 0.16] <- "color=dblue,fill_color=lblue"
  BAF_Segments$Param[BAF_Segments[,4] > 0.4] <- "color=dred,fill_color=lred"
  write.table(BAF_Segments,
              file = paste(OUTPUT_dir, ".circos.baf.segments.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")

  Scatter_Conf <- Generate_Scatter_CONF(Circos_Path = Circos_Path,
                        Circos_Input_DIR = OUTPUT_dir)

  print("# Running Circos")
  circos_scatter_command <- paste("circos -conf ", OUTPUT_dir, "Circos_Scatter.conf -outputdir ",OUTPUT_dir, " -outputfile ", OUTPUT_dir, ".circos.scatter", sep = "")
  print(circos_scatter_command)
  # system(circos_scatter_command)

  Segments_Conf <- Generate_Segment_CONF(Circos_Path = Circos_Path,
                                        Circos_Input_DIR = OUTPUT_dir)

  circos_segments_command <- paste("circos -conf ", OUTPUT_dir, "Circos_Segments.conf -outputdir ",OUTPUT_dir, " -outputfile ", OUTPUT_dir, ".circos.segments", sep = "")
  print(circos_segments_command)
  # system(circos_segments_command)
}

Dataframe_to_Circos <- function(Input, Chrom_column = 1, Start_column = 2, End_column = 3, Value_column = 4, Downsample = 0){
  Output <- Input[,c(Chrom_column, Start_column, End_column, Value_column)]
  Output[,Chrom_column] <- paste("hs", Output[,Chrom_column], sep = "")
  if(Downsample > 0){
    Output <- Output[seq(1, nrow(Output), Downsample),]
  }
  return(Output)
}

Overlap_CNVs <- function(CNVs, Data, Gain = "dblue", Loss = "dred", LOH = "dorange", Normal = "dgrey_a3"){
  Output <- Data
  Data_g <- GRanges(seqnames = Data[,1], IRanges(Data[,2], Data[,3]))
  Output$Param <- paste("color=", Normal, sep = "")
  CNVs_g <- GRanges(seqnames = paste("hs", CNVs[,1], sep = ""), IRanges(CNVs[,2], CNVs[,3]), State = CNVs[,4])

  Olap_Gains<- findOverlaps(Data_g, CNVs_g[CNVs_g$State == "Gain",])
  Output$Param[queryHits(Olap_Gains)] <- paste("color=", Gain, sep = "")
  Olap_Losses <- findOverlaps(Data_g, CNVs_g[CNVs_g$State == "Loss",])
  Output$Param[queryHits(Olap_Losses)] <- paste("color=", Loss, sep = "")
  Olap_LOH <- findOverlaps(Data_g, CNVs_g[CNVs_g$State == "LOH",])
  Output$Param[queryHits(Olap_LOH)] <- paste("color=", LOH, sep = "")

  return(Output)
}

SV_VCF_to_Links <- function(SVs){

  Links <- data.frame()
  Events <- unique(info(SVs)$EVENT)
  for(Event in Events){

    SV <- SVs[info(SVs)$EVENT == Event,]
    if(length(SV) == 2){

      Link <- data.frame(chrom1 = paste("hs", as.vector(seqnames(SV))[1], sep = ""),
                         start1 = start(rowRanges(SV)[1]),
                         end1 = end(rowRanges(SV)[1]),
                         chrom2 = paste("hs", as.vector(seqnames(SV))[2], sep = ""),
                         start2 = start(rowRanges(SV)[2]),
                         end2 = end(rowRanges(SV)[2]))

      Link$param <- ifelse(unique(info(SV)$SVTYPE) == "DEL", "color=dred", "color=dblue")
      Link$param <- ifelse(unique(info(SV)$SVTYPE) == "CTX", "color=dpurple", Link$param)
      Link$param <- ifelse(unique(info(SV)$SVTYPE) == "INV", "color=dgreen", Link$param)
    }
    Links <- rbind(Links, Link)
  }

  return(Links)
}


Generate_Scatter_CONF <- function(Circos_Path, Circos_Input_DIR){
  CONF <- paste(
  "<<include ",Circos_Path,"etc/colors_fonts_patterns.conf>>

show_ticks = no
show_tick_labels = no

<ideogram>
<spacing>
default = 0.002r
</spacing>

# Ideogram position, fill and outline
radius = 0.90r
thickness = 25p
fill = yes
stroke_color = dgrey
stroke_thickness = 2p

# Minimum definition for ideogram labels.

show_label       = yes
# see etc/fonts.conf for list of font names
label_font       = default
label_radius     = 1r + 30p
label_size       = 30
label_parallel   = yes
label_case       = upper

</ideogram>

# Decrease the resolution of the output file to improve readability
<image>
<<include ",Circos_Path,"etc/image.conf>>
radius* = 1000p
</image>

karyotype = ",Circos_Path,"data/karyotype/karyotype.human.hg38.txt

chromosomes_units = 1000000
chromosomes_display_default = yes
chromosomes = -hsZ

<plots>

## Somatic SNV Rainfall Scatter Plot
# This data is in log10 scale
<plot>
type = scatter
file = ",Circos_Input_DIR, ".circos.readcounts.1mb.txt
r0   = 0.7r
r1   = 0.975r
min = 0
max = 4
glyph = circle
glyph_size = 15

<axes>
<axis>
color     = lpurple
thickness = 1
spacing   = 0.25r
</axis>
</axes>

<backgrounds>
<background>
color = vlred_a5
y0 = 0
y1 = 2
</background>
<background>
color = vlblue_a5
y0 = 2
y1 = 4
</background>
</backgrounds>
</plot>

<plot>
type = scatter
file = ",Circos_Input_DIR, ".circos.baf.txt
r0   = 0.375r
r1   = 0.65r
min = 0
max = 0.5
glyph = circle
glyph_size = 15

<axes>
<axis>
color     = lpurple
thickness = 1
spacing   = 0.2r
</axis>
</axes>

<backgrounds>
<background>
color = vlgrey_a5
y0 = 0
y1 = 0.2
</background>
<background>
color = vlblue_a5
y0 = 0.2
y1 = 0.4
</background>
<background>
color = vlred_a5
y0 = 0.4
y1 = 0.5
</background>
</backgrounds>

</plot>

</plots>

<links>
<link>
file          = ",Circos_Input_DIR,".circos.links.txt
radius        = 0.35r
bezier_radius = 0.1r
thickness     = 10
ribbon        = yes
</link>
</links>

<<include ",Circos_Path,"etc/housekeeping.conf>>

<colors>
chr1 = 128,125,186
chr2 = 145,142,179
chr3 = 161,159,173
chr4 = 179,176,166
chr5 = 196,193,160
chr6 = 213,210,153

chr7 = 230,228,147
chr8 = 202,218,138
chr9 = 175,209,129
chr10 = 147,199,120
chr11 = 120,190,111
chr12 = 92,180,102

chr13 = 65,171,93
chr14 = 65,166,110
chr15 = 65,162,128
chr16 = 65,158,145
chr17 = 65,154,163
chr18 = 65,150,180

chr19 = 66,146,198
chr20 = 76,142,196
chr21 = 86,139,194
chr22 = 97,135,192
chrX = 107,132,190
chrY = 117,128,188
</colors>
", sep = "")
  print(paste(Circos_Input_DIR, ".circos.scatter.conf", sep = ""))
  write.table(CONF, paste(Circos_Input_DIR, ".circos.scatter.conf", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
  return(CONF)
}


Generate_Segment_CONF <- function(Circos_Path, Circos_Input_DIR){
  CONF <- paste(
    "<<include ",Circos_Path,"etc/colors_fonts_patterns.conf>>

show_ticks = no
show_tick_labels = no

<ideogram>
<spacing>
default = 0.002r
</spacing>

# Ideogram position, fill and outline
radius = 0.90r
thickness = 25p
fill = yes
stroke_color = dgrey
stroke_thickness = 2p

# Minimum definition for ideogram labels.

show_label       = yes
# see etc/fonts.conf for list of font names
label_font       = default
label_radius     = 1r + 30p
label_size       = 30
label_parallel   = yes
label_case       = upper

</ideogram>

# Decrease the resolution of the output file to improve readability
<image>
<<include ",Circos_Path,"etc/image.conf>>
radius* = 1000p
</image>

karyotype = ",Circos_Path,"data/karyotype/karyotype.human.hg38.txt

chromosomes_units = 1000000
chromosomes_display_default = yes
chromosomes = -hsZ

<plots>

<plot>
type = histogram
file = ",Circos_Input_DIR, ".circos.readcounts.segments.txt
r0   = 0.7r
r1   = 0.975r
min = -2
max = 2
stroke_type = outline
thickness   = 4
extend_bin  = no

<axes>
<axis>
color     = lgrey
thickness = 1
spacing   = 1
</axis>
<axis>
color     = dgrey
thickness = 3
spacing   = 2
</axis>
</axes>

<backgrounds>
<background>
color = vlred_a5
y0 = -2
y1 = 0
</background>
<background>
color = vlblue_a5
y0 = 0
y1 = 2
</background>
</backgrounds>
</plot>

<plot>
type = histogram
file = ",Circos_Input_DIR, ".circos.baf.segments.txt
r0   = 0.375r
r1   = 0.65r
min = 0
max = 0.5
stroke_type = outline
thickness   = 4
extend_bin  = no

<axes>
<axis>
color     = lgrey
thickness = 1
spacing   = 0.2r
</axis>
<axis>
color     = dgrey
thickness = 3
spacing   = 1r
</axis>
</axes>

<backgrounds>
<background>
color = vlgrey_a5
y0 = 0
y1 = 0.2
</background>
<background>
color = vlblue_a5
y0 = 0.2
y1 = 0.4
</background>
<background>
color = vlred_a5
y0 = 0.4
y1 = 0.5
</background>
</backgrounds>

</plot>

</plots>

<links>
<link>
file          = ",Circos_Input_DIR,".circos.links.txt
radius        = 0.35r
bezier_radius = 0.1r
thickness     = 10
ribbon        = yes
</link>
</links>

<<include ",Circos_Path,"etc/housekeeping.conf>>

<colors>
chr1 = 128,125,186
chr2 = 145,142,179
chr3 = 161,159,173
chr4 = 179,176,166
chr5 = 196,193,160
chr6 = 213,210,153

chr7 = 230,228,147
chr8 = 202,218,138
chr9 = 175,209,129
chr10 = 147,199,120
chr11 = 120,190,111
chr12 = 92,180,102

chr13 = 65,171,93
chr14 = 65,166,110
chr15 = 65,162,128
chr16 = 65,158,145
chr17 = 65,154,163
chr18 = 65,150,180

chr19 = 66,146,198
chr20 = 76,142,196
chr21 = 86,139,194
chr22 = 97,135,192
chrX = 107,132,190
chrY = 117,128,188
</colors>
", sep = "")
  print(paste(Circos_Input_DIR, ".circos.segments.conf", sep = ""))
  write.table(CONF, paste(Circos_Input_DIR, ".circos.segments.conf", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
  return(CONF)
}

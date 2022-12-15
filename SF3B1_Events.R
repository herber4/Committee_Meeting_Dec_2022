library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggcoverage)
library(rtracklayer)
library(graphics)
library(edgeR)
library(pheatmap)
library(ggbreak)
library(ggpubr)
library(patchwork)


#read in encode event frequencies
enc_freqs <- read.table(file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/enc_event_freqs.txt", sep = "\t", header = T)
brca <- read.table(file = "/Users/herber4/Desktop/brca_k700e_FDR.txt", header = T, sep = "\t")

brca %>%
  filter(deltapsi > 0) %>%
  count()
brca$Sample <- "BRCA"



#
#
#find shared C3SS in K700E mutant samples
#read in leafcutter output format long from Nalm-6, K562, and MDS patient data 
df <- read.table(file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/k700e_shared_events.txt", header = T, sep = "\t")
df <- rbind(df, brca)
write.table(df, file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/all_shared_events.txt", sep = "\t")

pos <- filter(df, deltapsi > 0)

shared <- df %>%
  filter(deltapsi > 0) %>%
  group_by(chr, start, end) %>%
  count()

shared <- shared %>%
  filter(n >= 2)

dpsi <- merge(shared, pos)
dpsi$end_minus <- dpsi$end-400
dpsi$end_plus <- dpsi$end+400
dpsi$joined_coords <- str_c(dpsi$end_minus, '-', dpsi$end_plus)
dpsi$SS_Coords <- str_c(dpsi$chr, ':', dpsi$joined_coords)
dpsi$joined_coords <- NULL
dpsi$end_plus <- NULL
dpsi$end_minus <- NULL
dpsi$k_end <- NULL
dpsi$k_start <- NULL
write.table(dpsi, file = "/Users/herber4/Desktop/SF3B1/all_k7_events_master.txt", sep = "\t")



genes <- read.table(file = "/Users/herber4/Desktop/SF3B1/BRCA/shared_splice_genes.txt", sep = "\t", header = T)
genes <- merge(genes, cpms, by.x = c("genes"), by.y = c("Gene.Symbol"))
new <- genes
new$sd <- NULL
new$ENSG.ID <- NULL
wide <- new %>%
  pivot_wider(names_from = Samples, values_from = median.cpm)
wide[is.na(wide)] <- 0
wide <- as.data.frame(wide)
rownames(wide) <- wide$genes
wide$genes <- NULL
wide <- log(wide, 10)

dpsi_cpms %>%
  filter(n >= 2) %>%
  ggplot(aes(x = Sample, y = median.cpm, fill = Sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 0.5, vjust = 0.5, size = 10)) +
  scale_y_break(c(1000, 1950)) +
  geom_hline(yintercept = median(genes$median.cpm), color = "black") +
  scale_fill_brewer(palette = "Pastel1")

ggplot(dpsi_cpms, aes(y = median.cpm, x = n, group = n)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank())
ggplot(filter(dpsi_cpms, n == 3), aes(x = n, y = median.cpm, group = n)) +
  geom_boxplot()

ggplot(dpsi_cpms, aes(x = Sample, y = median.cpm, fill = Sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 0.5, vjust = 0.5, size = 10)) +
  scale_y_break(c(1000, 1950)) +
  geom_hline(yintercept = median(genes$median.cpm), color = "black")

ggplot(dpsi_cpms, aes(x = verdict)) +
  geom_bar()

ggplot(dpsi_cpms, aes(x = Sample, y = deltapsi, fill = Sample)) +
  geom_violin() +
  theme_classic() +
  geom_hline(yintercept = median(dpsi_cpms$deltapsi), color = "black") +
  theme(axis.text.x = element_text(angle = 75, hjust = 0.5, vjust = 0.5, size = 10))

dpsi_cpms %>%
  distinct(verdict, event)%>%
  ggplot(aes(x = verdict, fill = verdict)) +
  geom_bar() +
  scale_fill_brewer(palette = "Pastel1") +
  theme_classic() +
  coord_flip()
  
  ggarrange(cpm_plot, dpsi_plot, ncol = 2, nrow = 1, align = c("h"), legend = FALSE)



ggplot(dpsi_cpms, aes(x = Sample, y = deltapsi, fill = Sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 0.5, vjust = 0.5, size = 10))



test_merge <- merge(all, genes, by.x = c("gene.y", "Sample"), by.y = c("genes", "Samples"))
test_merge <- dplyr::left_join(genes, all, by.x)
write.table(test_merge, file = "/Users/herber4/Desktop/SF3B1/dpsi_cpms_with_zeros.txt", sep = "\t")
dpsi_cpms <- test_merge

ggplot(dpsi_cpms, aes(x = Sample, y = deltapsi, fill = Sample)) +
  geom_violin() +
  theme_classic()

cpm_mat <- cpms
cpm_mat <- cpms[!(cpms$ENSG.ID == "ENSG00000124333.16" | cpms$ENSG.ID == "ENSG00000002586.20"),]


ENSG00000124333.16
ENSG00000002586.20
cpm_mat$Gene.Symbol <- NULL
cpm_mat$sd <- NULL
cpm_mat <- cpm_mat %>%
  pivot_wider(names_from = Samples, values_from = median.cpm)

cpm_mat <- as.data.frame(cpm_mat)
rownames(cpm_mat) <- cpm_mat$ENSG.ID
cpm_mat$ENSG.ID <- NULL
cpm_mat[is.na(cpm_mat)] <- 0

tmat <- t(cpm_mat)
tvar <- apply(tmat, MARGIN = 2, FUN = var)
tsorted <- sort(tvar, decreasing = TRUE, index.return = TRUE)$ix[1:1000]
tdat.highvariance <- tmat[, tsorted]
tdat.highvariance
log <- -log(tdat.highvariance, 10)
pheatmap(tdat.highvariance, cluster_rows = F, cluster_cols = T, clustering_method = "complete", show_rownames = T, show_colnames = F, fontsize_col = 3,
         color = colorRampPalette(c("red", "black", "green"))(100),
         clustering_distance_cols = "euclidean",
         scale = "column")

pheatmap(cpm_mat, cluster_rows = T, cluster_cols = T, clustering_method = "complete", show_colnames = T, show_rownames = F)



c3SS <- dpsi %>%
  filter(verdict == "cryptic_threeprime")


#write out cryptic 3'SS
write.table(c3SS, file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/cryptic_3SS.txt", sep = "\t")
#
#
#

k7_introns <- read.table(file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/k700e_ints.txt", sep = "\t", header = T)
ints_pubbed <- merge(pos, k7_introns, by = c("gene.y"))
write.table(ints_pubbed, file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/ints_pubbed_dpsi.txt", sep = "\t")


#
#
#assess C3SS freq in encode db
#read in encode event frequencies
enc_freqs <- read.table(file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/enc_event_freqs.txt", sep = "\t", header = T)
enc_sums <- enc_freqs %>%
  group_by(event) %>%
  summarise(sum = sum(n))

ggplot(enc_sums, aes(x = sum)) +
  geom_histogram(binwidth = 100)
#these c3SS events are completely unique to SF3B1 K700E and do not appear in ENCODE
c3SS_freqs <- merge(dpsi_cpms, enc_sums, by = c("event"))
c3SS_freqs <- dplyr::left_join(dpsi_cpms, enc_sums, by = "event")
c3SS_freqs[is.na(c3SS_freqs)] <- 0


write.table(c3SS_freqs, file = "/Users/herber4/Desktop/SF3B1/shared_event_CPM_enc_freq.txt", sep = "\t")

uniq <- c3SS_freqs %>%
  group_by(event, sum) %>%
  distinct(event)

ggplot(uniq, aes(x = sum)) +
  geom_histogram(binwidth = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank())



ggplot(uniq, aes(reorder(event, -sum), y = sum)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  xlab("") +
  ylab("")
  
#TRIP12 shows up one time
#
#
#

#cryptic splice site identifier 
gtf <- rtracklayer::import("/Users/herber4/gencode_annos/gencode.v39.annotation.gtf")
gtf_file <- as.data.frame(gtf)

#' Annotate alternative splicing events
#' 
#' This function takes a data frame containing coordinate and orientation information
#' for splice junctions, and returns a data frame with additional columns indicating
#' the type of alternative splicing event that occurs at each junction.
#' 
#' @param junctionCoords a data frame with columns for chromosome, start position,
#'   end position, clusterID, and orientation of each splice junction
#' @return a data frame with the same columns as junctionCoords, plus additional
#'   columns for the type of alternative splicing event at each junction
#'
annotateAltSplice <- function(junctionCoords) {
  # Extract the strand information from the clusterID column
  junctionCoords$orientation <- gsub("clu_[0-9]+_", "", junctionCoords$clusterID)
  
  # Correct the start and end positions for each junction based on its orientation
  junctionCoords$start <- ifelse(junctionCoords$orientation == "+", junctionCoords$start, junctionCoords$end)
  junctionCoords$end <- ifelse(junctionCoords$orientation == "+", junctionCoords$end, junctionCoords$start)
  
  # Add a column for the end position of the previous junction
  junctionCoords$prevEnd <- c(NA, junctionCoords$end[-nrow(junctionCoords)])
  
  # Add a column for the type of alternative splicing event
  junctionCoords$altSpliceType <- "none"
  
  # Identify alternative splicing events
  junctionCoords$altSpliceType[junctionCoords$start != junctionCoords$prevEnd] <- "exon skip"
  junctionCoords$altSpliceType[junctionCoords$start < junctionCoords$prevEnd & junctionCoords$strand == "+"] <- "alternative 3' splice site"
  junctionCoords$altSpliceType[junctionCoords$start > junctionCoords$prevEnd & junctionCoords$strand == "-"] <- "alternative 5' splice site"
  junctionCoords$altSpliceType[junctionCoords$start > junctionCoords$prevEnd] <- "intron retention"
  
  # Return the annotated data frame
  return(junctionCoords)
}


junctionCoords <- dpsi_cpms
junctionCoords$strand <- gsub("clu_[0-9]+_", "", junctionCoords$clusterID)

#ggcoverage for looking at coverage over splice junctions
gtf.gr <- rtracklayer::import.gff(con = "/Users/herber4/gencode_annos/gencode.v29.annotation.gtf", format = 'gtf')
track.df <- LoadTrackFile(track.folder = "/Users/herber4/Desktop/SF3B1/", format = "bw")

form <- FormatTrack(track.df, 
                    region = "chr5:864581-870245",
                    gtf.gr = gtf.gr,
                    gene.name = "CRNDE")

cov <- ggplot() +
  geom_coverage(form, color = "auto")
cov

cov +
  geom_gene(gtf.gr = gtf.gr)

basic.coverage <- ggcoverage(data = track.df,
                             region = "chr5:864581-870245",
                             range.position = "out",
                             gtf.gr = gtf.gr,
                             gene.name = "CRNDE",
                             gene.name.type = c("gene_name"),
                             group.key = "Group")
basic.coverage 
basic.coverage +
  geom_gene(gtf.gr = gtf.gr)




uniq <- shared %>%
  group_by(Sample) %>%
  subset(deltapsi > 0) %>%
  distinct(gene.y, chr, start, end, verdict)


ggplot(uniq, aes(x = verdict, group = Sample, fill = Sample)) +
  geom_bar(position = "dodge") +
  theme_classic()



samtools view -L regions.bed -o C_k5_regions.bam C_merged_K562.bam

samtools index X.bam

bamCoverage -b test.bam -o k5_k7_regions.bw -p 16



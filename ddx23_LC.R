library(ggplot2)
library(dplyr)
library(stringr)

gtex_ggc_merged <- read.table(file = "/Users/herber4/Desktop/DDX23_Master/DEG_FINAL/gtex_ggc_medians_merged_long_MASTER.txt", sep = "\t", header = T)

p1 <- merge(clusters, introns, by = "clusterID")
p2 <- merge(clusters, introns, by = "clusterID")

p1$Sample <- "P1"
p2$Sample <- "P2"
ddx_combo <- rbind(p1, p2)

ddx_combo %>%
  group_by(Sample) %>%
  filter(deltapsi > 0) %>%
  count() %>%
  ggplot(aes(x = Sample, y = n, fill = Sample)) +
  geom_col() +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  ylab("# AS Events")

ddx_combo$Sample[ddx_combo$Sample == "P1"] <- "DDX23P1"
ddx_combo$Sample[ddx_combo$Sample == "P2"] <- "DDX23P2"
ddx_dpsi_cpms <- merge(ddx_combo, gtex_ggc_merged, by.x = c("Sample", "gene.y"), by.y = c("Sample", "gene_name"))
write.table(ddx_dpsi, file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/ddx23_dpsi.txt", sep = "\t")

ddx_dpsi_cpms %>%
  group_by(gene.y) %>%
  filter(deltapsi > 0) %>%
  ggplot(aes(x = deltapsi, y = median)) +
  geom_point(size = 3) +
  theme_bw() +
  geom_hline(yintercept = median(ddx_dpsi_cpms$median),
             color = "orange",
             linetype = "dotdash") +
  geom_vline(xintercept = median(ddx_dpsi_cpms$deltapsi),
             color = "purple",
             linetype = "dotdash") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  ylab("median CPM (Expression)")


ddx_dpsi_cpms %>%
  filter(deltapsi > 0) %>%
  group_by(Sample) %>%
  ggplot(aes(x = Sample, y = median, fill = Sample)) +
  geom_boxplot() +
  theme_classic()

ddx_GO <- read.table(file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/ddx_splice_GO.txt", sep = "\t", header = T)
ggplot(ddx_GO, aes(reorder(Pathway, Fold.Enrichment), y = Fold.Enrichment, fill = nGenes)) +
  geom_col() +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  coord_flip() +
  xlab("") +
  ylab("Enrichment")




d_shared <- ddx_combo %>%
  filter(deltapsi > 0) %>%
  group_by(gene.y, chr, start, end) %>%
  count()
d_shared <- filter(d_shared, n == 2)

ddx_dpsi <- merge(ddx_combo, d_shared)
tmp <- merge(ddx_dpsi, gtex_ggc_merged, by.x = c("Sample", "gene.y"), by.y = c("Sample", "gene_name"))

tmp %>%
  group_by(Sample) %>%
  ggplot(aes(x = Sample, y = median, fill = Sample)) +
  geom_boxplot() +
  theme_classic()

tmp %>%
  group_by(gene.y) %>%
  filter(deltapsi > 0) %>%
  ggplot(aes(x = deltapsi, y = median)) +
  geom_point(size = 3) +
  theme_bw() +
  geom_hline(yintercept = median(tmp$median),
             color = "orange",
             linetype = "dotdash") +
  geom_vline(xintercept = median(tmp$deltapsi),
             color = "purple",
             linetype = "dotdash") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  ylab("median CPM (Expression)")


#rmats

dmats <- read.table(file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/ddx_rmats.txt", sep = "\t", header = T)

dmats %>%
  group_by(Patient) %>%
  count() %>%
  ggplot(aes(x = Patient, y = n, fill = Patient)) +
  geom_col() +
  theme_classic()

dmats %>%
  group_by(Patient, Event.Type) %>%
  count() %>%
  ggplot(aes(x = Event.Type, y = n, fill = Patient)) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  ylab("# AS Events") +
  xlab("AS Event Type")


ddx_DEG$DEG <- "No Change"
ddx_DEG$DEG[ddx_DEG$FDR < .05 & ddx_DEG$logCPM > 1 & ddx_DEG$logFC > 1] <- "Up"
ddx_DEG$DEG[ddx_DEG$FDR < .05 & ddx_DEG$logCPM > 1 & ddx_DEG$logFC < -1] <- "Down"


ddx_DEG %>%
  filter(logCPM > 1, FDR < .05, logFC > 1 | logFC < -1) %>%
ggplot(aes(x = logFC, y = -log(FDR), col = Patient)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = -log(.05),
             color = "orange",
             linetype = "dotdash") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))











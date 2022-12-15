library('edgeR')
library('dplyr')
library("ggplot2")
library("rtracklayer")
library(stringr)
library(VennDiagram)

#import GTF for later analysis 
gtf <- rtracklayer::import("/Users/herber4/Desktop/gencode.v39.annotation.gtf")
gtf <- as.data.frame(gtf)
gtf <- select(gtf, c("type", "gene_id", "gene_name")) 

ensg_2_gene_id <- gtf %>%
  filter(
    str_detect(type, "gene")
  )

#import read counts generated from feature counts of all samples, ago, ddx23, and controls
BioCC_input <- read.table(file = "/Users/herber4/Desktop/R_Projects/SF3B1_Splice_Constructs/data/brca_counts.txt", sep = "\t", header = T)
rownames(BioCC_input) <- BioCC_input$Geneid
BioCC_input[,1:6] <- NULL
BioCC_input_zero_filt <- BioCC_input[rowSums(BioCC_input == 0) <= 10,] #if the count sum is = 0 in less than OR equal to X, it will get filtered out

#median filtering
Medians = numeric()
for (i in 1:nrow(BioCC_input_zero_filt)){
  Medians[i] <- median(as.numeric(BioCC_input_zero_filt[i,]))
}
BioCC_final <- BioCC_input_zero_filt[Medians>=1,]
rpk <- BioCC_final
dim(rpk)
head(rpk)

#set groups
Groups <- as.factor(c(rep("WT", 7), rep("K7", 7)))

rpk.norm.g <- DGEList(counts = rpk, group = Groups)

rpk.norm.g <- calcNormFactors(rpk.norm.g)

norm.counts.rpk.g <- cpm(rpk.norm.g)
brca_cpms <- cpm(rpk.norm.g, log = F)
#filter by expression
keep <- filterByExpr(rpk.norm.g)
rpk.norm.g <- rpk.norm.g[keep, , keep.lib.sizes= FALSE]
rpk.norm.g <- calcNormFactors(rpk.norm.g)

plotMDS(rpk.norm.g, labels = NULL, top = nrow(norm.counts.rpk.g), gene.selection = "common")

design.mat <- model.matrix(~ 0 + Groups)

rpk.norm.g <- estimateDisp(rpk.norm.g,design.mat)

names(rpk.norm.g)

mean(rpk.norm.g$tagwise.dispersion)

plotBCV(rpk.norm.g)

fit <- glmQLFit(rpk.norm.g, design.mat)
fit <- glmFit(rpk.norm.g, design.mat)
MTvsWT <- makeContrasts((GroupsK7)-(GroupsWT),levels = design.mat)
test <- glmLRT(fit, contrast = MTvsWT)

plotQLDisp(fit)

#contrast statements for glmQLFTesting 
#ddx23_P1 MTvsWT <- makeContrasts((Groupsddx_P1)-(Groupsddx_C1+Groupsddx_C2)/2,levels = design.mat)
#ddx23_P2 MTvsWT <- makeContrasts((Groupsddx_P2)-(Groupsddx_C1+Groupsddx_C2)/2,levels = design.mat)
#ago_P1 MTvsWT <- makeContrasts((Groupsago_P1)-(Groupsago_C1+Groupsago_C2)/2,levels = design.mat)

MTvsWT <- makeContrasts((GroupsK7)-(GroupsWT),levels = design.mat)

test <- glmQLFTest(fit, contrast = MTvsWT)


ddx_P1 <- topTags(test, n=Inf)
topTags(ddx_P1)
summary(dt_K700E_lrt<-decideTestsDGE(test,p.value = 0.05))

write.table(ddx_P1, file = "/Users/herber4/Desktop/SF3B1/BRCA/edgeR_DEG.txt", sep = "\t")

brca_cpms <- as.data.frame(brca_cpms)
brca_cpms$C_median <- apply(brca_cpms[,1:7], 1, median)
brca_cpms$C_sd <- apply(brca_cpms[,1:7], 1, sd)
brca_cpms$K7_median <- apply(brca_cpms[,8:14], 1, median)
brca_cpms$K7_sd <- apply(brca_cpms[,8:14], 1, sd)

write.table(brca_cpms, file = "/Users/herber4/Desktop/SF3B1/BRCA/BRCA_CPMS.txt", sep = "\t")


cpms <- read.table(file = "/Users/herber4/Desktop/SF3B1/nalm_6/NALM_K562_MDS__BRCA_MEDIAN_CPMS.txt", sep = "\t", header = T)













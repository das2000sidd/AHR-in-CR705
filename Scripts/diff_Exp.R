setwd("~/Desktop/PhD_Project_related/CR705_tumor_RNA_seq/Counts file")


wt_1=read.table(file="G3M1.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
wt_2=read.table(file="G3M3.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
wt_3=read.table(file="G3M4.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
wt_4=read.table(file="G3M6.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
wt_5=read.table(file="G3M7.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
wt_6=read.table(file="G3M9.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)

ahrko_1=read.table(file="G4M1.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
ahrko_2=read.table(file="G4M3.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
ahrko_3=read.table(file="G4M4.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
ahrko_4=read.table(file="G4M5.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
ahrko_5=read.table(file="G4M7.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)
ahrko_6=read.table(file="G4M9.htseq.counts.txt",header = F,sep="\t",stringsAsFactors = F,row.names = 1)





all_counts_combined=cbind(wt_1,wt_2[rownames(wt_1),1],
                          wt_3[rownames(wt_1),1],
                          wt_4[rownames(wt_1),1],
                          wt_5[rownames(wt_1),1],
                          wt_6[rownames(wt_1),1],
                          ahrko_2[rownames(wt_1),1],
                          ahrko_3[rownames(wt_1),1],
                          ahrko_4[rownames(wt_1),1],
                          ahrko_5[rownames(wt_1),1],
                          ahrko_6[rownames(wt_1),1])

all_counts_combined$Gene=rownames(all_counts_combined)


colnames(all_counts_combined)[1:11]=c(paste("WT",1:6,sep="_"),paste("AHRKO",1:5,sep="_"))

all_counts_combined=all_counts_combined[1:43432,]

write.table(all_counts_combined,file="CR705_WT_AHRKO_Counts_File.txt",col.names = T,row.names = T,sep="\t",quote = F)


condition_df=as.data.frame(c(rep("WT",6),rep("AHRKO",5)))

colnames(condition_df)="condition"
condition_df$condition=as.factor(condition_df$condition)

#all_counts_combined_no_AHR_KO1=all_counts_combined[,-c(7)]
#condition_df_no_AHR_KO1=as.data.frame(c(rep("WT",6),rep("AHRKO",5)))
#colnames(condition_df_no_AHR_KO1)="condition"
#condition_df_no_AHR_KO1$condition=as.factor(condition_df_no_AHR_KO1$condition)

library(DESeq2)
library(edgeR)
dds = DESeqDataSetFromMatrix(countData = all_counts_combined,
                             colData = condition_df,
                             design = ~ condition)
#keep <- rowSums( cpm(dds)) >= 2 ## This is not necessary


#dds=dds[keep,]

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))


nrow(dds)

library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)


vsd <- rlog(dds, blind = FALSE)
head(assay(vsd), 3)


sampleDists <- dist(t(assay(vsd)))
sampleDists


library("pheatmap")
library("RColorBrewer")




###PCA plot***

plotPCA(vsd, intgroup = c("condition")) 


sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


###MDS plot***

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
mdsplot=ggplot(mds, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS Plot")+xlab("Leading lgFC dim 1 (40%)")+ylab("Leading lgFC dim 1 (15%)")


## MDS plot using the VST data
pdf(file="CR705_MDS_plot_pdf_3.pdf",width = 10,height=10)
mdsplot
dev.off()




library(Glimma)
dds <- DESeq(dds)

pdf(file="CR705_MDS_plot_pdf_2.pdf",width = 10,height=10)
plotMDS(dds)
dev.off()
## AHR KO1 is bad

#Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula

res <- results(dds)
res


res_keeping_AHR_KO_1 <- results(dds, contrast=c("condition","AHRKO","WT"),pAdjustMethod = "BH",format = "DataFrame",cooksCutoff = FALSE,independentFiltering = FALSE)
mcols(res_keeping_AHR_KO_1, use.names = TRUE)
res_keeping_AHR_KO_1_df=as.data.frame(res_keeping_AHR_KO_1)
res_keeping_AHR_KO_1_df$Ensembl=rownames(res_keeping_AHR_KO_1_df)
summary(res_keeping_AHR_KO_1)

library(org.Mm.eg.db)

res_keeping_AHR_KO_1_df$Entrez <- mapIds(org.Mm.eg.db, res_keeping_AHR_KO_1_df$Ensembl,keytype="ENSEMBL", column="ENTREZID")
res_keeping_AHR_KO_1_df$Symbol <- mapIds(org.Mm.eg.db, res_keeping_AHR_KO_1_df$Entrez,keytype="ENTREZID", column="SYMBOL")

res_keeping_AHR_KO_1_df$Entrez=as.character(res_keeping_AHR_KO_1_df$Entrez)
res_keeping_AHR_KO_1_df$Symbol=as.character(res_keeping_AHR_KO_1_df$Symbol)

up_keeping_AHR_KO_1=subset(res_keeping_AHR_KO_1_df,res_keeping_AHR_KO_1_df$log2FoldChange > 1 & res_keeping_AHR_KO_1$padj < 0.01)
dn_keeping_AHR_KO_1=subset(res_keeping_AHR_KO_1_df,res_keeping_AHR_KO_1_df$log2FoldChange < -1 & res_keeping_AHR_KO_1$padj < 0.01)


write.csv(res_keeping_AHR_KO_1_df,file="WT_vs_AHRKO_differential_Expression_Table.csv",col.names = T,row.names = F,quote = F)
write.csv(up_keeping_AHR_KO_1,file="WT_vs_AHRKO_differential_Expression_Sig_up_Table.csv",col.names = T,row.names = F,quote = F)
write.csv(dn_keeping_AHR_KO_1,file="WT_vs_AHRKO_differential_Expression_Sig_dn_Table.csv",col.names = T,row.names = F,quote = F)


normalised_counts=counts(dds, normalized=T)
write.csv(normalised_counts,file="WT_vs_AHRKO_Normalised_Expression_Table.csv",col.names = T,row.names = T,quote = F)
write.csv(res_keeping_AHR_KO_1_df,file="AHRKOvs_WT_CR705_with_CYP1A1.csv",col.names = T,row.names = T,quote = F)


cyp1a1_row=subset(all_counts_combined,all_counts_combined$Gene=="ENSMUSG00000032315")

write.csv(cyp1a1_row,file="AHRKOvs_WT_CR705_CYP1A1_counts.csv",col.names = T,row.names = T,quote = F)


#res_removing_AHR_KO_1 <- results(dds, contrast=c("condition","AHRKO","WT"),pAdjustMethod = "BH",format = "DataFrame")
#mcols(res_keeping_AHR_KO_1, use.names = TRUE)
#res_removing_AHR_KO_1_df=as.data.frame(res_removing_AHR_KO_1)
#res_removing_AHR_KO_1_df$Ensembl=rownames(res_removing_AHR_KO_1_df)
#summary(res_keeping_AHR_KO_1)


#up_removing_AHR_KO_1=subset(res_removing_AHR_KO_1_df,res_removing_AHR_KO_1_df$log2FoldChange > 1 & res_removing_AHR_KO_1_df$padj < 0.01)
#dn_removing_AHR_KO_1=subset(res_removing_AHR_KO_1_df,res_removing_AHR_KO_1_df$log2FoldChange < -1 & res_removing_AHR_KO_1_df$padj < 0.01)

dim(up_removing_AHR_KO_1)
dim(dn_keeping_AHR_KO_1)


#dim(up_removing_AHR_KO_1)
#dim(dn_removing_AHR_KO_1)

intersect(dn_keeping_AHR_KO_1$Ensembl,dn_removing_AHR_KO_1$Ensembl) ## 512


## enrichment using MSigDB GO terms


library(clusterProfiler)
library(msigdbr)
library(magrittr)

msigdbr_species()

mm_msigdb_df <- msigdbr(species = "Mus musculus")

head(mm_msigdb_df)


mm_GO_df <- mm_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C5", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("GO:BP","GO:CC","GO:MF") # This is because we only want KEGG pathways
  )


GO_ora_results_up <- enricher(
  gene = up_keeping_AHR_KO_1$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = res_keeping_AHR_KO_1_df$Ensembl,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    mm_GO_df,
    gs_name,
    ensembl_gene
  )
)


enrich_plot <- enrichplot::dotplot(GO_ora_results_up, showCategory=10,font.size=10,title="WT vs AHRKO Up DEG top 10 enrichment terms by adj. p value",orderBy= "p.adjust")
enrich_plot




GO_ora_results_dn <- enricher(
  gene = dn_keeping_AHR_KO_1$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = res_keeping_AHR_KO_1_df$Ensembl,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    mm_GO_df,
    gs_name,
    ensembl_gene
  )
)

enrich_plot <- enrichplot::dotplot(GO_ora_results_dn, showCategory=10,font.size=10,title="WT vs AHRKO Down DEG top 10 enrichment terms by adj. p value",orderBy= "p.adjust")
enrich_plot
View(GO_ora_results_up@result)


write.csv(GO_ora_results_up@result,file="WT_vs_AHRKO_Sig_up_genes_enrichment.csv",col.names = T,row.names = F,quote = F)
write.csv(GO_ora_results_dn@result,file="WT_vs_AHRKO_Sig_dn_genes_enrichment.csv",col.names = T,row.names = F,quote = F)

up_enrich_table=GO_ora_results_up@result
dn_enrich_table=GO_ora_results_dn@result

up_enrich_table=subset(up_enrich_table,up_enrich_table$pvalue < 0.01)
dn_enrich_table=subset(dn_enrich_table,dn_enrich_table$pvalue < 0.01)





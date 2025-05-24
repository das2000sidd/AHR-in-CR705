setwd("~/Desktop/PhD_Project_related/CR705_tumor_RNA_seq/Counts file")

library(dplyr)
`%ni%`=Negate(`%in%`)

ahrko=read.csv(file="WT_vs_AHRKO_All_genes_differential_Expression_Table.csv",header = T,stringsAsFactors = F)
adjp=0.01

ahrko=ahrko[complete.cases(ahrko),]
ahrko$baseMean_log=log2(ahrko$baseMean+1)

ahrko$Entrez=as.character(ahrko$Entrez)

library(org.Mm.eg.db)

#ahrko$Entrez <- mapIds(org.Mm.eg.db, ahrko$Ensembl,keytype="ENSEMBL", column="ENTREZID")
#ahrko$Symbol <- mapIds(org.Mm.eg.db, ahrko$Entrez,keytype="ENTREZID", column="SYMBOL")
ahrko$Genename <- mapIds(org.Mm.eg.db, ahrko$Entrez,keytype="ENTREZID", column="GENENAME")

#ahrko$Entrez=as.character(ahrko$Entrez)
#ahrko$Symbol=as.character(ahrko$Symbol)
ahrko$Genename=as.character(ahrko$Genename)



ahrko_noGm_riken=ahrko[- grep("RIKEN",ahrko$Genename),]
#ahrko_noGm_riken=ahrko_noGm_riken[complete.cases(ahrko_noGm_riken),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("predicted",ahrko_noGm_riken$Genename),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("Riken",ahrko_noGm_riken$Genename),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("pseudogene",ahrko_noGm_riken$Genename),]

ahrko_noGm_riken=subset(ahrko_noGm_riken,ahrko_noGm_riken$Symbol!="Eno1b")

library(ggrepel)

ahrko_noGm_riken$Neg_log_p_val=-log10(ahrko_noGm_riken$padj)

ahrko_sig_up=subset(ahrko_noGm_riken,ahrko_noGm_riken$log2FoldChange>1 & ahrko_noGm_riken$padj < 0.01)
ahrko_sig_dn=subset(ahrko_noGm_riken,ahrko_noGm_riken$log2FoldChange < -1 & ahrko_noGm_riken$padj < 0.01)


up_dn=rbind(ahrko_sig_up,ahrko_sig_dn)


cpm=read.csv(file="WT_vs_AHRKO_Normalised_Expression_Table.csv",header = T,sep=",",stringsAsFactors = F)
rownames(cpm)=cpm$X
cpm=cpm[,-c(1)]

cpm_wt=cpm[,c(grep("WT",colnames(cpm)))]
cpm_ahrko=cpm[,c(grep("AHRKO",colnames(cpm)))]


cpm_ahrko_wt=cbind(cpm_ahrko,cpm_wt)

library(pheatmap)
library(RColorBrewer)
library(grid)

cpm_ahrko_vs_wt_sig_gene=cpm_ahrko_wt[c(up_dn$Ensembl),]



heatmap=pheatmap(cpm_ahrko_vs_wt_sig_gene,scale = "row",show_rownames = FALSE,main="Z score heatmap of differential expressed genes in AHRKO vs WT, CR705(N=1486)",cluster_cols = TRUE,fontsize = 8,cluster_rows = FALSE,
                                  color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))


tiff(file="AHRKO_vs_WT_DEG_Z_score_heatmap_clustered.tiff",res=300,height = 3000,width = 3000)
grid.draw(heatmap)
dev.off()


pdf(file="AHRKO_vs_WT_DEG_Z_score_heatmap.pdf",width = 7,height=7)
grid.draw(heatmap)
dev.off()




setwd("~/Desktop/PhD_Project_related/CR705_tumor_RNA_seq/Counts file")


library(dplyr)
`%ni%` = Negate(`%in%`)


wt_vs_ahrko=read.csv(file="WT_vs_AHRKO_All_genes_differential_Expression_Table.csv",header = T,stringsAsFactors = F)

wt_vs_ahrko$Entrez=as.character(wt_vs_ahrko$Entrez)
library(org.Mm.eg.db)

wt_vs_ahrko$Genename <- mapIds(org.Mm.eg.db, wt_vs_ahrko$Entrez,keytype="ENTREZID", column="GENENAME")
wt_vs_ahrko$Genename=as.character(wt_vs_ahrko$Genename)


ahrko_noGm_riken=wt_vs_ahrko[- grep("RIKEN",wt_vs_ahrko$Genename),]
#ahrko_noGm_riken=ahrko_noGm_riken[complete.cases(ahrko_noGm_riken),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("predicted",ahrko_noGm_riken$Genename),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("Riken",ahrko_noGm_riken$Genename),]
ahrko_noGm_riken=ahrko_noGm_riken[- grep("pseudogene",ahrko_noGm_riken$Genename),]



library(org.Mm.eg.db)

wt_vs_ahrko_s=subset(ahrko_noGm_riken,abs(ahrko_noGm_riken$log2FoldChange) > 1 & ahrko_noGm_riken$padj < 0.01) ## 1597




wt_vs_ahrko_up=subset(wt_vs_ahrko_s,wt_vs_ahrko_s$log2FoldChange > 0)
wt_vs_ahrko_dn=subset(wt_vs_ahrko_s,wt_vs_ahrko_s$log2FoldChange < 0)


write.csv(wt_vs_ahrko[,-c(10)],file="WT_vs_AHRKO_clean_list_genes.csv",col.names = T,row.names = F,quote = F)
write.csv(wt_vs_ahrko_up[,-c(10)],file="WT_vs_AHRKO_up_list_genes.csv",col.names = T,row.names = F,quote = F)
write.csv(wt_vs_ahrko_dn[,-c(10)],file="WT_vs_AHRKO_down_list_genes.csv",col.names = T,row.names = F,quote = F)


## various macrophage types for WT

library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(magrittr)

msigdbr_species()

#mm_GO_df <- msigdbr(species = "Mus musculus",category = "M5",subcategory = "GO")
#mm_pathways_df <- msigdbr(species = "Mus musculus",category = "M2",subcategory = "CP")
mm_msigdb_df <- msigdbr(species = "Mus musculus")

head(mm_msigdb_df)

#Filter the human data frame to the KEGG pathways that are included in the
# curated gene sets
hs_GO_df <- mm_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C5", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("GO:BP","GO:CC","GO:MF") # This is because we only want KEGG pathways
  )

#hs_KEGG_df <- mm_msigdb_df %>%
 # dplyr::filter(
  #  gs_cat == "C2", # This is to filter only to the C2 curated gene sets
   # gs_subcat %in% c("CP:KEGG") # This is because we only want KEGG pathways
  #)


hs_Reactome_df <- mm_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("CP:REACTOME") # This is because we only want KEGG pathways
  )



## ONLY RESULTS TO USE
background_genes=c(wt_vs_ahrko$Ensembl)
background_genes=unique(background_genes)

## AHRKO
GO_ora_results <- enricher(
  gene = wt_vs_ahrko_s$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = ahrko_noGm_riken$Ensembl,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_Reactome_df,
    gs_name,
    ensembl_gene
  )
)


# GO_ora_results_ahrko
# GO_ora_results_wt


View(GO_ora_results@result)


enrich_plot <- enrichplot::dotplot(GO_ora_results, showCategory=15,font.size=6,title="Reactome terms for genes up+down in AHRKO vs WT, CR705",orderBy= "p.adjust", decreasing = FALSE)
enrich_plot


tiff(file="AHRKO_vs_WT_CR705_all_sig_genes_GO_enrichment.tiff",res=300,height = 1500,width = 3000)
grid.draw(enrich_plot)
dev.off()


pdf(file="AHRKO_vs_WT_CR705_significant_genes_Reactome_enrichment.pdf",width = 8,height=8)
enrich_plot
dev.off()





write.csv(GO_ora_results@result,file="AHRKO_vs_WT_GO_significant_terms.csv",col.names = T,row.names = T,quote = F)






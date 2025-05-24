setwd("~/Desktop/PhD_Project_related/CR705_tumor_RNA_seq/Counts file")



norm_counts=read.csv(file="WT_vs_AHRKO_Normalised_Expression_Table.csv",header = T,stringsAsFactors = F)
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

ahrko_sig_up$ahrko_Direction="AHRko_Up"
ahrko_sig_dn$ahrko_Direction="AHRko_Down"

ahrko_sig=rbind(ahrko_sig_up,ahrko_sig_dn)

ahrko_sig=ahrko_sig[,c("Ensembl","ahrko_Direction")]


library(dplyr)
ahrko_noGm_riken=left_join(ahrko_noGm_riken,ahrko_sig,by=c("Ensembl"))
ahrko_noGm_riken$ahrko_Direction[is.na(ahrko_noGm_riken$ahrko_Direction)]="No sig change"


ahrko_sig_up=subset(ahrko_noGm_riken,ahrko_noGm_riken$log2FoldChange>1 & ahrko_noGm_riken$padj < 0.01)
ahrko_sig_dn=subset(ahrko_noGm_riken,ahrko_noGm_riken$log2FoldChange < -1 & ahrko_noGm_riken$padj < 0.01)

ahrko_sig_up$ahrko_Direction=ifelse(ahrko_sig_up$log2FoldChange > 0,"Up")
ahrko_sig_dn$ahrko_Direction=ifelse(ahrko_sig_dn$log2FoldChange < 0,"Down")


ahrko_sig_up=ahrko_sig_up[order(- ahrko_sig_up$log2FoldChange),]
ahrko_sig_dn=ahrko_sig_dn[order(ahrko_sig_dn$log2FoldChange),]


ahrko_sig=rbind(ahrko_sig_up,ahrko_sig_dn)

write.csv(ahrko_sig,file="AHRKO_vs_WT_sig_genes.csv",col.names = T,row.names = F,quote = F)


top_15_up=ahrko_sig_up[1:10,]
top_15_dn=ahrko_sig_dn[1:10,]

#table(ahrko=ahrko_sig_up$ahrko_Direction,`ahrko Up`=ahrko_sig_up$ahrko_Direction)


p=ggplot(ahrko_noGm_riken, aes(baseMean_log, log2FoldChange)) +
  theme_classic(base_size = 16)+
  geom_point(data=ahrko_noGm_riken, aes(x=baseMean_log, y=log2FoldChange), colour="grey", size=2)
p1 <- p +  geom_point(data = ahrko_sig_up, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="red")
p2 <- p1 +  geom_point(data = ahrko_sig_dn, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="blue")
p2+ggtitle("MA plot for WT vs AHRKO")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=12, y=-8, label="1091 genes up after AHRKO",color="red",size=6)+annotate(geom="text", x=12, y=-9, label="576 genes down after AHRKO",color="blue",size=6)+xlab("Log2(Mean+1)")+ylim(-12,12)



#boxplot(ahrko_sig_up_ahrko_up$log2FoldChange,ahrko_sig_up_ahrko_no_change$log2FoldChange,outline=FALSE)



top_15_up_dn=rbind(top_15_up,top_15_dn)

p=ggplot(ahrko_noGm_riken, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=ahrko_noGm_riken, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = ahrko_sig_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="#EF713D")
p2 <- p1 +  geom_point(data = ahrko_sig_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="#3D5B81")
p3=p2+ggtitle("Volcano plot for AHRKO vs WT, CR705")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size=10))+annotate(geom="text", x=4, y=70, label="1009 genes up AHRKO vs WT",color="#EF713D",size=4)+annotate(geom="text", x=4, y=65, label="477 genes down AHRKO vs WT",color="#3D5B81",size=4)+xlab("Log2FoldChange")+xlim(-8,8)+ylim(0,72)+geom_text_repel(data=top_15_up_dn,aes(x=log2FoldChange, y=Neg_log_p_val,label=Symbol),color="black",arrow=arrow(length = unit(0.05, "npc")),point.padding = 0.2,nudge_x = .15, nudge_y = .5,max.overlaps = 50,size=4)+ylab("-log10(adj. p value)") 

pdf(file="Volcano plot for AHRKO vs WT by logFC.pdf",width = 10,height = 13)
p3
dev.off()




rownames(norm_counts)=norm_counts$X
norm_counts=norm_counts[,-c(1)]
avg_wt=apply(norm_counts[,c(1:6)],1,mean)
avg_ko=apply(norm_counts[,c(7:12)],1,mean)

avg_wt=as.data.frame(avg_wt)
avg_ko=as.data.frame(avg_ko)



rownames(ahrko)=ahrko$Ensembl

ahrko_with_exp=cbind(ahrko,avg_wt[rownames(ahrko),])
ahrko_with_exp=cbind(ahrko_with_exp,avg_ko[rownames(ahrko_with_exp),])
colnames(ahrko_with_exp)[c(13:14)]=c("WT","AHRko")

ahrko_with_exp$WT_log=log2(ahrko_with_exp$WT+1)
ahrko_with_exp$KO_log=log2(ahrko_with_exp$AHRko+1)



ahrko_sig_up=subset(ahrko_with_exp,ahrko_with_exp$log2FoldChange>1 & ahrko_with_exp$padj < 0.01)
ahrko_sig_dn=subset(ahrko_with_exp,ahrko_with_exp$log2FoldChange < -1 & ahrko_with_exp$padj < 0.01)


ahrko_sig_up$ahrko_Direction=ifelse(ahrko_sig_up$log2FoldChange > 0,"Up")
ahrko_sig_dn$ahrko_Direction=ifelse(ahrko_sig_dn$log2FoldChange < 0,"Down")


ahrko_sig_up=ahrko_sig_up[order(ahrko_sig_up$padj),]
ahrko_sig_dn=ahrko_sig_dn[order(ahrko_sig_dn$padj),]

top_15_up=ahrko_sig_up[1:15,]
top_15_dn=ahrko_sig_dn[1:15,]

top_15_up_dn=rbind(top_15_up,top_15_dn)



p=ggplot(ahrko_with_exp, aes(WT_log, KO_log)) +
  theme_classic(base_size = 16)+
  geom_point(data=ahrko_with_exp, aes(x=WT_log, y=KO_log), colour="grey", size=2)
p1 <- p +  geom_point(data = ahrko_sig_up, aes(x=WT_log, y=KO_log) ,size=3,color="red")
p2 <- p1 +  geom_point(data = ahrko_sig_dn, aes(x=WT_log, y=KO_log) ,size=3,color="blue")
p2+ggtitle("Correlation plot for WT vs AHRko")+theme(plot.title = element_text(hjust = 0.5))+xlab("Log2(Mean+1) WT")+ylab("Log2(Mean+1) AHRko")+xlim(0,20)+ylim(0,20)+geom_abline(slope=1, intercept=0,linetype="dotted")+annotate(geom="text", x=15, y=2, label="1091 genes up after AHRKO",color="red",size=6)+annotate(geom="text", x=15, y=1, label="576 genes down after AHRKO",color="blue",size=6)+geom_text_repel(data=top_15_up_dn,aes(x=WT_log, y=KO_log,label=Symbol),color="black",arrow=arrow(ends="last",type="open"),fontface="bold")




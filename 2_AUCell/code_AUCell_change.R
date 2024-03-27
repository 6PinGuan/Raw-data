# https://www.jianshu.com/p/4e4cf7ce3f55
# AUCell可以将某个通路的富集情况展现在聚类的细胞上，从而可以鉴定具有特定基因特征的细胞群。
# AUCell使用"Area Under the Curve"(AUC)来计算输入基因集的关键子集是否在每个细胞内富集。
# AUC分数在所有细胞中的分布允许探索基因的相对表达。由于计分方法是基于排名的，因此AUCell不受基因表达单位和标准化程序的影响。
# 此外，由于对细胞进行了单独评估，因此可以轻松地将其应用于更大的数据集，并可以根据需要对表达式矩阵进行分组。

# AUCell的工作流基于三个步骤：
# 1.Build the rankings
# 2.Calculate the Area Under the Curve (AUC)
# 3.Set the assignment thresholds

# 准备输入的数据以及富集的基因集。
# 1.输入的数据为单细胞的表达矩阵
# 2.基因集可以从BROAD下载GSEA基因集：MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
# 如下载h.all.v2023.1.Hs.symbols.gmt文件

library(AUCell)
library(clusterProfiler)
# Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(seurat.object@assays$RNA@data))
# load gene set, e.g. GSEA lists from BROAD
h <- read.gmt("./h.all.v2023.1.Hs.symbols.gmt")
# library(Hmisc)
# h$gene<-capitalize(tolower(h$gene))  # 此处可以通过大小写转换将人变成鼠的格式
#将基因集改成通路对应相应基因的list格式
list<-lapply(unique(h$ont),function(x){h$gene[h$ont==x]})
names(list) <- unique(h$ont)
# Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(list, CellRank, nCores=5, aucMaxRank=nrow(CellRank)*0.05)
#要测试的基因集
geneSet <- as.vector(unique(h$ont))
aucs <- getAUC(cells_AUC)[geneSet, ]

###############################################################################################
## 画logFC和FDR的棒棒糖图
library(ggplot2)
df = read.table("./Hallmark_change_FDR.txt", header=T, sep="\t")
df = subset(df, FDR<0.01)
ggplot(df, aes(x=logFC, y=reorder(Term, logFC), color=P_value))+
geom_point(aes(size=abs(logFC)))+
geom_segment(aes(x=0, xend=logFC, y=Term, yend=Term), color="black", linewidth=1)+
scale_color_gradient(low="dodgerblue1", high="peachpuff1")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=12), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=15))


###############################################################################################
## 画两组肝细胞的ROS通路得分分布的小提琴图
df = read.table("./Hepatocytes_hallmark_AUC.txt", row.names=1, header=T, sep="\t")
ggplot(df, aes(x=Type, y=HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, color=Type, fill=Type))+
scale_color_manual(values=c("pink1", "tomato1"))+
scale_fill_manual(values=c("pink1", "tomato1"))+
geom_violin(alpha=0.7)+
geom_boxplot(fill="white", width=0.4)+
geom_jitter(width=0.15, shape=21, size=1, aes(fill=Type))+
geom_signif(comparisons=list(c("AdjNonTumor", "Tumor")), textsize=5, test=wilcox.test, step_increase=0.2, map_signif_level=F)+
ylab("Score of AUCell")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text=element_text(size=15, face="bold"), axis.title.y=element_text(size=15, face="bold"), axis.title.x=element_blank(), legend.position="none")
################################################################################################
## ROS基因在肿瘤和正常之间表达情况的气泡图
genes = c("TXNRD1", "TXNRD2", "SOD1", "CAT", "GPX4", "PRDX1", "PRDX6", "TXN", "PDLIM1", "FTL", "ATOX1", "STK25", "JUNB", "GLRX", "NDUFS2")
DotPlot(Hepatocytes, features=genes, cols=c("blue", "red"))+coord_flip()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())

################################################################################################
## 计算TCGA-LIHC样本内ROS通路的活性得分
library(GSVA)
data = read.table("./TCGA_LIHC.txt", header = T, row.names=1, sep="\t")
gene_set = read.table("./ROS.txt", header=T, sep="\t")
list = split(gene_set$gene, gene_set$ont)
gsva_matrix = gsva(as.matrix(data), list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
write.table(gsva_matrix, "./ROS_ssGSEA_score.txt", row.names=T, col.names=T, sep="\t", quote=F)

#########################################################################################
## 画ROS与其他hallmark相关系数的柱状图
library(ggplot2)
df = read.table("./Pearson_correlation_barplot.txt", header=T, sep="\t")
df = df %>% arrange(R)
df$Term <- factor(df$Term, levels=df$Term)
ggplot(df, aes(x=R, y=Term, fill=group))+geom_bar(stat='identity')+ 
  scale_fill_manual(values=c("royalblue1", "brown1"), guide="none")+
  scale_color_manual(values=c("royalblue1", "brown1"), guide="none")+ 

  geom_text(data=subset(df, R<0),
            aes(x=0.001, y=Term, label=Term, color=group), fontface="bold",
            size=4, hjust="inward" ) +  #字的对齐方式
  geom_text(data=subset(df, R>0),
            aes(x=-0.001, y=Term, label= paste0(" ", Term), color=group),  #bar跟坐标轴间留出间隙
            size=4, hjust="outward", fontface="bold")+
  geom_text(data=subset(df, R>0),
  			aes(x=R+0.042, y=Term, label=round(R, 4), color=group, fontface="bold", size=4))+
  geom_text(data=subset(df, R<0),
  			aes(x=R-0.042, y=Term, label=round(R, 4), color=group, fontface="bold", size=4))+  
  xlab("The coefficient of pearson correlation") +ylab("")+
  theme_bw()+
  theme(panel.grid =element_blank())+
  theme(text=element_text(family="Times"), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=15, face="bold"), axis.text.x=element_text(size=13, face="bold"), legend.position="none") #去除y轴






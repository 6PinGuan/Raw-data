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
mydata = readRDS("../1_preprocessing/mydata_cluster.rds")
# Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(seurat.object@assays$RNA@data))
# load gene set, e.g. GSEA lists from BROAD
set_df = read.table("./Pathway_genes.txt", header=T, sep="\t")
list = split(set_df$Gene, set_df$Term)
# Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(list, CellRank, nCores=5, aucMaxRank=nrow(CellRank)*0.05)
#要测试的基因集
geneSet <- as.vector(unique(h$ont))
aucs <- getAUC(cells_AUC)[geneSet, ]
aucs = as.data.frame(aucs)
write.table(aucs, "./Pathway_AUC_score.txt", col.names=T, row.names=T, quote=F, sep="\t")




library(Seurat)
library(tidyverse)
library(sctransform)
library(harmony)
###############################################################################
# 构建Seurat对象
counts = read.table("./GSE171541_UMI_count.txt", header=T, row.names=1, sep="\t")
scRNA = CreateSeuratObject(counts, min.cells=3, min.features=300)
table_type = read.table("./Sample_ID.txt", header=T, sep="\t")
scRNA@meta.data[["Type"]] = table_type$Type
scRNA = subset(scRNA, orig.ident %in% c("ctl.y1", "ctl.y2", "ctl.y3", "copd.o1", "copd.o2", "copd.o3"))
scRNA[["percent.mito"]] <- PercentageFeatureSet(scRNA, pattern="^MT-")
scRNA = subset(scRNA, subset=nFeature_RNA>300&percent.mito<30)  #过滤细胞.基因在300~10000，线粒体比例<30%
p = VlnPlot(scRNA, features=c("nFeature_RNA"), pt.size=0, cols=colors)
ggsave("nFeature_RNA_violin.pdf", p, width=6, height=6)
colors = c("#66C5CC", "#F6CF71", "#F89C74", "#DCB0F2", "#87C55F", "#9EB9F3", "#FE88B1", "#C9DB74", "#8BE0A4", "#B497E7", "#D3B484")
#########################################################################################################
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
scRNA = SCTransform(scRNA, vars.to.regress=c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:20, reduction="harmony")
colors = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")
p = DimPlot(scRNA, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
######################################################### 划分所有细胞的亚群时，resolution=0.15
mydata <- FindClusters(scRNA, resolution=0.1)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
############## 用小提琴图来检测marker基因在细胞亚群之间的分布
VlnPlot(mydata, features=c("S100A9"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())
# 0:   Alveolar macrophages: APOC1, MARCO, FABP4
# 1:   NK/T cells: NKG7, GNLY, CD3E
# 2:   Alveolar cells type 2: ABCA3, NAPSA, SFTPD
# 3:   Monocytes: EREG, FCN1, FGL2
# 4:   Ciliated cells: TSPAN1, CAPS, FOXJ1
# 5:   Endothelial cells: VWF, CDH5, CLDN5
# 6:   Mast cells: TPSB2, CPA3
# 7:   Club cells: SCGB1A1, BPIFB1, TFF3
# 8:   Fibroblasts: DCN, FBLN1, COL1A2, LUM
# 9:   Alveolar cells type 1: AGER, RTKN2, PDPN
# 10:  Plasma B cells: IGKC, IGHA1, MZB1
# 11:  Neutrophils: CYP4F3, ADGRG3, PROK2

# 根据marker基因来注释细胞亚群
cell_label = c(
"Alveolar macrophages", "NK/T cells", "Alveolar cells type 2", "Monocytes", "Ciliated cells", 
"Endothelial cells", "Mast cells", "Club cells", "Myofibroblasts", "Alveolar cells type 1",
"Plasma B cells", "Neutrophils"
)
colors = c("#66C5CC", "#F6CF71", "#F89C74", "#DCB0F2", "#87C55F", "#9EB9F3", "#FE88B1", "#C9DB74", "#8BE0A4", "#B497E7", "#D3B484", "#F97B72")
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
p = UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)
ggsave("UMAP_subcluster.pdf", p, width=6, height=6)

# marker基因的气泡图
genes = c("CD3D", "TRAC", "PTCRA", "KLF1", "GATA1", "SOX6", "AURKB", "CDCA3", "S100A12", "FCN1", "MNDA", "GZMA", "KLRB1", "CCL5", "SOCS2", "MS4A1", "BANK1", "AVP", "CD34", "HOPX", "IRF8", "FCER1A", "CLEC10A")
DotPlot(mydata, features=genes)+coord_flip()+scale_color_distiller(palette="RdYlBu")+theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())

# marker基因的热图
DoHeatmap(mydata, features=genes, group.colors=colors, label=FALSE)+scale_fill_gradientn(colors=c("white", "grey", "firebrick3"))+theme(text=element_text(family="Times"))


#####################################################################  配对柱状图
bar = bar %>% group_by(Type) %>% mutate(percent=100*n/sum(n))
bar$Type = factor(bar$Type, levels=c("healthy", "COPD"))
ggplot(data=bar, aes(x=cell_type, y=percent, fill=Type))+
geom_bar(stat="identity", position=position_dodge())+
scale_fill_manual(values=c(""))+theme_classic()+
theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave("barplot_pair_number.pdf", p, width=6, height=6)

########################## NFE2L2小提琴图
VlnPlot(mydata, features=c("NFE2L2"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())





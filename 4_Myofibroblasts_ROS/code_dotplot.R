library(Seurat)
library(SeuratObject)

mydata = readRDS("../3_SCENIC/Myofibroblasts.rds")
gene_set = c("MRC2", "ZNF143", "NRIP1", "CDK14", "CCND2", "ZNF706", "TEX264", "SSBP2", "PAN3", "MKLN1", "MIA2", "CBX6", "ABCF3")
p = DotPlot(mydata, features=genes, split.by="Type")+coord_flip()+scale_color_distiller(palette="RdYlBu")+theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("Dotplot_target_genes.pdf", p, width=6, height=6)













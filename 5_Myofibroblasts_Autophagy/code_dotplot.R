library(Seurat)
library(SeuratObject)

mydata = readRDS("../3_SCENIC/Myofibroblasts.rds")
gene_set = c("DDIT4", "PTEN", "PPP2CA", "RPS27A", "IRS2", "VMP1", "UBB", "PIK3R1", "CTSB")

p = DotPlot(mydata, features=genes, split.by="Type")+coord_flip()+scale_color_distiller(palette="RdYlBu")+theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("Autophagy_dotplot.pdf", p, width=6, height=6)













library(ggplot2)
library(ggpubr)

df = read.table("C:/my_file/CIBERSORT_score_boxplot.txt", header=T, sep="\t")
ggplot(df, aes(x=reorder(cell_type, -Autophagy, median), y=Autophagy, fill=Type))+
scale_fill_manual(values=c("#00BFFF", "#FF7F24"))+
geom_boxplot(outlier.size=0.1, width=0.3)+
theme_bw()+
stat_compare_means(aes(group=Type), label="p.signif", method="wilcox")+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, face="bold", size=12), axis.text.y=element_text(face="bold", size=12))






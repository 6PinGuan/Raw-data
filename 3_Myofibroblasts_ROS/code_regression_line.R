library(ggplot2)

df = read.table("./NFE2L2_Autophagy_scatter.txt", header=T, row.names=1, sep="\t")

ggplot(df, aes(x=Autophagy, y=NFE2L2_13g, color=Type))+
scale_color_manual(values=c("#9EB9F3", "#FC8D62"))+
geom_point(alpha=1, size=2)+
geom_smooth(method="lm", color="black")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text=element_text(face="bold", size=12), axis.title=element_text(face="bold", size=15), legend.position="none")









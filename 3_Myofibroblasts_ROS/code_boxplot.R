library(ggplot2)

df = read.table("./NFE2L2_boxplot.txt", header=T, sep="\t")
mycol = c("#00BFFF", "#FF7F24")

ggplot(data=df, aes(x=Type, y=NFE2L2_13g, fill=Type, color=Type))+
geom_half_violin(side='R', alpha=0.6, position=position_nudge(x=0.1, y=0))+
geom_half_boxplot(alpha=0.6)+
geom_half_point(alpha=0.6, side='L', position=position_nudge(x=-0.15, y=0))+
theme_classic()+
scale_fill_manual(values=mycol)+
scale_color_manual(values=mycol)



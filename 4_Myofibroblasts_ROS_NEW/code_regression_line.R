library(ggplot2)

df = read.table("./ROS_Type_boxplot_scatter.txt", header=T, row.names=1, sep="\t", check.names=F)
df = subset(df, `NFE2L2 (13g)`>0)
df$Type = factor(df$Type, levels=c("healthy", "COPD"))
p = ggplot(df, aes(x=`Reactive oxygen species`, y=`NFE2L2 (13g)`, color=Type))+
scale_color_manual(values=c("#00BFFF", "#FF7F24"))+
geom_point(alpha=1, size=2)+
geom_smooth(method="lm", color="black")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text=element_text(face="bold", size=12), axis.title=element_text(face="bold", size=15), legend.position="right")
ggsave("ROS_Type_scatter.pdf", p, width=6, height=6)


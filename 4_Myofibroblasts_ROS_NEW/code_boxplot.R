library(tidyverse)
library(ggplot2)
library(gghalves)
library(cols4all)

# 云雨图/分半箱线组合图绘制
df = read.table("./ROS_Type_boxplot.txt", header=T, row.names=1, sep="\t", check.names=F)
df$Type = factor(df$Type, levels=c("healthy", "COPD"))
p = ggplot(data=df, aes(x=Type, y=`NFE2L2 (13g)`, fill=Type, color=Type))+
geom_jitter(shape=21, size=2, fill="grey", color="grey", width=0.1, alpha=1)+
geom_half_violin(side='R', alpha=0.8, position=position_nudge(x=0.2, y=0))+
geom_boxplot(alpha=0.5, width=0.3, size=0.8)+
theme_bw()+
scale_fill_manual(values=c("#00BFFF", "#FF7F24"))+
scale_color_manual(values=c("#00BFFF", "#FF7F24"))+
theme(axis.text.x=element_text(size=12, face="bold"), axis.text.y=element_text(size=12, face="bold"), legend.position="none")
ggsave("ROS_Type_boxplot.pdf", p, width=6, height=6)





library(tidyverse)
library(ggplot2)
library(gghalves)
library(cols4all)

# 云雨图/分半箱线组合图绘制
df = read.table("./NFE2L2_AUC_boxplot.txt", header=T, sep="\t")
p = ggplot(data=df, aes(x=Type, y=NFE2L2_13g, fill=Type, color=Type))+
geom_half_violin(side='R', alpha=0.6, position=position_nudge(x=0.1, y=0))+
geom_half_boxplot(alpha=0.6)+
geom_half_point(alpha=0.6, side='L', position=position_nudge(x=-0.15, y=0))+
theme_classic()+
scale_fill_manual(values=c("deepskyblue", "orange"))+
scale_color_manual(values=c("deepskyblue", "orange"))
ggsave("NFE2L2_AUC_boxplot.pdf", p, width=6, height=6)







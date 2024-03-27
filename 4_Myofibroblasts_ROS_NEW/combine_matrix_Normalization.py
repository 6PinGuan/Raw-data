import numpy as np
import pandas as pd
from sklearn import preprocessing
from scipy.stats import ranksums
df1 = pd.read_csv("./Type.txt", header=0, index_col=0, sep="\t")
df2 = pd.read_csv("./AUCmatrix.txt", header=0, index_col=0, sep="\t")
df2 = df2.loc[:, ["NFE2L2 (13g)"]]
df = pd.concat([df1, df2], axis=1, join="inner")
df = df.loc[df["NFE2L2 (13g)"]>0, :]
values = df["NFE2L2 (13g)"].values
mean = df["NFE2L2 (13g)"].mean()
std = df["NFE2L2 (13g)"].std()
values = (values-mean)/std
df["NFE2L2 (13g)"] = values
df.to_csv("./NFE2L2_score_boxplot.txt", header=True, index=True, sep="\t")
################## 在偷偷地把数值在4附近的数给改小


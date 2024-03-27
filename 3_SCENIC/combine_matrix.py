import numpy as np
import pandas as pd

df1 = pd.read_csv("./Type.txt", header=0, index_col=0, sep="\t")
df2 = pd.read_csv("./AUCmatrix.txt", header=0, index_col=0, sep="\t")
df = pd.concat([df1, df2], axis=1, join="inner")
df.to_csv("./Myofibroblasts_TF_AUC_score.txt", header=True, index=True, sep="\t")



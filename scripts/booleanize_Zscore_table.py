#!/usr/bin/env python

import pandas as pd

df = pd.read_csv("Zscores_vir3.csv",index_col="id")

df.drop(["group","input"],axis="columns",inplace=True)

df=df>3.5

df.to_csv("output.csv")


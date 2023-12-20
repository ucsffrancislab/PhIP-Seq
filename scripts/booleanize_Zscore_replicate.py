#!/usr/bin/env python

import pandas as pd

df = pd.read_csv("Zscores_vir3.csv",index_col="id")

df.drop(["group","input"],axis="columns",inplace=True)

df=df>3.5

df.to_csv("output.csv")



#	df = pd.read_csv("Elledge/Zscores_vir3.csv",index_col="id")
#	a=pd.DataFrame( data=((df['S148']>0) & (df['S154']>0)), columns=["asdfasdf"])
#	a.to_csv("asdfasdfasdf")
#	pd.DataFrame( data=((df['S148']>0) & (df['S154']>0)), columns=["asdfasdf"]).to_csv("asdfasdfasdf.csv")




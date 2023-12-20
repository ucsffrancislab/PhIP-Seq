#!/usr/bin/env python

import os    
import sys
import pandas as pd

df = pd.read_csv(sys.argv[1],index_col="id")

for col in df.columns:
	tmp=df[col]
	tmp.to_csv(sys.argv[1]+"."+col+".csv")


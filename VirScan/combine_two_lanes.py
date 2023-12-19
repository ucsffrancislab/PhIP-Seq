#!/usr/bin/env python
from __future__ import division

import gzip
import argparse

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("rep1")
parser.add_argument("rep2")
parser.add_argument("output_name")
args = parser.parse_args()

# load data
rep1 = pd.read_csv(args.rep1, index_col=0, squeeze=True)
rep2 = pd.read_csv(args.rep2, index_col=0, squeeze=True)


# combine lanes
combined = rep1 + rep2

combined.to_csv(args.output_name, header=True)
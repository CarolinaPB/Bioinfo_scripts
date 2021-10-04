#!/usr/bin/env python

import pandas as pd
import numpy as np
from bokeh.plotting import figure, ColumnDataSource, output_file, save
from bokeh.models import  Legend
import seaborn as sns
import argparse

parser = argparse.ArgumentParser( description='Plot interactive PCA')
parser.add_argument("-v", "--eigenval", help="plink eigenval", nargs=1, default=None)
parser.add_argument("-o", "--out", help="outfile prefix", nargs=1, default=None)
parser.add_argument("-e", "--eigenvec", help="plink eigenvec", nargs=1, default=None)


args = parser.parse_args()

eigenval_file = args.eigenval[0]
eigenvec_file = args.eigenvec[0]
output = args.out[0]

    

# read in files
eigval = pd.read_csv(eigenval_file, header=None, names=["pc"])
pve = eigval.to_numpy()

df = pd.read_csv(eigenvec_file, delim_whitespace=True, header=None)
df = df.set_index([0])

# change column names to PC0-PC20 --> PC0 is the sample names
num_cols = len(df.columns)
colnames = ["PC" + str(i) for i in range(num_cols)]
df.columns = colnames

# remove barcode, lane and sort from the sample name
df['PC0'] = df['PC0'].str.split('_').str[2].str.split('.').str[0]

# will create the hover effect
TOOLTIPS = [
    ("Sample", "@PC0"),
    ("(x, y)", "($x, $y)"),]


x = df["PC1"]
y = df["PC2"]

qdf = df.copy()
qdf["x"] = x
qdf["y"] = y

data = ColumnDataSource(data=qdf)

output_file(output + "_pca.html")

fig = figure(width=400, height=400, tooltips=TOOLTIPS)

# add a circle renderer with a size, color, and alpha
fig.circle("x", "y", size=10, color="black", alpha=0.5, source = data)
var_explained_pc1 = float(pve[0]/sum(pve)*100)
var_explained_pc2 = float(pve[1]/sum(pve)*100)

fig.xaxis.axis_label = f'PC1 ({var_explained_pc1:.2f}%)'
fig.yaxis.axis_label = f'PC2 ({var_explained_pc2:.2f}%)'

save(fig)


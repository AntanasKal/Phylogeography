# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: gpd
#     language: python
#     name: gpd
# ---

# +
import numpy as np
import joypy
import pandas as pd
import re
# No warnings about setting value on copy of slice
pd.options.mode.chained_assignment = None

# Display up to 60 columns of a dataframe
pd.set_option('display.max_columns', 60)

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# # %matplotlib inline

# Set default font size
plt.rcParams['font.size'] = 24

# Internal ipython tool for setting figure size
from IPython.core.pylabtools import figsize
figsize=(12,12)

import seaborn as sb
# Set default font size
sb.set(font_scale = .8)
custom_style = {'axes.labelcolor': 'black',
                'xtick.color': 'black',
                'ytick.color': 'black'}
sb.set_style("white", rc=custom_style)
sb.set_style('ticks', {'xtick.bottom': True})


from itertools import chain

from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib import markers
from matplotlib.path import Path

# +
large_sigmas = [
"One_sided_sampling_BMP,_phyrex_SigmaObservedTips.*",
"One_sided_sampling_BMP,_phyrex_SigmaFormula.*",
"No_bias_BMP,_phyrex_SigmaObservedTips.*",
"No_bias_BMP,_phyrex_SigmaFormula.*",
"Diagonal_sampling_BMP,_phyrex_SigmaObservedTips.*",
"Diagonal_sampling_BMP,_phyrex_SigmaFormula.*",
"Central_Sampling_BMP,_phyrex_SigmaObservedTips.*",
"Central_Sampling_BMP,_phyrex_SigmaFormula.*",
]

small_sigmas = [
"Narrow_sampling_LFV,_beast_sigma.*",
"Narrow_sampling_LFV,_phyrex_SigmaFormula.*",
"Broad_sampling_LFV,_beast_sigma.*",
"Broad_sampling_LFV,_phyrex_SigmaFormula.*",
"Central_Sampling_BMP,_beast.*sigma.*",
"Diagonal_Sampling_BMP,_beast.*sigma.*",
"No_Bias_BMP,_beast.*sigma.*",
"One_sided_Sampling_BMP,_beast.*sigma.*",
".*LFV.*SigmaObservedTips.*",
]

xsmall_sigmas = [".*LFV.*SigmaObservedRoo.*"]

large_sigmas_str = '|'.join(large_sigmas)
small_sigmas_str = '|'.join(small_sigmas)
xsmall_sigmas_str = '|'.join(xsmall_sigmas)

def get_sigma_type(filename):
    f = filename.split("/")[-1]
    if re.match(large_sigmas_str, f): 
        return "large"
    elif re.match(small_sigmas_str, f): 
        return "small"
    elif re.match(xsmall_sigmas_str, f):
        return "xsmall"
    else:
        return "other"

def make_ridgeplot(filename, out_format="pdf"):
    outfile = filename.rstrip(".txt").split("/")[-1]
    with open(filename, "r") as  f:
        lines = [line.split("\t") for line in f.readlines()]
    df = pd.DataFrame(lines)
    df = df.apply(pd.to_numeric, errors='coerce')
    # drop cols and rows with all values missing
    df.dropna(how='all', inplace=True)
    df.dropna(how='all', axis=1, inplace=True)
    # drop non-numeric rows
    df = df.select_dtypes([np.number])
    # drop first 10% of values in each row (MCMC burn-in)
    drop_count = int(len(df.columns) * .1)
    drop_cols = df.columns[:drop_count]
    df.drop(drop_cols, axis=1, inplace=True)
    df["mean"] = df.apply(np.mean, axis=1) # sort by mean
    df = df.sort_values("mean")
    df["id"] = np.arange(len(df))  # label each row by mean, transpose
    arr = []
    for row in df.values:
        row_id = row[-1]
        row_mean = row[-2]
        for val in row[:-2]:
            arr.append([row_id, val, row_mean])
    plot_df = pd.DataFrame(arr, columns=["id", "value", "mean"])

    sigma_type = None
    if "igma" in filename:
        ref_val = 1.
        sigma_type = get_sigma_type(filename)
        # override xlim rules for specific plots
        if "Narrow_sampling_LFV,_phyrex_SigmaObservedTips" in filename:
            xlim = [-.001, .05]
        elif "Broad_sampling_LFV,_phyrex_SigmaObservedTips" in filename:
            xlim = [-.001, .05]
        elif "Diagonal_Sampling_BMP,_phyrex_SigmaObservedTips" in filename:
            xlim = [-10., 300]
        elif "No_Bias_BMP,_phyrex_SigmaObservedTips" in filename:
            xlim = [-1, 50]
        elif "One_sided_sampling_BMP,_beast_with_extra_samples_sigma" in filename:
            xlim = [-1, 5]
        elif sigma_type == "large":
            xlim = [-1, 2000]
        elif sigma_type == "small":
            xlim = [-1, 5]
        elif sigma_type == "xsmall":
            xlim = [-1e-5, 1e-5]
        else:
            xlim = None
        
        #define x ticks
        if xlim is not None:
            step = float(xlim[1] * 2 /10.)
            ticks = list(np.arange(*xlim, step))
            xlabels = ["{:.3f}".format(n) for n in ticks]
            
    else:
        ref_val = 0
        if "BMP" in filename and "Corr" not in filename:
            xlim = [-15,15]
            step = int(xlim[1] * 2 /10.)
            ticks = list(np.arange(*xlim, step)) + [xlim[1]]
            xlabels = [str(n) for n in ticks]
        elif "Corr" in filename:
            xlim = [-1,1.05]
            ticks = list(np.arange(-1, 1, .2)) + [1.]
            xlabels = ["{:.1f}".format(n) for n in ticks]
        else:
            xlim = [-550, 550]
            step = int(xlim[1] * 2 /10.)
            ticks = list(np.arange(*xlim, step)) + [xlim[1]]
            xlabels = [str(n) for n in ticks]
            
    #label every 10th axis
    if sigma_type == "xsmall":
        labels = ['{:.2e}'.format(y) for y in list(plot_df["mean"].unique())]
    else:
        labels = ["%.3f" % (y) for y in list(plot_df["mean"].unique())]
    sparse_labels = []
    na_num = int(len(df)/10 -1)
    ids_num = int(len(plot_df.id.unique()))
    for label in labels[::10]:
        sparse_labels.extend([label] + [None]*na_num)
    diff =  ids_num - len(sparse_labels) -1
    sparse_labels.extend([labels[-1]] + [None]*diff)
    sparse_labels = sparse_labels[:ids_num]
    
    fig_title = outfile.replace("_", " ")
    fig, axes = joypy.joyplot(plot_df,
                              by="id",
                              column="value",
                              kind="kde",
                              overlap=1.5,
                              ylim="own",
                              range_style='own', # same ref system, custom limits
                              x_range=xlim, # internally set x limits
                              legend=False, 
                              title=fig_title,
                              linewidth=.3,
                              labels=sparse_labels,
                              xlabelsize="small",
                              ylabelsize="small",
                              colormap=cm.autumn_r)
    for ax in axes:
        # adding reference would increase x range too much for these
        if sigma_type != "xsmall": 
            # draw reference for true param value
            ax.axvline(x=ref_val, color="k", linewidth=.8, zorder=10000)
    
    last_ax = axes[-1]
    last_ax.axhline(y=last_ax.get_ylim()[0], linewidth=1, color='k') # draw x axis
    if sigma_type == "xsmall": 
        ticks = last_ax.get_xticks()
        xlabels = ['{:.2e}'.format(n) for n in ticks]

    last_ax.tick_params(axis='both', which='minor', bottom=False, labelbottom=False) 
    last_ax.tick_params(axis='x', which='major', length=3, width=.8)
    
    if xlim is not None:
        ax.set_xticks(ticks)
        ax.set_xticklabels(xlabels)
    
    fig.savefig("figures/" + outfile + "_v2.{}".format(out_format), 
                format=out_format, dpi="figure", bbox_inches="tight")
    plt.close()


# -

filename = "data/plots_v2/One_sided_sampling_BMP,_beast_without_extra_samples_sigmaX.txt"
make_ridgeplot(filename)

#filenames = !ls data/plots_v2/*.txt
filenames = [
    "data/plots_v2/Diagonal_Sampling_BMP,_phyrex_rootX.txt",
    "data/plots_v2/Diagonal_Sampling_BMP,_phyrex_rootY.txt",
    #"data/plots_v2/Diagonal_Sampling_BMP,_phyrex_SigmaFormula.txt"
]
for filename in filenames:
    make_ridgeplot(filename)

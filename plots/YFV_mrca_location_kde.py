# -*- coding: utf-8 -*-
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
import pandas as pd
# No warnings about setting value on copy of slice
pd.options.mode.chained_assignment = None

# Display up to 60 columns of a dataframe
pd.set_option('display.max_columns', 60)

from io import StringIO
from dendropy import Tree, TaxonNamespace

import geopandas as gpd
from shapely.geometry import Point
import contextily as ctx
import geoplot as gplt
import geoplot.crs as gcrs

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# # %matplotlib inline

# Set default font size
plt.rcParams['font.size'] = 24

# Internal ipython tool for setting figure size
from IPython.core.pylabtools import figsize

import seaborn as sb
# Set default font size
sb.set(font_scale = .8)
custom_style = {'axes.labelcolor': 'black',
                'xtick.color': 'black',
                'ytick.color': 'black'}
sb.set_style("white", rc=custom_style)

from itertools import chain

from matplotlib.collections import LineCollection
from matplotlib import markers
from matplotlib.path import Path

import numpy as np

from ete3 import Tree as eteTree
from ete3 import TreeStyle, NodeStyle

full_filename = "data/YFV.trees"
downsamp_filename = "data/YFV_downSamp.trees"
full_tip_count = 65
downsamp_tip_count = 47


# +
def read_ids(filename, tip_count):
    # read ids for trees sampled from an MCMC process from a nexus file
    tree_ids_str = !tail  -n +{tip_count * 2 + 13} {filename} | cut -d " "  -f2  | grep "STATE_"| tr "\n" " "
    tree_ids = tree_ids_str[0].split()
    return tree_ids

def read_trees(filename, tip_count):
    # read the nexus strings
    tree_str = !tail  -n +{tip_count * 2 + 13} {filename} | grep "STATE_"  | cut -d " " -f3- 
    return tree_str

def read_header(filename, tip_count):
    header = !head -n +{tip_count * 2 + 12} {filename} # read the nexus header
    return header

tree_ids = read_ids(full_filename, full_tip_count)
full_trees = read_trees(full_filename, full_tip_count)
downsamp_trees = read_trees(downsamp_filename, downsamp_tip_count)

full_header = read_header(full_filename, full_tip_count)
downsamp_header = read_header(downsamp_filename, downsamp_tip_count)


# +
def write_and_read_nexus(filename, header, tree_id, tree_str):
    tns = TaxonNamespace(is_case_sensitive=True)
    # write a temp file containing  tree
    with open(filename, "w") as f:
        for line in header + ["tree " + tree_id + " " + tree_str]:
            f.write(line + "\n");
    # read tree as dendropy tree
    tree = Tree.get(path=filename, schema="nexus",
                    taxon_namespace=tns, case_sensitive_taxon_labels=True, 
                    suppress_internal_node_taxa=False)
    return tree

def get_locations(tree_id, full_str, downsamp_str):
    full_tree = write_and_read_nexus('full.tmp', full_header, tree_id, full_str)
    downsamp_tree = write_and_read_nexus('downsamp.tmp', downsamp_header, tree_id, downsamp_str)
    
    downsamp_tip_ids = list([n.taxon.label for n in downsamp_tree.leaf_node_iter()])
    
    # get the root node's inferred location in the full tree
    full_tree_root_xy = full_tree.seed_node.annotations["location"].value
    
    # get the downsample tree root node location
    downsamp_tree_root_xy = downsamp_tree.seed_node.annotations["location"].value

    output = [full_tree_root_xy, downsamp_tree_root_xy]
    flat_output = [tree_id] + [item for sublist in output for item in sublist]
    return flat_output


# +
arr = []
for tree_id, full_str, downsamp_str in zip(tree_ids, full_trees, downsamp_trees):
    arr.append(get_locations(tree_id, full_str, downsamp_str))
mrca_xy_df = pd.DataFrame(arr, columns=["tree_id", 
                                        "full_tree_root_x", 
                                        "full_tree_root_y", 
                                        "downsamp_tree_root_x", 
                                        "downsamp_tree_root_y"])
# reshape df 
df1 = mrca_xy_df[["tree_id", "full_tree_root_x", "full_tree_root_y"]]
df1["type"] = ["full"]*len(df1)
df1.columns = ["tree_id", "Latitude", "Longitude", "type"]
df2 = mrca_xy_df[["tree_id", "downsamp_tree_root_x", "downsamp_tree_root_y"]]
df2["type"] = ["downsamp"]*len(df2)
df2.columns = ["tree_id", "Latitude", "Longitude", "type"]
mrca_xy_reshaped_df = pd.concat([df1, df2], ignore_index=True)

# make coords numeric
mrca_xy_reshaped_df[["Longitude", "Latitude"]] = mrca_xy_reshaped_df[["Longitude", "Latitude"]].\
apply(pd.to_numeric, errors='coerce')
# drop illegal coordinates
mrca_xy_reshaped_df = mrca_xy_reshaped_df[(mrca_xy_reshaped_df["Latitude"] > -90) & \
                                          (mrca_xy_reshaped_df["Latitude"] < 90) & \
                                         (mrca_xy_reshaped_df["Longitude"] > -180) & \
                                          (mrca_xy_reshaped_df["Longitude"] < 180)]

# covert to geodataframe, generate Point objects
gdf = gpd.GeoDataFrame(mrca_xy_reshaped_df)
gdf["geometry"] = [Point(x, y) for x,y in zip(mrca_xy_reshaped_df.Longitude, mrca_xy_reshaped_df.Latitude)]
# specify epsg coordinate system for lon, lat coods
# this is the WGS 84 -- WGS84 - World Geodetic System 1984, used in GPS
gdf.crs = {'init': 'epsg:4326'}
# transform from long, lat to projection coordinates (not necessary for geoplot)
# this is  WGS 84 / Pseudo-Mercator - Spherical Mercator(Google Maps, OpenStreetMap) 
gdf = gdf.to_crs(epsg=3857)
# get back transformed xy coords for plotting
# in a projection system
mrca_xy_crs_df = pd.DataFrame(gdf["geometry"].\
apply(lambda point: [float(s.replace("(", "").replace(")", "")) \
for s in str(point).split()[1:]]).tolist(), columns=["Longitude", "Latitude"])
mrca_xy_crs_df[["tree_id", "type"]] = gdf[["tree_id", "type"]] 
mrca_xy_crs_df.head();
mrca_xy_crs_df.to_csv("data/YFV_ALLsubset_vs_fullmrca_location_crs.tsv", sep="\t")


# +
# get all observed tip node locations
def get_tip_locations(tree_id, full_str, downsamp_str):
    full_tree = write_and_read_nexus('full.tmp', full_header, tree_id, full_str)
    downsamp_tree = write_and_read_nexus('downsamp.tmp', downsamp_header, tree_id, downsamp_str)
    
    full_tip_ids = list([n.taxon.label for n in full_tree.leaf_node_iter()]) 
    downsamp_tip_ids = list([n.taxon.label for n in downsamp_tree.leaf_node_iter()])
    output = []
    
    for tip_id in full_tip_ids:
        if tip_id in downsamp_tip_ids:
            tip_type = "downsamp"            
        else:
            tip_type = "removed"
        tip_node = full_tree.find_node_with_taxon_label(tip_id)
        tip_x, tip_y = tip_node.annotations["location"].value
        output.append([tip_id, tip_type, tip_x, tip_y])
    return output


for tree_id, full_str, downsamp_str in zip(tree_ids[0], full_trees, downsamp_trees):
    arr = get_tip_locations(tree_id, full_str, downsamp_str)
    
# order in node annotation is lat;lon
tip_df = pd.DataFrame(arr, columns=["tip_id", "tip_type", "Latitude", "Longitude"])
tip_df[["Longitude", "Latitude"]] = tip_df[["Longitude", "Latitude"]].\
apply(pd.to_numeric)

# covert to geodataframe, generate Point objects
tip_gdf = gpd.GeoDataFrame(tip_df)
# Point object is specified with lon;lat
tip_gdf["geometry"] = [Point(x, y) for x,y in zip(tip_df.Longitude, tip_df.Latitude)]
# specify epsg coordinate system for lon, lat coods
tip_gdf.crs = {'init': 'epsg:4326'}

# transform from lon,lat to projection coordinates (not necessary for geoplot)
tip_gdf = tip_gdf.to_crs(epsg=3857)

# get back transformed xy coords for plotting
# in a projection system
tip_crs_df = pd.DataFrame(tip_gdf["geometry"].\
apply(lambda point: [float(s.replace("(", "").replace(")", "")) \
for s in str(point).split()[1:]]).tolist(), columns=["Longitude", "Latitude"])
tip_crs_df[["tip_id", "tip_type"]] = tip_gdf[["tip_id", "tip_type"]]
tip_crs_df.to_csv("data/YFV_ALLsubset_vs_fulltip_location.tsv", sep="\t")

# +
tip_crs_df = pd.read_csv("data/YFV_ALLsubset_vs_fulltip_location.tsv", sep="\t")

# cities from which samples have been removed for downsample analysis,
# truncate city names b/c some are misspelled in labels
downsamp_cities = ["Ladainha",
"TeofiloOtoni",
"DomingosMartins",
"Caria",
"NovoCruzeiro",
"Itamba"]

def match_city(s):
    for c in downsamp_cities:
        if c in s:
            return c
        
# get the tip node's location city name from the label      
tip_crs_df["downsamp_city"] = tip_crs_df.tip_id.apply(match_city)
# marker size = how many samples from each city
marker_size = 20
city_count_dict = tip_crs_df.tip_id.apply(match_city).value_counts().to_dict()
tip_crs_df["downsamp_city_count"] = tip_crs_df["downsamp_city"].\
apply(lambda x: marker_size * city_count_dict.get(x, 1.))

# save df
tip_crs_df.to_csv("data/YFV_ALLsubset_vs_fulltip_location.tsv", sep="\t")


# +
def calc_summary_tree(mcmc_tree_filename):
    # create a summary tree from allsampled trees using Bsouth's treeannotator
    !//anaconda3/envs/gpd/bin/treeannotator -burninTrees 1000 {mcmc_tree_filename} {mcmc_tree_filename}_summary.tree

def summary_to_nw_str(mcmc_tree_filename):
    calc_summary_tree(mcmc_tree_filename)
    # convert summary nexus tree to newick for ete3
    tns = TaxonNamespace(is_case_sensitive=True)
    filename = mcmc_tree_filename + "_summary.tree"
    dp_tree = Tree.get(path=filename, 
                       schema="nexus",
                       taxon_namespace=tns,
                       case_sensitive_taxon_labels=True,
                       suppress_internal_node_taxa=False)

    # drop all annotations and illegal characters
    return dp_tree.as_string('newick', suppress_annotations=True)[5:].rstrip("\n")                

# compute summary tree, write as newick string
full_nw_str = summary_to_nw_str(full_filename)
downsamp_nw_str = summary_to_nw_str(downsamp_filename)
with open("full_nw_str", "w") as f:
    f.write(full_nw_str)
with open("downsamp_nw_str", "w") as f:
    f.write(downsamp_nw_str)
# -

mrca_xy_crs_df

# +
# read inferred and observed postion dfs
mrca_xy_crs_df = pd.read_csv("data/YFV_ALLsubset_vs_fullmrca_location_crs.tsv", sep="\t")
tip_crs_df = pd.read_csv("data/YFV_ALLsubset_vs_fulltip_location.tsv", sep="\t")

# read the full tree
ete_tree = eteTree("./full_nw_str")
ts = TreeStyle()
ete_tree.set_style(ts) #set deafult node and tree style

style_dict = {
    "downsamp" : ("blue", "square", "b", "s", sb.light_palette("blue", as_cmap=True)),
    "removed" : ("orange", "circle", "o", "o", sb.light_palette("orange", as_cmap=True)),
}

# create custom node styles
nstyle_downsamp = NodeStyle()
nstyle_downsamp['fgcolor'] = style_dict["downsamp"][0]
nstyle_downsamp['size'] = 4
nstyle_downsamp['shape'] = style_dict["downsamp"][1]

nstyle_removed = NodeStyle()
nstyle_removed['fgcolor'] = style_dict["removed"][0]
nstyle_removed['size'] = 4
nstyle_removed['shape'] = style_dict["removed"][1]

nstyle_root = NodeStyle()
nstyle_root['fgcolor'] = 'red'
nstyle_root['size'] = 8
nstyle_root['shape'] = '^'

nstyle_int = NodeStyle()
nstyle_int['fgcolor'] = 'black'
nstyle_int['size'] = 2
nstyle_int['shape'] = "circle"

# get the full set of tip node names
ete_full_tips = [n for n in ete_tree.iter_leaf_names()]
# map tip names to numbers for visualization (full labels are too long)
name_to_num = dict(zip(ete_full_tips, range(len(ete_full_tips))))
root_node = ete_tree.get_tree_root()
# get subset of tip node names from the downsampled tree
ete_subset_tree = eteTree("./downsamp_nw_str")
ete_subset_tips =  [n for n in ete_subset_tree.iter_leaf_names()]

# assign custom styles to nodes in full tree
for n in ete_tree.traverse():
    if n.is_leaf() and n.name in ete_subset_tips:
        n.set_style(nstyle_downsamp)
    elif n.is_leaf():
        n.set_style(nstyle_removed)
    elif n == root_node:
        n.set_style(nstyle_root)
    else:
        n.set_style(nstyle_int)

# create a custom grid of subplots
fig = plt.figure(figsize=(10, 12))
widths = [1, 1]
heights = [30, 1, 30]
gs = fig.add_gridspec(ncols=2, nrows=3,
                      width_ratios=widths,
                      height_ratios=heights,
                     hspace=0.01)
tree_ax = fig.add_subplot(gs[0, :])
map_ax = fig.add_subplot(gs[2, :])
cax_green = fig.add_subplot(gs[1, 0])
cax_blue = fig.add_subplot(gs[1, 1])

def round_sig(x, sig=2):
    return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)

def to_coord(x, y, xmin, xmax, ymin, ymax, plt_xmin, plt_ymin, plt_width, plt_height):
    x = (x - xmin) / (xmax - xmin) * plt_width  + plt_xmin
    y = (y - ymin) / (ymax - ymin) * plt_height + plt_ymin
    return x, y

cstyle = NodeStyle()
align_names=False
name_offset=None
max_dist=None
font_size=3

def __draw_edge_nm(c, y):
    h = node_pos[c]
    vlinec.append(((h, y), (h, y + c.dist)))
    vlines.append(cstyle)
    return (h, y + c.dist)

def __draw_edge_md(c, x):
    h = node_pos[c]
    if c in cut_edge:
        offset = max_x / 600.
        hlinec.append(((x, h), (x + c.dist / 2 - offset, h)))
        hlines.append(cstyle)
        hlinec.append(((x + c.dist / 2 + offset, h), (x + c.dist, h)))
        hlines.append(cstyle)
        hlinec.append(((x + c.dist / 2, h - 0.05), (x + c.dist / 2 - 2 * offset, h + 0.05)))
        hlines.append(cstyle)
        hlinec.append(((x + c.dist / 2 + 2 * offset, h - 0.05), (x + c.dist / 2, h + 0.05)))
        hlines.append(cstyle)
        tree_ax.text(x + c.dist / 2, h - 0.07, '+%g' % max_dist, va='top', 
                 ha='center', rotation=90, size=2. * font_size / 3)
    else:
        hlinec.append(((x, h), (x + c.dist, h)))
        hlines.append(cstyle)
    return (x + c.dist, h)

__draw_edge = __draw_edge_nm if max_dist is None else __draw_edge_md

vlinec = []
vlines = []
hlinec = []
hlines = []
nodes = []
nodex = []
nodey = []
ali_lines = []

# to align leaf names
tree = ete_tree
max_y = max(n.get_distance(tree) for n in tree.iter_leaves())
# extra pad
max_y = max_y + max_y * .1

coords = {}
# position on x axis
node_pos = dict((n2, i) for i, n2 in enumerate(tree.get_leaves()[::-1]))
node_list = tree.iter_descendants(strategy='postorder')
node_list = chain(node_list, [tree])

# reduce branch length
cut_edge = set()
if max_dist is not None:
    for n in tree.iter_descendants():
        if n.dist > max_dist:
            n.dist -= max_dist
            cut_edge.add(n)

if name_offset is None:
    name_offset = max_y / 30.
# draw tree
for n in node_list:
    style = n._get_style()
    # distance from root on y-axis
    y = __builtin__.sum(n2.dist for n2 in n.iter_ancestors()) + n.dist
    if n.is_leaf():
        x = node_pos[n]
        tree_ax.text(x,
                     y + name_offset,
                     name_to_num[n.name],
                     ha="center",
                     va='center',
                     rotation=90,
                     size=8)
    else:
        x = np.mean([node_pos[n2] for n2 in n.children])
        node_pos[n] = x

        # draw horizontal line
        hlinec.append(((node_pos[n.children[0]], y), (node_pos[n.children[-1]], y)))
        hlines.append(style)

        # draw vertical lines
        for child in n.children:
            cstyle = child._get_style()
            coords[child] = __draw_edge(child, y)
    nodes.append(style)
    nodex.append(x)
    nodey.append(y)

# draw root
__draw_edge(tree, 0)

lstyles = ['-', '--', ':']
hline_col = LineCollection(hlinec, colors=[l['hz_line_color'] for l in hlines], 
                          linestyle=[lstyles[l['hz_line_type']] for l in hlines],
                         linewidth=[(l['hz_line_width'] + 1.) / 2 for l in hlines])
vline_col = LineCollection(vlinec, colors=[l['vt_line_color'] for l in vlines], 
                         linestyle=[lstyles[l['vt_line_type']] for l in vlines],
                          linewidth=[(l['vt_line_width'] + 1.) / 2 for l in vlines])
ali_line_col = LineCollection(ali_lines, colors='k')

tree_ax.add_collection(hline_col)
tree_ax.add_collection(vline_col)
tree_ax.add_collection(ali_line_col)

nshapes = dict((('circle', 'o'), ('square', 's'), ('sphere', 'o')))
shapes = set(n['shape'] for n in nodes)
for shape in shapes:
    indexes = [i for i, n in enumerate(nodes) if n['shape'] == shape]
    scat = tree_ax.scatter([nodex[i] for i in indexes], 
                       [nodey[i] for i in indexes], 
                       s=0, marker=nshapes.get(shape, shape))
    scat.set_sizes([(nodes[i]['size'])**2 / 2 for i in indexes])
    scat.set_color([nodes[i]['fgcolor'] for i in indexes])
    scat.set_zorder(10)

# custom y axis
tree_ax.set_frame_on(False)    
tree_ax.get_xaxis().set_visible(False)
tree_ax.invert_yaxis()
xmin, xmax = tree_ax.get_xaxis().get_view_interval()
ymin, ymax = tree_ax.get_yaxis().get_view_interval()
tree_ax.add_artist(plt.Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=2))
tree_ax.tick_params(labelsize=10, length=3, left=True)
tree_ax.set_ylabel("Years")

#from http://www.naturalearthdata.com/downloads/10m-cultural-vectors/
brazil_states = gpd.read_file("data/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")
brazil_states = brazil_states[brazil_states["name_en"].isin(["Rio de Janeiro", "Minas Gerais","EspÃ\xadrito Santo"])].to_crs(epsg=3857)

brazil_states.plot(linewidth=1,
        edgecolor='white',
        facecolor='lightgray',
        ax=map_ax);


# drop ~10% of sampled trees as MCMC burn-in
drop_num = 1000
# inferred root location in full tree
full_df = mrca_xy_crs_df[mrca_xy_crs_df["type"] == "full"][drop_num:] 
# root location in downsampled tree
subset_df = mrca_xy_crs_df[mrca_xy_crs_df["type"] == "downsamp"][drop_num:]   

# draw 2D KDE plot
# have to flip lat;lon -> lon;lat when plotting on X;Y axis system
sb.kdeplot(full_df.Longitude, 
           full_df.Latitude, 
           cmap=sb.light_palette("green", as_cmap=True), 
           shade=True, 
           shade_lowest=False, 
           alpha=.7, 
           cbar=True, 
           cbar_ax=cax_green,
           cbar_kws={"orientation": "horizontal"},
           zorder=5,
           ax=map_ax);

sb.kdeplot(subset_df.Longitude,
           subset_df.Latitude,
           cmap=sb.light_palette("blue", as_cmap=True),
           shade=True, 
           shade_lowest=False,  
           alpha=.7, 
           cbar=True, 
           cbar_ax=cax_blue,
           cbar_kws={"orientation": "horizontal"},
           ax=map_ax);

# set colormap axes style
cax_green.xaxis.set_ticks_position('top')
cax_green.xaxis.set_label_position('top')
cax_green.tick_params(labelsize=10, length=3, left=True, labelleft=True, 
                    right=False, labelright=False)
cax_green.set_xlabel("root location (KDE) in full analysis")
cax_blue.xaxis.set_ticks_position('top')
cax_blue.xaxis.set_label_position('top')
cax_blue.tick_params(labelsize=10, length=3, right=True, labelright=True)
cax_blue.set_xlabel("root location (KDE) in downsampled analysis")

# show tip location w/t scatterplot over geo map
downsamp_df = tip_crs_df[tip_crs_df["tip_type"] == "downsamp"]
downsamp_df.downsamp_city_count = [20] * len(downsamp_df)
sb.scatterplot(x=downsamp_df.Longitude, 
               y=downsamp_df.Latitude,
               hue=downsamp_df.tip_type,
               palette = {"downsamp" : "blue"},
               size=downsamp_df.downsamp_city_count,
               markers = {"downsamp" : style_dict["downsamp"][3]},
               #legend = False,
               ax=map_ax);

# subset markers on top
removed_df = tip_crs_df[tip_crs_df["tip_type"] == "removed"]
sb.scatterplot(x=removed_df.Longitude, 
               y=removed_df.Latitude,
               hue=removed_df.tip_type,
               palette = {"removed" : "orange"},
               size=removed_df.downsamp_city_count,
               markers = {"removed" : style_dict["removed"][3]},
               #legend = False,
               ax=map_ax);


# annotate scatter with tip_idx label
name_to_num = dict([(k.replace("'", ""), name_to_num[k]) for k in name_to_num.keys()])
# generic annotation
generic_df =  tip_crs_df[tip_crs_df["downsamp_city_count"] ==  marker_size]
generic_labels = [name_to_num[l] for l in generic_df["tip_id"]]
X=generic_df.Longitude.values
Y=generic_df.Latitude.values
for i, label in enumerate(generic_labels):
    map_ax.annotate(label, 
                   (X[i], 
                     Y[i]),
                  fontsize=6)

# downsampled cities annotation with arrow
for city in downsamp_cities:
    city_df = tip_crs_df[tip_crs_df["downsamp_city"] ==  city]
    downsamp_city_df = city_df[city_df["tip_type"] == "downsamp"]
    removed_city_df = city_df[city_df["tip_type"] == "removed"]
    downsamp_city_labels = [name_to_num[l] for l in downsamp_city_df["tip_id"]]
    removed_city_labels = [name_to_num[l] for l in removed_city_df["tip_id"]]
    # list all labels for city, put parentheses around removed samples
    annot_str = "(" + ",".join(map(str, removed_city_labels)) + "), " +\
    ",".join(map(str, downsamp_city_labels))
    # use the first tip coords
    x = city_df.Longitude.values[0]
    y = city_df.Latitude.values[0]
    map_ax.annotate(annot_str, 
                    xy=(x, y),
                    fontsize=6, 
                    xytext=(50,-10), 
                    textcoords='offset points', 
                    ha='center',
                    va='bottom', 
                    arrowprops=dict(arrowstyle='-', 
                                    color='k', 
                                    linewidth=.5))

# custom legend
handles, labels = map_ax.get_legend_handles_labels()
map_ax.legend([handles[i] for i in [1,6]], 
              ["retained","removed"], 
              loc='lower right', 
              title="Tip subset")
#[-14201398.96929843, -7133189.60868023, 2712783.1402848554, 6513955.476226102]
map_ax.set_xlim(left=-5400000, right=-4250000)
map_ax.set_ylim(bottom=-2700000, top=-1600000) 
map_ax.set_axis_off()
gs.tight_layout(fig)
fig.show();
fig.savefig("figures/YFV_mixed1_vs_fullmrca_location_kde_tree.pdf", 
            format="pdf", dpi="figure", bbox_inches="tight")
fig.savefig("figures/YFV_mixed1_vs_fullmrca_location_kde_tree.svg", 
            format="svg", dpi="figure", bbox_inches="tight")

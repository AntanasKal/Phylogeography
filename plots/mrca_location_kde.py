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

full_filename = "data/WNV1.trees"
west_filename = "data/WNV_west1.trees"
#east_filename = "data/WNV_east1.trees"
mixed_filename = "data/WNV_west+20east.trees"
mixed_alt_filename = "data/WNV_west+52east.trees"
full_tip_count = 104
west_tip_count = 52
mixed_tip_count = 72
mixed_alt_tip_count = full_tip_count


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
    # read the nexus header
    header = !head -n +{tip_count * 2 + 12} {filename}
    return header

tree_ids = read_ids(full_filename, full_tip_count)
full_trees = read_trees(full_filename, full_tip_count)
west_trees = read_trees(west_filename, west_tip_count)
#east_trees = read_trees(east_filename, east_tip_count)
mixed_trees = read_trees(mixed_filename, mixed_tip_count)
mixed_alt_trees = read_trees(mixed_alt_filename, mixed_alt_tip_count)

full_header = read_header(full_filename, full_tip_count)
west_header = read_header(west_filename, west_tip_count)
#east_header = read_header(east_filename, east_tip_count)
mixed_header = read_header(mixed_filename, mixed_tip_count)
mixed_alt_header = read_header(mixed_alt_filename, mixed_alt_tip_count)
# -

mixed_alt_tip_count


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

def get_locations(tree_id, full_str, west_str, #east_str, 
                  mixed_str, mixed_alt_str):
    full_tree = write_and_read_nexus('full.tmp', full_header, tree_id, full_str)
    west_tree = write_and_read_nexus('west.tmp', west_header, tree_id, west_str)
    #east_tree = write_and_read_nexus('east.tmp', east_header, tree_id, east_str)
    mixed_tree = write_and_read_nexus('mixed.tmp', mixed_header, tree_id, mixed_str)
    mixed_alt_tree = write_and_read_nexus('mixed_alt.tmp', mixed_alt_header, tree_id, mixed_alt_str)
    west_tip_ids = list([n.taxon.label for n in west_tree.leaf_node_iter()])
    #east_tip_ids = list([n.taxon.label for n in east_tree.leaf_node_iter()])
    mixed_tip_ids = list([n.taxon.label for n in mixed_tree.leaf_node_iter()])
    mixed_alt_tip_ids = list([n.taxon.label for n in mixed_alt_tree.leaf_node_iter()])
    
    # get the most recent common ancestor (mrca) for the tip subset
    west_mrca_node = full_tree.mrca(taxon_labels=west_tip_ids)
    #east_mrca_node = full_tree.mrca(taxon_labels=east_tip_ids)
    # for the "mixed" case we compare against the "west" subset mrca location
    mixed_mrca_node = west_mrca_node
    mixed_alt_mrca_node = west_mrca_node
    
    # get the mrca inferred location in the full tree
    west_full_tree_xy = west_mrca_node.annotations["location"].value
    #east_full_tree_xy = east_mrca_node.annotations["location"].value
    mixed_full_tree_xy = mixed_mrca_node.annotations["location"].value
    mixed_alt_full_tree_xy = mixed_alt_mrca_node.annotations["location"].value
    
    # get the subset tree root node location
    west_tree_xy = west_tree.seed_node.annotations["location"].value
    #east_tree_xy = east_tree.seed_node.annotations["location"].value
    # for the "mixed" case get the west subset mrca location in the mixed subset tree
    mixed_mixed_mrca_node = mixed_tree.mrca(taxon_labels=west_tip_ids)
    mixed_mixed_alt_mrca_node = mixed_alt_tree.mrca(taxon_labels=west_tip_ids)
    mixed_tree_xy = mixed_mixed_mrca_node.annotations["location"].value
    mixed_alt_tree_xy = mixed_mixed_alt_mrca_node.annotations["location"].value

    output = [west_full_tree_xy,
              west_tree_xy,
              #east_full_tree_xy,
              #east_tree_xy,
             mixed_full_tree_xy,
             mixed_tree_xy,
             mixed_alt_full_tree_xy,
             mixed_alt_tree_xy]
    flat_output = [tree_id] + [item for sublist in output for item in sublist]
    return flat_output


# +
arr = []
for tree_id, full_str, west_str, mixed_str, mixed_alt_str in zip(tree_ids, 
                 full_trees, 
                 west_trees, 
                 #east_trees, 
                 mixed_trees,
                 mixed_alt_trees):
    arr.append(get_locations(tree_id, 
                             full_str, 
                             west_str, 
                             #east_str, 
                             mixed_str,
                            mixed_alt_str))
mrca_xy_df = pd.DataFrame(arr, columns=["tree_id", 
                                        "west_full_tree_x", 
                                        "west_full_tree_y", 
                                        "west_tree_x", 
                                        "west_tree_y",
                                       #"east_full_tree_x", 
                                        #"east_full_tree_y", 
                                        #"east_tree_x", 
                                        #"east_tree_y",
                                       "mixed_full_tree_x", 
                                        "mixed_full_tree_y", 
                                        "mixed_tree_x", 
                                        "mixed_tree_y",
                                       "mixed_alt_full_tree_x", 
                                        "mixed_alt_full_tree_y", 
                                        "mixed_alt_tree_x", 
                                        "mixed_alt_tree_y"])
# reshape df 
df1 = mrca_xy_df[["tree_id", "west_full_tree_x", "west_full_tree_y"]]
df1["type"] = ["west_full"]*len(df1)
df1.columns = ["tree_id", "Latitude", "Longitude", "type"]
df2 = mrca_xy_df[["tree_id", "west_tree_x", "west_tree_y"]]
df2["type"] = ["west"]*len(df2)
df2.columns = ["tree_id", "Latitude", "Longitude", "type"]
#df3 = mrca_xy_df[["tree_id", "east_full_tree_x", "east_full_tree_y"]]
#df3["type"] = ["east_full"]*len(df3)
#df3.columns = ["tree_id", "Latitude", "Longitude", "type"]
#df4 = mrca_xy_df[["tree_id", "east_tree_x", "east_tree_y"]]
#df4["type"] = ["east"]*len(df4)
#df4.columns = ["tree_id", "Latitude", "Longitude", "type"]
df5 = mrca_xy_df[["tree_id", "mixed_full_tree_x", "mixed_full_tree_y"]]
df5["type"] = ["mixed_full"]*len(df5)
df5.columns = ["tree_id", "Latitude", "Longitude", "type"]
df6 = mrca_xy_df[["tree_id", "mixed_tree_x", "mixed_tree_y"]]
df6["type"] = ["mixed"]*len(df6)
df6.columns = ["tree_id", "Latitude", "Longitude", "type"]
df7 = mrca_xy_df[["tree_id", "mixed_alt_full_tree_x", "mixed_alt_full_tree_y"]]
df7["type"] = ["mixed_alt_full"]*len(df7)
df7.columns = ["tree_id", "Latitude", "Longitude", "type"]
df8 = mrca_xy_df[["tree_id", "mixed_alt_tree_x", "mixed_alt_tree_y"]]
df8["type"] = ["mixed_alt"]*len(df8)
df8.columns = ["tree_id", "Latitude", "Longitude", "type"]
mrca_xy_reshaped_df = pd.concat([df1, 
                                 df2, 
                                 #df3, 
                                 #df4, 
                                 df5, 
                                 df6,
                                df7,
                                df8], ignore_index=True)

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
mrca_xy_crs_df.to_csv("data/WNV_ALLsubset_vs_fullmrca_location_crs.tsv", sep="\t")


# +
# get all the tip locations
def get_tip_locations(tree_id, full_str, west_str, mixed_str, mixed_alt_str):
    full_tree = write_and_read_nexus('full.tmp', full_header, tree_id, full_str)
    west_tree = write_and_read_nexus('west.tmp', west_header, tree_id, west_str)
    mixed_tree = write_and_read_nexus('mixed.tmp', mixed_header, tree_id, mixed_str)
    mixed_alt_tree = write_and_read_nexus('mixed_alt.tmp', mixed_alt_header, tree_id, mixed_alt_str)
    full_tip_ids = list([n.taxon.label for n in full_tree.leaf_node_iter()]) 
    west_tip_ids = list([n.taxon.label for n in west_tree.leaf_node_iter()])
    mixed_tip_ids = list([n.taxon.label for n in mixed_tree.leaf_node_iter()])
    mixed_alt_tip_ids = list([n.taxon.label for n in mixed_alt_tree.leaf_node_iter()])
    output = []
    
    for tip_id in full_tip_ids:
        if tip_id in west_tip_ids:
            tip_type = "west"            
        else:
            tip_type = "east"
        tip_node = full_tree.find_node_with_taxon_label(tip_id)
        tip_x, tip_y = tip_node.annotations["location"].value
        output.append([tip_id, tip_type, tip_x, tip_y])
        
    for tip_id in full_tip_ids:
        if tip_id in mixed_tip_ids:
            tip_type = "mixed"
        else:
            tip_type = "non_mixed"
        tip_node = full_tree.find_node_with_taxon_label(tip_id)
        tip_x, tip_y = tip_node.annotations["location"].value
        output.append([tip_id, tip_type, tip_x, tip_y])
        
    for tip_id in full_tip_ids:
        if tip_id in mixed_alt_tip_ids:
            tip_type = "mixed_alt"
        else:
            tip_type = "non_mixed_alt"
        tip_node = full_tree.find_node_with_taxon_label(tip_id)
        tip_x, tip_y = tip_node.annotations["location"].value
        output.append([tip_id, tip_type, tip_x, tip_y])
    return output

for tree_id, full_str, west_str, mixed_str, mixed_alt_str in zip(tree_ids[0],
                                                                 full_trees,
                                                                 west_trees,
                                                                 mixed_trees,
                                                                mixed_alt_trees):
    arr = get_tip_locations(tree_id, full_str, west_str, mixed_str, mixed_alt_str)
    
# order in node annotation is lat;lon
tip_df = pd.DataFrame(arr, columns=["tip_id", "tip_type", "Latitude", "Longitude"])
tip_df[["Longitude", "Latitude"]] = tip_df[["Longitude", "Latitude"]].\
apply(pd.to_numeric, errors='coerce')

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
tip_crs_df.to_csv("data/WNV_ALLsubset_vs_fulltip_location.tsv", sep="\t")


# +
def calc_summary_tree(mcmc_tree_filename):
    # create a summary tree from allsampled trees using BEAST's treeannotator
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


# +
full_nw_str = summary_to_nw_str(full_filename)
west_nw_str = summary_to_nw_str(west_filename)
#east_nw_str = summary_to_nw_str(east_filename)
mixed_nw_str = summary_to_nw_str(mixed_filename)
mixed_alt_nw_str = summary_to_nw_str(mixed_alt_filename)


mrca_xy_crs_df = pd.read_csv("data/WNV_ALLsubset_vs_fullmrca_location_crs.tsv", sep="\t")
tip_crs_df = pd.read_csv("data/WNV_ALLsubset_vs_fulltip_location.tsv", sep="\t")

# +
# custom for "mixed" cases
# plot_map_figure("mixed", "non_mixed", full_nw_str, mixed_nw_str, 
#                mrca_str=west_nw_str, mrca_label="west");
ete_tree = eteTree(full_nw_str)
ts = TreeStyle()

style_dict = {
    "west" : ("red", "square", "r", "s", sb.light_palette("red", as_cmap=True)),
    "east" : ("green", "circle", "g", "o", sb.light_palette("green", as_cmap=True)),
    "non_mixed" : ("green", "circle", "b", "o", sb.light_palette("magenta", as_cmap=True))
}

# set node style
nstyle_west = NodeStyle()
nstyle_west['fgcolor'] = style_dict["west"][0]
nstyle_west['size'] = 4
nstyle_west['shape'] = style_dict["west"][1]

nstyle_east = NodeStyle()
nstyle_east['fgcolor'] = style_dict["east"][0]
nstyle_east['size'] = 4
nstyle_east['shape'] = style_dict["east"][1]

nstyle_cmpl = NodeStyle()
nstyle_cmpl['fgcolor'] = style_dict["non_mixed"][0]
nstyle_cmpl["shape"] = style_dict["non_mixed"][1]
nstyle_cmpl['size'] = 4

nstyle_mrca = NodeStyle()
nstyle_mrca['fgcolor'] = 'red'
nstyle_mrca['size'] = 8
nstyle_mrca['shape'] = '^'

nstyle_int = NodeStyle()
nstyle_int['fgcolor'] = 'black'
nstyle_int['size'] = 2
nstyle_int['shape'] = "circle"

subset_nw_str = mixed_alt_nw_str
ete_full_tips = [n for n in ete_tree.iter_leaf_names()]
name_to_num = dict(zip(ete_full_tips, range(len(ete_full_tips))))
ete_subset_tree = eteTree(subset_nw_str)
ete_subset_tips =  [n for n in ete_subset_tree.iter_leaf_names()]
ete_tree.set_style(ts)

mrca_ete_tree = eteTree(west_nw_str)
mrca_subset = [n for n in mrca_ete_tree.iter_leaf_names()]
mrca = ete_tree.get_common_ancestor(*mrca_subset)

for n in ete_tree.traverse():
    if n.is_leaf() and n.name in ete_subset_tips:
        if n.name in mrca_subset:
            n.set_style(nstyle_west)
        else:
            n.set_style(nstyle_east)
    elif n.is_leaf():
        n.set_style(nstyle_cmpl)
    elif n == mrca:
        n.set_style(nstyle_mrca)
    else:
        n.set_style(nstyle_int)

name_to_num = dict(zip(ete_full_tips, range(len(ete_full_tips))))
fig = plt.figure(figsize=(10, 12))
widths = [1, 1, 1]
heights = [30, 1, 30]
gs = fig.add_gridspec(ncols=3, nrows=3,
                      width_ratios=widths,
                      height_ratios=heights,
                     hspace=0.01)
tree_ax = fig.add_subplot(gs[0, :])
map_ax = fig.add_subplot(gs[2, :])
cax_red = fig.add_subplot(gs[1, 0])
cax_org = fig.add_subplot(gs[1, 1])
cax_blue = fig.add_subplot(gs[1, 2])

# add tree
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

contiguous_usa = gpd.read_file(gplt.datasets.get_path('contiguous_usa'))
contiguous_usa = contiguous_usa.to_crs(epsg=3857)

contiguous_usa.plot(linewidth=1,
    edgecolor='white',
    facecolor='lightgray',
    ax=map_ax);    


# drop ~10% of sampled trees (10k) as MCMC burn-in
drop_num = 1000
full_df = mrca_xy_crs_df[mrca_xy_crs_df["type"] == "west_full"][drop_num:] # west MRC location in full tree
subset_df = mrca_xy_crs_df[mrca_xy_crs_df["type"] == "west"][drop_num:]   # west tree root location
rescue_df = mrca_xy_crs_df[mrca_xy_crs_df["type"] == "mixed_alt"][drop_num:]  # west MRC location in WNV_west+52east.trees

# have to flip lat;lon -> lon;lat when plotting
sb.kdeplot(full_df.Longitude, 
           full_df.Latitude, 
           cmap=sb.light_palette("blue", as_cmap=True), 
           shade=True, 
           shade_lowest=False, 
           alpha=.7, 
           cbar=True, 
           cbar_ax=cax_blue,
           cbar_kws={"orientation": "horizontal"},
           zorder=5,
           ax=map_ax);

sb.kdeplot(subset_df.Longitude,
           subset_df.Latitude,
           cmap=sb.light_palette("red", as_cmap=True),
           shade=True, 
           shade_lowest=False,  
           alpha=.7, 
           cbar=True, 
           cbar_ax=cax_red,
           cbar_kws={"orientation": "horizontal"},
           ax=map_ax);

sb.kdeplot(rescue_df.Longitude,
           rescue_df.Latitude,
           cmap=sb.light_palette("orange", as_cmap=True),
           shade=True, 
           shade_lowest=False,  
           alpha=.7, 
           cbar=True, 
           cbar_ax=cax_org,
           cbar_kws={"orientation": "horizontal"},
           ax=map_ax);

cax_red.xaxis.set_ticks_position('top')
cax_red.xaxis.set_label_position('top')
cax_red.tick_params(labelsize=10, length=3, left=True, labelleft=True, 
                    right=False, labelright=False)
cax_red.set_xlabel("west MRCA location (KDE) in west analysis")


cax_org.xaxis.set_ticks_position('top')
cax_org.xaxis.set_label_position('top')
cax_org.tick_params(labelsize=10, length=3, left=True, labelleft=True, 
                    right=False, labelright=False)
cax_org.set_xlabel("west MRCA location (KDE) in extra samples analysis")

cax_blue.xaxis.set_ticks_position('top')
cax_blue.xaxis.set_label_position('top')
cax_blue.tick_params(labelsize=10, length=3, right=True, labelright=True)
cax_blue.set_xlabel("west MRCA location (KDE) in full analysis")

# show tip location w/t scatterplot
tip_df = tip_crs_df[(tip_crs_df["tip_type"] == "west") |\
                    (tip_crs_df["tip_type"] == "east")]
sb.scatterplot(x=tip_df.Longitude, 
               y=tip_df.Latitude,
               hue=tip_df.tip_type,
               palette = {"east" : "g",
                          "west" : "r"},
               style=tip_df.tip_type,
               markers ={"east" : style_dict["east"][3],
                         "west" : style_dict["west"][3]},
               s=20,
               ax = map_ax)

# annotate scatter with tip_idx label
X=tip_crs_df.Longitude.values
Y=tip_crs_df.Latitude.values
name_to_num = dict([(k.replace("'", ""), name_to_num[k]) for k in name_to_num.keys()])
labels = [name_to_num[l.replace(" ", "_")] for l in tip_crs_df["tip_id"]]
for i, label in enumerate(labels):
    map_ax.annotate(label, 
                    (X[i], 
                     Y[i]),
                   fontsize=6)


# custom legend
handles, labels = map_ax.get_legend_handles_labels()
map_ax.legend(handles[1:],  labels[1:], 
              loc='lower left', 
              title="Tip subset")
#[-14201398.96929843, -7133189.60868023, 2712783.1402848554, 6513955.476226102]
map_ax.set_xlim(left=-14201398.969, right=-7133189.608)
map_ax.set_ylim(bottom=2712783.140, top=6513955.476)
map_ax.set_axis_off()
gs.tight_layout(fig)
fig.show();
fig.savefig("figures/WNV_west1_vs_fullmrca_vs_west+52east_location_kde_tree.pdf", 
            format="pdf", dpi="figure", bbox_inches="tight")
fig.savefig("figures/WNV_west1_vs_fullmrca_vs_west+52east_location_kde_tree.svg", 
            format="svg", dpi="figure", bbox_inches="tight")

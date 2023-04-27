# generates the final circular dendrogram + heatmap that's present in the paper

## uncomment following to install pre-requisites to the libraries
# install.packages("BiocManager")
# BiocManager::install(version = "3.17")
# BiocManager::install("BiocUpgrade")
# BiocManager::install("ggtreeExtra")

library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(reshape)

setwd("~/code/plant-chemical-space/notebooks/")

tree <- read.newick("../data/taxonomy_tree_1B.nwk")

dat2 <-
  read.table(
    "../data/chemical_classes_family.tsv",
    sep = "\t",
    header = TRUE,
    row.names = 1,
    fill = TRUE
  )

column_order <-
  c(
    "Ursane.and.Taraxastane.triterpenoids",
    'Oleanane.triterpenoids',
    'Lupane.triterpenoids',
    "Menthane.monoterpenoids",
    'Acyclic.monoterpenoids',
    'Monocyclic.monoterpenoids',
    "Eudesmane.sesquiterpenoids",
    'Flavones',
    'Flavanones',
    'Flavonols',
    'Flavan.3.ols',
    'Dihydroflavonols',
    'Shikimic.acids.and.derivatives',
    'Cinnamic.acids.and.derivatives',
    "Stigmastane.steroids",
    "Gallotannins",
    "Chalcones",
    "Furofuranoid.lignans",
    "Neolignans",
    "Simple.coumarins"
  )

colors = c(
  '#fd8d3c',
  '#fdae6b',
  '#fdd0a2',

  '#3182bd',
  '#6baed6',
  '#9ecae1',

  "#e6550d",

  '#2ca02c',
  '#98df8a',

  '#d62728',
  '#ff9896',

  '#756bb1',

  '#e377c2',
  '#f7b6d2',
  '#7f7f7f',
  '#c7c7c7',
  '#bcbd22',
  '#dbdb8d',
  '#8c564b',
  '#c49c94'
)

dat2 <- dat2[, column_order]
# dat2 <- data.frame(lapply(dat2, function(x)
#   scale(
#     x,
#     center = min(x, na.rm = FALSE),
#     scale = diff(range(x, na.rm = FALSE))
#   )))

dat2$rowname <- rownames(dat2)
dat2 <- melt(dat2) # melting helps with plotting in geom_fruit()

# rescale for plotting
dat2$value <- dat2$value - min(dat2$value) / diff(range(dat2$value))


nodeids <- nodeid(tree, tree$node.label)
nodedf <- data.frame(node = nodeids)
labdf <-
  data.frame(node = nodeids, label = tree$node.label)

colors_from_tab20c = c(
  '#3182bd',
  '#6baed6',
  '#9ecae1',
  '#c6dbef',
  '#e6550d',
  '#fd8d3c',
  '#fdae6b',
  '#fdd0a2',
  '#31a354',
  '#74c476',
  '#a1d99b',
  '#c7e9c0',
  '#756bb1',
  '#9e9ac8',
  '#bcbddc',
  '#dadaeb',
  '#636363',
  '#969696',
  '#bdbdbd',
  '#d9d9d9'
)


colors_from_tab20 = c(
  '#1f77b4',
  '#aec7e8',
  '#ff7f0e',
  '#ffbb78',
  '#2ca02c',
  '#98df8a',
  '#d62728',
  '#ff9896',
  '#9467bd',
  '#c5b0d5',
  '#8c564b',
  '#c49c94',
  '#e377c2',
  '#f7b6d2',
  '#7f7f7f',
  '#c7c7c7',
  '#bcbd22',
  '#dbdb8d',
  '#17becf',
  '#9edae5'
)

# The circular layout tree.
p <-
  ggtree(
    tree,
    layout = "fan",
    # comment out for rectangular graph
    branch.length = 'none',
    size = 0.1,
    # line thickness
    open.angle = 90
  ) + # gives it a nice gap
  geom_tiplab(# align = TRUE, # only use if you need dotted / dashed lines from label to heatmap column
    geom = 'text',
    size = 0.5,
    # text size of labels
    offset = 0.04) + # of the labels from the tree)
  new_scale_fill() +
  geom_fruit(
    # lets you map whatever plot to a histogram
    data = dat2,
    geom = geom_tile,
    # basically a heat map
    mapping = aes(
      y = rowname,
      x = variable,
      fill = variable,
      alpha = value
    ),
    # if you set fill = "variable", then each column gets a different color
    color = "#666666",
    # of the bounding boxes in the grid
    offset = 0.05,
    # of the heat map from the
    size = 0.02,
    # width of lines
    # lwd = 0.1, # line width
    axis.params = list(
      # add labels to each column in the circular heat map
      axis = "x",
      # title = "disease area", # use illustrator for this instead
      text.angle = 45,
      # angle of text
      text.size = 1,
      # size of text
      line.size = 0,
      vjust = 0,
      hjust = 1
    )
  ) + scale_fill_manual(values = colors,
                        guide = guide_legend(
                          keywidth = 0.3,
                          keyheight = 0.3,
                          order = 4
                        ))
p

library(tidyr)


p <- ggplot(dat2, aes(x = value)) + geom_histogram()
p



# colors from seaborn's "deep" color palette (python)
pal = c(
  '#4c72b0',
  '#dd8452',
  '#55a868',
  '#c44e52',
  '#8172b3',
  '#937860',
  '#da8bc3',
  '#8c8c8c',
  '#ccb974',
  '#64b5cd'
)

### Other scale functions to try out.

## continuous color palette
# scale_color_continuous(high = "#4c72b0", low = "#FFFFFF") +

## specific viridis color palette that pre-bins your data before coloring
# scale_fill_viridis_b(option = "viridis", name = "value") #+ #bins before coloring

## can directly set a theme instead for lots of different plot elements
# theme_minimal() #+

## use this with the fill="variable" setting in the mapping for the heat map. This will allow for a gradient of
## transparency for each column.
# scale_alpha_continuous(range = c(0.0, 20),
#                        guide = 'none',guide_legend(
#                          keywidth = 0.1,
#                          keyheight = 0.3,
#                          order = 5
# ))

## i think this adds a scale for the tree
# geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1)

p
# add this option if you don't want the scales plotted
# + guides(fill = 'none', alpha = 'none')

## Save files as SVG using this package
#install.packages("svglite")
library(svglite)
ggsave(
  "fig_1_b.svg",
  #infinite zoom; good for editing in illustrator
  plot = p,
  width = 12,
  height = 12,
  units = 'in',
  dpi = 300 # good for publishing
)

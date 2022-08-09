# Visualising the Pangenome in Cytoscape

### Loading

Panaroo produces a GML file (`final_graph.gml`) for input to [Cytoscape](https://cytoscape.org/). 

This graph can be loaded by clicking on the 'Import Network' button at the top right of the Cytoscape application and navigating to the `final_graph.gml` file. We have indicated this in the figure below.

<img src="_figures/import_net.png" width="500">

### Layout

We recommend the [yFiles Organic Layout](http://apps.cytoscape.org/apps/yfileslayoutalgorithms) algorithm to arrange the nodes of the graph for visualisation. Once the [yFiles](http://apps.cytoscape.org/apps/yfileslayoutalgorithms) plugin has been installed this layout can be called by selecting it from the dropdown menu in Cytoscape. This is indicated below.

<img src="_figures/organic_layout.png" width="500">


### Annotation

The style of the nodes and edges in the graph can be used to convey extra information. These attributes can be controlled using the 'style' tab in cytoscape. Below we indicate how to colour the nodes based on the number of genomes that they contain. This helps to visualise the core and accessory genome based upon the colour of each node.

<img src="_figures/style.png" width="500">

The style tabs for edges, nodes and the network can be toggled at the bottom of the style tab as indicated below.

<img src="_figures/style_tab.png" width="500">

### Filtering

Sometimes it is useful to focus on a subset of the graph. Nodes can be selected using the mouse or rules in the 'select' tab. These can then either be converted into a subgraph for further analysis or filtered out from the main graph.

The figure below indicates how to select nodes using a set of rules. Here we are selecting all nodes that appear in between 1 and 3 genomes.

<img src="_figures/select.png" width="1000">

The buttons in the figure below can then be used to either create a new graph from the selected nodes or to filter them out.

<img src="_figures/filter.png" width="500">


### Additional Metadata

Additional metadata can be imported to cytoscape using the "Import table from file" button, next to the "Import network" button. These data can be associated with a graph and then user for additional filtering or annotation. The `name` attribute in the pangenome graph is a good option to index the imported spreadsheet for associating additional gene data with the pangenome graph. See [Meta-data Annotation](vis/metadata.md) for additional details. 

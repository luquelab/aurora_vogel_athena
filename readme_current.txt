The main bulk of the code is in "potential_unfolding_creator.py".

The "usual numbering" of vertices of icosahedral graph is
   0    1    2    3     4
 5   6     7   8     9    10
11 12   13  14   15    16
  17  18  19   20   21

Below line 388 are examples using the classes and functions imported/written from above.

The following is the basic structure of an example:
- filt1 (this tells PolyhedralGraph which vertices are glued together based on
         the numbering given above)
  pgraph (gives the connectivity of the polyhedral graph)
  top_vector_sequence (this is the top triangle path)
  bottom_vector_sequence (this is the bottom triangle path)
  delta (this is the separation vector)
  upper_path_data (inputs vector sequence and calculates relevant data like in thesis)
  lower_path_data (inputs vector sequence and calculates relevant data like in thesis)
  lumper (tells get_stats which vertices to treat as identical since some of the vertices
          on the graph are the same for the 3d structure, for example with the "usual
          numbering" vertex 0, 1, 2, 3, and 4 of the icosahedral graph are identical)
  cluster_analysis (does cluster analysis, !!! Comment out brute_modularity_optim in this function
                    unless you want the brute algorithm to run, which may take 5-10 minutes !!!)
  get_stats (gives some basic stats of the unfolding using lumper)
  graph_unfolding (graphs the unfolding on hexagonal lattice)

--------------------
"triangle_paths.py"
--------------------
Class Triangle
- creates triangle, used in TrianglePathData class, with nice properties (like being able to get all the sides
  angles, or area))
Class TrianglePathData
- creates triangle paths with their data (like in the thesis)

--------------------
--------------------
"polyhedral_graph.py"
--------------------
Class PolyhedralGraph
--------------------
- builds the graph of the unfolding using vertices and given connectivity
- records faces based edges and connectivity given from filter (filt1 in example)
--------------------
Class edge
--------------------
- sometimes used to represent edges (like in the PolyhedralGraph class)
- allows one to check whether a vertex is in an edge and whether two edges have a common vertex
--------------------
--------------------
"potential_unfolding_creator.py"
--------------------
Class PotentialSpacedUnfolding
--------------------
- creates an unfolding from triangle path data, separation vector, and polyhedral graph
--------------------
get_easy_stats(array)
- returns max, min, mean, median, std, mode from numpy array
get_stats(array, kind, lumper)
- array is numpy array, kind indications how the "lumping" is being done, where
  lumper indicates which vertices should be treated as the same vertex
- uses

brute_modularity_optim(G)
- calculates the modularity of our graph G using brute force (since our graphs are small
  we don't needs algorithms for estimating modularity, we can just calculate it)
- uses the functions get_all_possible_partitions, _get_all_possible_partitions_rec,
  _get_all_possible_subsets_rec (code cited from stacke exchange, see comment on lines 247-249)

graph_unfolding(unfolding, export, path, name)
- creates a a sufficiently large lattice (lat.make_tiling()) using information from the unfolding and
  graphs it
- iterates through faces from unfolding (unfolding.polyhedral_graph.faces) and their colors and plots them
- exports to path with file name if requested

--------------------
--------------------
"my_library.py"
--------------------
Class HexagonalPoint
--------------------
Implements addition (+), subtraction (-), multiplication (*) by a scalar.

reduce():
- makes 3rd coord 0, i.e. (h,k,l)=(h-l,k+l,0)

rotate_n_pi_3(n):
- rotates by n multiple of pi/3

get_cartesian():
- returns conversion of hexagonal point (h,k,l) to cartesian.

--------------------
Class HexagonalPoint
--------------------
Implements addition (+), subtraction (-), multiplication (*) by a scalar.

get_hexagonal():
- returns conversion of cartesian point [a,b] to hexagonal point.

--------------------
"polyhedral_graph.py"
--------------------
Used to create the excel sheets you have seen as well as the graphs. Currently some of this looks to be bugged as
the graphs don't look quite right. This code is lengthy and repetitive.

add_model_data()
- given the necessary inputs, it creates a PotentialSpacedUnfolding and then formats the relevant excel sheet
  extracts relevant stats, and creates graphs

lines 251 onwards goes through examples using the add_model_data function.
At the very end the excel sheet is closed and final figures are created / saved

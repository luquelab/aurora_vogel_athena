import xlsxwriter as xlsx
from potential_unfolding_creator import *
import os

###################################################################################
###################################################################################
edge_lengths_master = []
node_degrees_master = []
labels_master = []
workbook = xlsx.Workbook('../hiv_model/hiv_models.xlsx')
bold = workbook.add_format({'bold': True})
colors = ['red', 'blue', 'lime', 'tan', 'orange', 'purple', 'yellow', 'brown', 'black', 'cyan', 'lightcoral', 'peru']
###################################################################################
main_col_names = ["PDB_ID", "Structure Name", "Concise Description from Study", "Reference", "DOI",
                  "Visual Interpretation", "Proteins Calc", "Hexamers Calc", "Pentamers Calc", "Unfolding Match Experiment",
                  "Vertices", "Faces", "Edges", "Number of degree 4 vertices", "Number of degree 5 vertices",
                  "Number of degree 6 vertices", "Number of degree 7 vertices", "Max edge length", "Min edge length", "Mean edge length",
                  "Median edge length", "Std edge length", "Number of clusters", "Cluster sizes"]
# 16 to start proper
mast_row = 0
mast_col = 0
master_sheet = workbook.add_worksheet('all models')
# make names of master sheet
for i in range(len(main_col_names)):
    master_sheet.write(0, i, main_col_names[i], bold)
mast_row = 1

master_sheet.write(mast_row, 0, "3J3Y")
master_sheet.write(mast_row, 1, "supp_fig_4d")
master_sheet.write(mast_row, 2, "fullerene cone, atomic-level structure of entire HIV-1 capsid")
master_sheet.write(mast_row, 3, "Zhao et al 2013")
master_sheet.write(mast_row, 4, "https://doi.org/10.1038/nature12162")
master_sheet.write(mast_row, 5, "conical shape with 7 and 5 pentamers")


def add_model_data(monicer, current_sheet, pgraph, upper_path_data, lower_path_data, lumper, pdb='NA',
                   structure_name='', concise_descr='', reference='', doi='', visual_interp=''):
    global workbook, master_sheet, mast_row, mast_col, edge_lengths_master, node_degrees_master, labels_master
    # Make titles
    ###################################################################################
    # PDB 3J3Y #
    ###################################################################################
    STRUCT_NAME = monicer
    sheet_row = 0
    sheet_col = 0
    #############################
    unfolding = PotentialSpacedUnfolding(upper_path_data=upper_path_data, lower_path_data=lower_path_data, delta=delta,
                                         polyhedral_graph=pgraph)
    # community_data = cluster_analysis(unfolding, lumper, data=True)
    community_data = cluster_analysis(unfolding, lumper, data=True, calc_type='inv')
    stats_lengths = get_stats(unfolding.get_edge_lengths())
    stats_degrees = get_stats(unfolding.polyhedral_graph.get_degrees(), kind="degree list", lumper=lumper)
    try:
        os.makedirs("../hiv_model/")
    except FileExistsError:
        pass
    try:
        os.makedirs("../hiv_model/unfoldings/")
    except FileExistsError:
        pass
    try:
        os.makedirs("../hiv_model/histograms/")
    except FileExistsError:
        pass
    try:
        os.makedirs("../hiv_model/histograms/"+STRUCT_NAME)
    except FileExistsError:
        pass
    graph_unfolding(unfolding=unfolding, export=True, path="../hiv_model/", name=monicer+"_unfolding")
    #############################

    master_sheet.write(mast_row, 0, pdb)
    master_sheet.write(mast_row, 1, structure_name)
    master_sheet.write(mast_row, 2, concise_descr)
    master_sheet.write(mast_row, 3, reference)
    master_sheet.write(mast_row, 4, doi)
    master_sheet.write(mast_row, 5, visual_interp)
    master_sheet.write(mast_row, 6, unfolding.get_total_subunits())
    master_sheet.write(mast_row, 7, unfolding.get_total_subunits()-12*5)
    master_sheet.write(mast_row, 8, 12*5)
    master_sheet.write(mast_row, 9, "Yes")
    master_sheet.write(mast_row, 10, 12)
    master_sheet.write(mast_row, 11, 20)
    master_sheet.write(mast_row, 12, 30)
    # degree 4
    master_sheet.write(mast_row, 13, 0)
    for dat in stats_degrees["Frequencies"]:
        if dat[0] == 4:
            master_sheet.write(mast_row, 13, dat[1])
            written = True
    # degree 5
    master_sheet.write(mast_row, 14, 0)
    for dat in stats_degrees["Frequencies"]:
        if dat[0] == 5:
            master_sheet.write(mast_row, 14, dat[1])
            written = True
    # degree 6
    master_sheet.write(mast_row, 15, 0)
    for dat in stats_degrees["Frequencies"]:
        if dat[0] == 6:
            master_sheet.write(mast_row, 15, dat[1])
            written = True
    # degree 7
    master_sheet.write(mast_row, 16, 0)
    for dat in stats_degrees["Frequencies"]:
        if dat[0] == 7:
            master_sheet.write(mast_row, 16, dat[1])
            written = True

    master_sheet.write(mast_row, 17, stats_lengths["Max"])
    master_sheet.write(mast_row, 18, stats_lengths["Min"])
    master_sheet.write(mast_row, 19, stats_lengths["Mean"])
    master_sheet.write(mast_row, 20, stats_lengths["Median"])
    master_sheet.write(mast_row, 21, stats_lengths["Std"])
    master_sheet.write(mast_row, 22, len(community_data["brute"]["community"]))
    comm_sizes = []
    for comm in community_data["brute"]["community"]:
        comm_sizes.append(len(comm))
    master_sheet.write(mast_row, 23, str(comm_sizes)[1:-1])
    #############################
    sheet_col_names = ["Community edge max", "Community edge min", "Community edge mean", "Community edge median", "Community edge std",
                       "Community degree max", "Community degree min", "Community degree mean", "Community degree median", "Community degree std",
                       "Type", "Community"]
    for i in range(len(sheet_col_names)):
        current_sheet.write(sheet_row, i, sheet_col_names[i], bold)
    sheet_row = 1
    #############################
    G = community_data["graph"]
    edge_collection_brute = []
    node_collection_brute = []
    edge_collection_greedy = []
    node_collection_greedy = []
    comm_counter = 0
    # adding to cumulative edges
    edge_lengths_to_add = []
    node_degrees_to_add = []
    for edge in list(G.edges()):
        edge_lengths_to_add.append(G[edge[0]][edge[1]]['length'])
    for node in list(G.nodes()):
        node_degrees_to_add.append(G.degree[node])
    edge_lengths_master.append(edge_lengths_to_add)
    node_degrees_master.append(node_degrees_to_add)
    labels_master.append(monicer)
    # for communities
    for comm in community_data["brute"]["community"]:
        comm_counter += 1
        edge_lengths = []
        node_degrees = []
        G_sub = G.subgraph(comm)
        for edge in list(G_sub.edges()):
            edge_lengths.append(G[edge[0]][edge[1]]['length'])
        for node in list(G_sub.nodes()):
            node_degrees.append(G.degree[node])
        edge_stats = get_easy_stats(np.array(edge_lengths))
        edge_collection_brute.append(np.array(edge_lengths))
        node_stats = get_easy_stats(np.array(node_degrees))
        node_collection_brute.append(np.array(node_degrees))
        #
        current_sheet.write(sheet_row, 0, edge_stats["Max"])
        current_sheet.write(sheet_row, 1, edge_stats["Min"])
        current_sheet.write(sheet_row, 2, edge_stats["Mean"])
        current_sheet.write(sheet_row, 3, edge_stats["Median"])
        current_sheet.write(sheet_row, 4, edge_stats["Std"])
        current_sheet.write(sheet_row, 5, node_stats["Max"])
        current_sheet.write(sheet_row, 6, node_stats["Min"])
        current_sheet.write(sheet_row, 7, node_stats["Mean"])
        current_sheet.write(sheet_row, 8, node_stats["Median"])
        current_sheet.write(sheet_row, 9, node_stats["Std"])
        current_sheet.write(sheet_row, 10, "brute")
        current_sheet.write(sheet_row, 11, str(comm)[1:-1])
        #
        plt.hist(np.array(edge_lengths))
        plt.title("BC"+str(comm_counter)+" Edge Length Distribution")
        plt.xlabel("Edge Length")
        plt.ylabel("Number of Edges")
        plt.savefig('../hiv_model/histograms/'+STRUCT_NAME+'/'+STRUCT_NAME+'_edge_length_brute_comm_'+str(comm_counter)+'_histo')  # save the figure to file
        plt.close()
        current_sheet.insert_image(sheet_row, 12, '../hiv_model/histograms/'+STRUCT_NAME+'/'+STRUCT_NAME+'_edge_length_brute_comm_' + str(comm_counter)
                                   + '_histo.png')
        #
        plt.hist(np.array(node_degrees))
        plt.title("BC"+str(comm_counter)+" Degree Distribution")
        plt.xlabel("Degree")
        plt.ylabel("Number of Nodes")
        plt.savefig('../hiv_model/histograms/'+STRUCT_NAME+'/'+STRUCT_NAME+'_node_degree_brute_comm_'+str(comm_counter)+'_histo')  # save the figure to file
        plt.close()
        current_sheet.insert_image(sheet_row, 13, '../hiv_model/histograms/'+STRUCT_NAME+'/'+STRUCT_NAME+'_node_degree_brute_comm_'+str(comm_counter) + '_histo.png')
        sheet_row += 1

    comm_counter = 0
    for comm in community_data["greedy"]:
        comm_counter += 1
        edge_lengths = []
        node_degrees = []
        G_sub = G.subgraph(comm)
        for edge in list(G_sub.edges()):
            edge_lengths.append(G[edge[0]][edge[1]]['length'])
        for node in list(G_sub.nodes()):
            node_degrees.append(G.degree[node])
        edge_stats = get_easy_stats(np.array(edge_lengths))
        edge_collection_greedy.append(np.array(edge_lengths))
        node_stats = get_easy_stats(np.array(node_degrees))
        node_collection_greedy.append(np.array(node_degrees))
        #
        current_sheet.write(sheet_row, 0, edge_stats["Max"])
        current_sheet.write(sheet_row, 1, edge_stats["Min"])
        current_sheet.write(sheet_row, 2, edge_stats["Mean"])
        current_sheet.write(sheet_row, 3, edge_stats["Median"])
        current_sheet.write(sheet_row, 4, edge_stats["Std"])
        current_sheet.write(sheet_row, 5, node_stats["Max"])
        current_sheet.write(sheet_row, 6, node_stats["Min"])
        current_sheet.write(sheet_row, 7, node_stats["Mean"])
        current_sheet.write(sheet_row, 8, node_stats["Median"])
        current_sheet.write(sheet_row, 9, node_stats["Std"])
        current_sheet.write(sheet_row, 10, "greedy")
        current_sheet.write(sheet_row, 11, str(list(comm))[1:-1])
        #
        plt.hist(np.array(edge_lengths))
        plt.title("GC" + str(comm_counter)+" Edge Length Distribution")
        plt.xlabel("Edge Length")
        plt.ylabel("Number of Edges")
        plt.savefig('../hiv_model/histograms/'+STRUCT_NAME+'/'+STRUCT_NAME+'_edge_length_greedy_comm_' + str(comm_counter) + '_histo')  # save the figure to file
        plt.close()
        #
        plt.hist(np.array(node_degrees))
        plt.title("GC" + str(comm_counter) + " Node Degree Distribution")
        plt.xlabel("Degree")
        plt.ylabel("Number of Nodes")
        plt.savefig('../hiv_model/histograms/'+STRUCT_NAME+'/'+STRUCT_NAME+'_node_degree_greedy_comm_' + str(comm_counter) + '_histo')  # save the figure to file
        plt.close()
        sheet_row += 1

    current_sheet.autofit()
    #############################
    plt.hist(edge_collection_brute, color=colors[0:len(edge_collection_brute)])
    plt.title("BC Edge Length Distribution")
    plt.xlabel("Edge Length")
    plt.ylabel("Number of Edges")
    plt.savefig('../hiv_model/histograms/'+STRUCT_NAME+'_edge_length_mast_histo')  # save the figure to file
    plt.close()

    plt.hist(node_collection_brute, color=colors[0:len(node_collection_brute)])
    plt.title("BC Node Degree Distribution")
    plt.xlabel("Node Degree")
    plt.ylabel("Number of Nodes")
    plt.savefig('../hiv_model/histograms/'+STRUCT_NAME+'_node_degree_mast_histo')  # save the figure to file
    plt.close()

    mast_row += 1

###################################################################################
################################################################################### 1, PDB_3J3Y
###################################################################################
filt1 = [[0, 5], [1, 6], [2, 7], [3, 8], [4, 9],
             [5, 11],
             [11, 17], [12, 18], [13, 19], [14, 20], [15, 21]]
pgraph=pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                       [0, 11], [1, 13], [7, 8], [8, 9], [9, 10], [9, 15], [10, 15],
                                       [5, 11], [6, 11], [6, 12], [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [10, 16],
                                       [11, 12], [12, 13], [13, 14], [14, 15], [10, 21],
                                       [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                                filt=filt1)
# paths and delta
top_vector_sequence = [np.array([4, 1, 0]),
          np.array([1, 7, 0]),
          np.array([0, 6, 0]),
          np.array([1, 3, 0]),
          np.array([1, 4, 0]),
          np.array([4, 1, 0])]
bottom_vector_sequence = [np.array([1, 2, 0]),
          np.array([1, 2, 0]),
          np.array([1, 3, 0]),
          np.array([5, 4, 0]),
          np.array([0, 9, 1]),
          np.array([1, 2, 0])]
delta = lib.HexagonalPoint(h=-9, l=-1)
upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]

add_model_data(monicer="PDB_3J3Y",
               current_sheet=workbook.add_worksheet('PDB_3J3Y'),
               pgraph=pgraph,
               upper_path_data=upper_path_data,
               lower_path_data=lower_path_data,
               lumper=lumper,
               pdb="3J3Y",
               structure_name="supp_fig_4d",
               concise_descr="fullerene cone, atomic-level structure of entire HIV-1 capsid",
               reference="Zhao et al 2013",
               doi="https://doi.org/10.1038/nature12162",
               visual_interp="conical shape with 7 and 5 pentamers"
               )

###################################################################################
################################################################################### 2, PDB-3J3Q
###################################################################################
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [7, 8], [8, 9], [4, 15], [9, 15], [10, 15],
                                   [5, 11], [5, 12], [6, 12], [7, 12], [7, 13], [8, 13], [9, 13], [9, 14], [10, 16],
                                   [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)
# paths and delta
top_vector_sequence = [np.array([0, 5, 0]),
          np.array([0, 5, 0]),
          np.array([3, 2, 0]),
          np.array([0, 2, 0]),
          np.array([3, 7, 0]),
          np.array([0, 5, 0])]
bottom_vector_sequence = [np.array([3, 5, 0]),
          np.array([0, 9, 0]),
          np.array([1, 2, 0]),
          np.array([1, 2, 0]),
          np.array([1, 3, 0]),
          np.array([3, 5, 0])]
delta = lib.HexagonalPoint(h=-4, k=2)

upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# lumper
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
add_model_data(monicer="PDB_3J3Q",
               current_sheet=workbook.add_worksheet('PDB_3J3Q'),
               pgraph=pgraph,
               upper_path_data=upper_path_data,
               lower_path_data=lower_path_data,
               lumper=lumper,
               pdb="3J3Q",
               structure_name="supp_fig_7b",
               concise_descr="fullerene cone, atomic-level structure of entire HIV-1 capsid",
               reference="Zhao et al 2013",
               doi="https://doi.org/10.1038/nature12162",
               visual_interp="conical shape with 7 and 5 pentamers"
               )
###################################################################################
################################################################################### 3, mattei_1_1
###################################################################################
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [7, 8], [3, 14], [9, 10],
                                   [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 14], [9, 15], [9, 16], [10, 16],
                                   [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)
# paths and delta
top_vector_sequence = [np.array([0, 10, 2]),
          np.array([9, 4, 0]),
          np.array([5, 0, -8]),
          np.array([2, 0, 0]),
          np.array([0, 5, 0]),
          np.array([0, 10, 2])]
bottom_vector_sequence = [np.array([11, 0, 0]),
          np.array([2, 10, 0]),
          np.array([2, 1, 0]),
          np.array([4, 1, 0]),
          np.array([3, 1, 0]),
          np.array([11, 0, 0])]
delta = lib.HexagonalPoint(l=2.0)

upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# lumper
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
add_model_data(monicer="mattei_1_1",
               current_sheet=workbook.add_worksheet('mattei_1_1'),
               pgraph=pgraph,
               upper_path_data=upper_path_data,
               lower_path_data=lower_path_data,
               lumper=lumper,
               pdb="NA",
               structure_name="mov_6_grid_1_1",
               concise_descr="complete fullerene with twelve pentamers",
               reference="Mattei et al 2016",
               doi="https://doi.org/10.1126/science.aah4972",
               visual_interp="conical shape with 7 and 5 pentamers"
               )

###################################################################################
################################################################################### 4, mattei_1_2
###################################################################################
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [2, 14], [8, 9], [9, 10],
                                   [5, 11], [5, 12], [6, 12], [7, 12], [7, 13], [7, 14], [8, 14], [8, 15], [9, 15], [10, 15], [10, 16],
                                   [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)
# paths and delta
top_vector_sequence = [np.array([0, 2, 9]),
          np.array([0, 9, 2]),
          np.array([6, 5, 0]),
          np.array([0, 1, 2]),
          np.array([0, 1, 1]),
          np.array([0, 2, 9])]
bottom_vector_sequence = [np.array([0, 3, 1]),
          np.array([0, 1, 3]),
          np.array([0, 2, 0]),
          np.array([3, 9, 0]),
          np.array([0, 6, 7]),
          np.array([0, 3, 1])]
delta = lib.HexagonalPoint(h=-3.0, l=1.0)

upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# lumper
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
add_model_data(monicer="mattei_1_2",
               current_sheet=workbook.add_worksheet('mattei_1_2'),
               pgraph=pgraph,
               upper_path_data=upper_path_data,
               lower_path_data=lower_path_data,
               lumper=lumper,
               pdb="NA",
               structure_name="mov_6_grid_1_2",
               concise_descr="complete fullerene with twelve pentamers",
               reference="Mattei et al 2016",
               doi="https://doi.org/10.1126/science.aah4972",
               visual_interp="conical shape with 7 and 5 pentamers"
               )

###################################################################################
################################################################################### 5, mattei_1_3
###################################################################################
filt2 = [[0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [5, 11],
         [6, 13],
         [13, 18], [14, 19], [15, 20], [16, 21]]
pgraph2 = pg.PolyhedralGraph(edges=[[0, 6], [0, 7], [1, 7], [1, 8], [2, 8], [2, 9], [3, 9], [3, 10], [4, 10], [4, 11], [5, 11], [5, 12],
                                   [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12],
                                   [6, 13], [7, 13], [7, 14], [8, 14], [8, 15], [9, 15], [10, 15], [10, 16], [11, 16], [11, 17], [12, 17],
                                   [13, 14], [14, 15], [15, 16], [16, 17],
                                   [13, 18], [14, 18], [14, 19], [15, 19], [15, 20], [16, 20], [16, 21], [17, 21]],
                            filt=filt2)
# paths and delta
top_vector_sequence = [np.array([1, 2, 0]),
          np.array([0, 4, 0]),
          np.array([1, 4, 0]),
          np.array([3, 2, 0]),
          np.array([4, 0, 0]),
          np.array([2, 1, 0]),
          np.array([2, 0, -1])]
bottom_vector_sequence = [np.array([0, 1, 0]),
          np.array([1, 1, 0]),
          np.array([1, 1, 0]),
          np.array([1, 1, 0]),
          np.array([1, 0, 0])]
delta = lib.HexagonalPoint(h=-7.0, l=6.0)

upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# lumper
lumper = [[0, 1, 2, 3, 4, 5], [6, 6 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4], [-5, -5 - len(lower_path_data.face_triangles)]]
add_model_data(monicer="mattei_1_3",
               current_sheet=workbook.add_worksheet('mattei_1_3'),
               pgraph=pgraph2,
               upper_path_data=upper_path_data,
               lower_path_data=lower_path_data,
               lumper=lumper,
               pdb="NA",
               structure_name="mov_6_grid_1_3",
               concise_descr="complete fullerene with twelve pentamers",
               reference="Mattei et al 2016",
               doi="https://doi.org/10.1126/science.aah4972",
               visual_interp="conical shape with 7 and 5 pentamers"
               )

###################################################################################
################################################################################### 6, mattei_2_1
###################################################################################
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
                                   [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 14], [9, 15], [9, 16], [10, 16],
                                   [6, 17], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)
# paths and delta
top_vector_sequence = [np.array([0, 1, 1]),
          np.array([0, 2, 10]),
          np.array([0, 11, 1]),
          np.array([0, 1, 1]),
          np.array([0, 1, 1]),
          np.array([0, 1, 1])]
bottom_vector_sequence = [np.array([-2, 0, 7]),
          np.array([0, 2, 0]),
          np.array([0, 2, 1]),
          np.array([1, 8, 0]),
          np.array([0, 5, 5]),
          np.array([-2, 0, 7])]
delta = lib.HexagonalPoint(h=-5.0, k=-2.0)

upper_path_data = tpf.TrianglePathData(vector_sequence=top_vector_sequence)
lower_path_data = tpf.TrianglePathData(vector_sequence=bottom_vector_sequence, sigma=-1)
# lumper
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]
add_model_data(monicer="mattei_2_1",
               current_sheet=workbook.add_worksheet('mattei_2_1'),
               pgraph=pgraph,
               upper_path_data=upper_path_data,
               lower_path_data=lower_path_data,
               lumper=lumper,
               pdb="NA",
               structure_name="mov_6_grid_2_1",
               concise_descr="complete fullerene with twelve pentamers",
               reference="Mattei et al 2016",
               doi="https://doi.org/10.1126/science.aah4972",
               visual_interp="boomerang shape with 5 and 4 and 3 pentamers"
               )

###################################################################################
################################################################################### 7, mattei_2_2
###################################################################################
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
                                   [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 14], [9, 15], [10, 15], [10, 16],
                                   [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)

upper_path_data = tpf.TrianglePathData(vector_sequence=[[2, 0, 0], [1, 1, 0], [1, 2, 0], [2, 2, 0], [2, 1, 0], [2, 0, 0]])
lower_path_data = tpf.TrianglePathData(vector_sequence=[[2, 2, 0], [1, 2, 0], [1, 1, 0], [1, 1, 0], [3, 0, 0], [2, 2, 0]], sigma=-1)
delta = lib.HexagonalPoint(h=-1.0, l=11.0)
# lumper
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]

add_model_data(monicer="mattei_2_2",
               current_sheet=workbook.add_worksheet('mattei_2_2'),
               pgraph=pgraph,
               upper_path_data=upper_path_data,
               lower_path_data=lower_path_data,
               lumper=lumper,
               pdb="NA",
               structure_name="mov_6_grid_2_2",
               concise_descr="complete fullerene with twelve pentamers",
               reference="Mattei et al 2016",
               doi="https://doi.org/10.1126/science.aah4972",
               visual_interp="spherocylindrical shape with 6 and 6 pentamers"
               )

###################################################################################
################################################################################### 8, mattei_2_3
###################################################################################
pgraph = pg.PolyhedralGraph(edges=[[0, 5], [0, 6], [1, 6], [1, 7], [2, 7], [2, 8], [3, 8], [3, 9], [4, 9], [4, 10],
                                   [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
                                   [5, 11], [6, 11], [6, 12], [7, 12], [7, 13], [8, 13], [8, 14], [9, 14], [9, 15], [10, 15], [10, 16],
                                   [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
                                   [11, 17], [12, 17], [12, 18], [13, 18], [13, 19], [14, 19], [14, 20], [15, 20], [15, 21], [16, 21]],
                            filt=filt1)

# paths and delta
upper_path_data = tpf.TrianglePathData(vector_sequence=[[2, 1, 0], [3, 2, 0], [7, 3, 0], [1, 1, 0], [2, 1, 0], [2, 1, 0]])
lower_path_data = tpf.TrianglePathData(vector_sequence=[[3, 3, 0], [4, 2, 0], [2, 3, 0], [2, 0, 0], [4, 0, 0], [3, 3, 0]], sigma=-1)
delta = lib.HexagonalPoint(h=-1.0, l=7.0)
# lumper
lumper = [[0, 1, 2, 3, 4], [5, 5 + len(upper_path_data.face_triangles)], [-1, -2, -3, -4, -5], [-6, -6 - len(upper_path_data.face_triangles)]]

add_model_data(monicer="mattei_2_3",
               current_sheet=workbook.add_worksheet('mattei_2_3'),
               pgraph=pgraph,
               upper_path_data=upper_path_data,
               lower_path_data=lower_path_data,
               lumper=lumper,
               pdb="NA",
               structure_name="mov_6_grid_2_3",
               concise_descr="complete fullerene with twelve pentamers",
               reference="Mattei et al 2016",
               doi="https://doi.org/10.1126/science.aah4972",
               visual_interp="potato shape with 3 and 3 and 3 and 3 pentamers"
               )

###################################################################################
master_sheet.autofit()
workbook.close()

# Final histograms
# handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in [low,medium, high]]
plt.hist(edge_lengths_master, color=colors[0:len(edge_lengths_master)])
plt.legend(labels_master)
plt.title("Edge Length Distribution")
plt.xlabel("Edge Length")
plt.ylabel("Number of Edges")
plt.savefig('../hiv_model/edge_length_global_histo')  # save the figure to file
plt.close()

plt.hist(node_degrees_master, color=colors[0:len(node_degrees_master)])
plt.legend(labels_master)
plt.title("Node Degree Distribution")
plt.xlabel("Node Degree")
plt.ylabel("Number of Nodes")
plt.savefig('../hiv_model/node_degree_global_histo')  # save the figure to file
plt.close()

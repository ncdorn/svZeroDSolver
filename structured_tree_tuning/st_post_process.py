import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from struct_tree_utils import *
import json
from math import trunc
import random

# BinaryTree class for structured tree visualization (thank you chatGPT)
class BinaryTree:
    def __init__(self, value):
        self.value = value
        self.left = None
        self.right = None

# Function to create a binary tree from a list of values
def create_binary_tree(nodes):
    if not nodes:
        return None

    root = BinaryTree(nodes[0])
    queue = [root]
    i = 1

    while i < len(nodes):
        current_node = queue.pop(0)

        if nodes[i] is not None:
            current_node.left = BinaryTree(nodes[i])
            queue.append(current_node.left)

        i += 1

        if i >= len(nodes):
            break

        if nodes[i] is not None:
            current_node.right = BinaryTree(nodes[i])
            queue.append(current_node.right)

        i += 1



    return root


# Function to visualize the binary tree
def visualize_binary_tree(root,
                          labels,
                          vessel_ds,
                          vessel_lengths,
                          figure_number=0):
    G = nx.Graph()
    edges = [] # need to figure out how to add variable edge length and labels
    def traverse(node, parent='outlet', idx=None):
        if node is None:
            return

        G.add_node(node.value)
        if parent is not None:
            # set up edges dict to add in edge lengths later
            edges.append((parent, node.value))
            G.add_edge(parent, node.value)

        traverse(node.left, node.value)
        traverse(node.right, node.value)

    traverse(root)
    edge_lengths = {edges[i]: vessel_lengths[i] for i in range(len(edges))}
    edge_labels = {edges[i]: labels['edges'][i] for i in range(len(edges))}
    nx.set_edge_attributes(G, edge_lengths, name='weight')


    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
    # shift the first two nodes to the center
    pos[0] = (pos[0][0] + 125, pos[0][1])
    # move the outlet right above the first node
    pos['outlet'] = (pos[0][0], pos[0][1] + 50)

    plt.figure(figure_number)

    # need to add in total resistance
    nx.draw_networkx(G,
                     pos,
                     with_labels=True,
                     labels=labels['nodes'],
                     node_color="red",
                     node_shape='s',
                     node_size=1,
                     font_size=7,
                     font_weight="bold",
                     font_color="k",
                     width=[vessel_d * 10 for vessel_d in vessel_ds],
                     edge_color='red'
                     )
    nx.draw_networkx_edge_labels(G,
                                 pos,
                                 edge_labels=edge_labels,
                                 font_size=6
                                 )

def build_tree_figures(tree_config, fig):
    vessel_ids = []
    label_dict = {'nodes': {'outlet': 'outlet D = ' + str(round(tree_config["origin_d"], 3)) + '\n' +
                                      'tree D = ' + str(round(tree_config["vessels"][0]["vessel_D"], 3))},
                  'edges': {}}
    diameters = []
    lengths = []
    for vessel in tree_config["vessels"]:
        vessel_ids.append(vessel["vessel_id"])
        label_dict['edges'][vessel["vessel_id"]] = 'R = ' + str(trunc(vessel["zero_d_element_values"].get("R_poiseulle")))
        diameters.append(vessel["vessel_D"])
        lengths.append(vessel["vessel_length"])

    # print(resistances)
    # print(generations)
    structured_tree = create_binary_tree(vessel_ids)
    visualize_binary_tree(structured_tree,
                          label_dict,
                          diameters,
                          lengths,
                          fig)
    plt.title(tree_config["name"])

def visualize_trees(input_file):
    # load input file
    with open(input_file) as ff:
        config = json.load(ff)

    fig = 0
    for vessel_config in config["vessels"]:
        if "boundary_conditions" in vessel_config:
            if "outlet" in vessel_config["boundary_conditions"]:

                build_tree_figures(vessel_config["tree"], fig)
                fig += 1

    plt.show()

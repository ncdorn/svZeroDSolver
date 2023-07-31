import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from struct_tree_utils import *
import json
from math import trunc
import seaborn as sns
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

def visualize_binary_tree(root,
                          labels,
                          # vessel_ds,
                          vessel_lengths,
                          last_vessel=None,
                          edge_labeling=False):
    G = nx.Graph()
    edges = []  # need to figure out how to add variable edge length and labels
    vessel_ds = []
    print(root.d)
    def traverse(node, parent='outlet'):
        if node is None:
            return

        G.add_node(node.id)
        if parent is not None:
            # set up edges dict to add in edge lengths later
            edges.append((parent, node.id))
            vessel_ds.append(node.d)
            G.add_edge(parent, node.id)

        traverse(node.left, node.id)
        traverse(node.right, node.id)

    traverse(root)
    edges = edges[:last_vessel]
    # edge_lengths = {edges[i]: vessel_lengths[i] for i in range(len(edges))}
    # edge_labels = {edges[i]: labels['edges'][i] for i in range(len(edges))}
    edge_labels = {edges[i]: round(vessel_ds[i], 3) for i in range(len(edges))}
    # nx.set_edge_attributes(G, edge_lengths, name='weight')

    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
    # shift the first two nodes to the center
    pos[0] = (pos[0][0] + 125, pos[0][1])
    # move the outlet right above the first node
    pos['outlet'] = (pos[0][0], pos[0][1] + 50)

    plt.figure()

    # need to add in total resistance
    nx.draw_networkx(G,
                     pos,
                     with_labels=True,
                     labels = labels["nodes"],
                     node_color="red",
                     node_shape='s',
                     node_size=0,
                     font_size=7,
                     font_weight="bold",
                     font_color="k",
                     width=[vessel_d * 10 for vessel_d in vessel_ds],
                     edge_color='red'
                     )

    if edge_labeling:
        nx.draw_networkx_edge_labels(G,
                                     pos,
                                     edge_labels=edge_labels,
                                     font_size=6
                                     )


def build_tree_figure(tree_config, root, last_vessel=None, edge_labeling=False, fig_dir=None, fig_name=None):
    vessel_ids = []
    label_dict = {'nodes': {'outlet': 'outlet D = ' + str(round(tree_config["origin_d"], 3)) + '\n' +
                                      'tree D = ' + str(round(tree_config["vessels"][0]["vessel_D"], 3))},
                  'edges': {}}
    diameters = []
    lengths = []

    for vessel in tree_config["vessels"][:last_vessel]:
        vessel_ids.append(vessel["vessel_id"])
        label_dict['edges'][vessel["vessel_id"]] = 'R = ' + str(trunc(vessel["zero_d_element_values"].get("R_poiseulle")))
        diameters.append(vessel["vessel_D"])
        lengths.append(vessel["vessel_length"])
    # print(resistances)
    # print(generations)
    # structured_tree = create_binary_tree(vessel_ids)
    structured_tree = root
    visualize_binary_tree(structured_tree, # visualize the tree
                          label_dict,
                          # diameters,
                          lengths,
                          last_vessel,
                          edge_labeling=edge_labeling)

    plt.title(tree_config["name"] + '_'+ fig_name) # title the tree figure

    if fig_dir is not None: # save the figure if a directory is specified
        plt.savefig(str(fig_dir) + '/' + tree_config["name"] + '_' + str(fig_name) + '_visualized.png')
    else:
        plt.show()

def visualize_trees(config, roots, fig=None, fig_dir=None, fig_name=None):
    # method to visualze all trees
    fig_num = fig
    for vessel_config in config["vessels"]:
        if "boundary_conditions" in vessel_config:
            if "outlet" in vessel_config["boundary_conditions"]:
                for root in roots:
                    if root.name in vessel_config["tree"]["name"]:
                        print('building tree vis for ' + root.name)
                        build_tree_figure(vessel_config["tree"], root, fig_dir=fig_dir, fig_name=fig_name)
                        # plot_vessels_per_generation()
                        if fig is not None:
                            fig_num += 1

def plot_vessels_per_generation(tree_config: dict, name=None):
    '''
    plot a bar chart for the number of vessels in each tree generation level.
    Args:
        tree_config: config dict of tree
        name: extra naming convention to add onto the tree["name"]

    Returns:
        A bar chart of number of vessels plotted against tree generation

    '''
    gen_list = []
    for vessel in tree_config["vessels"]:
        gen_list.append(vessel["generation"])

    gen_count = []
    for i in range(max(gen_list) + 1):
        gen_count.append(gen_list.count(i))

    # bar chart comparing number of vessels per generation to generation
    plt.figure()
    fig, ax = plt.subplots()

    bars = ax.bar(range(max(gen_list) + 1),
            gen_count,
            tick_label=range(max(gen_list) + 1),
            log=True)
    # create the bar labels
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}',  # Text label
                    xy=(bar.get_x() + bar.get_width() / 2, height),  # Position of the label
                    xytext=(0, 3),  # Offset (to move the text above the bar)
                    textcoords="offset points",
                    ha='center', va='bottom')  # Alignment of the label
    # add in an annotation for d_out and d_min
    ax.annotate('D_out = ' + str(round(tree_config["origin_d"], 3)) + 'cm \n' + 'D_min = ' + str(tree_config["D_min"]) + 'cm',
                xy = (0, 512))
    # axis labels and title
    ax.set_xlabel('generation number')
    ax.set_ylabel('number of vessels')
    ax.set_title('Number of vessels per tree generation for ' + f'{tree_config["name"]}')


def plot_terminal_vessel_diameter(tree_config: dict, fig: int=1):
    terminal_dias = []
    terminal_gens = []
    for vessel in tree_config["vessels"]:
        if vessel["vessel_D"] < tree_config["r_min"]:
            terminal_dias.append(vessel["vessel_D"])
            terminal_gens.append(vessel["generation"])

    # scatter plot of terminal diameters
    plt.figure(fig)
    sns.swarmplot(terminal_dias)
    # axis labels and title
    plt.ylabel('terminal diameter (mm)')
    plt.title('Diameter of ' + str(len(terminal_dias)) + ' terminal vessels in a structured tree with r_min ' + str(tree_config["r_min"]))

    # plot terminal vessel diameter vs generation



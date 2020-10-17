"""Perform assembly based on debruijn graph."""


#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

import os
import sys
import random
import argparse
import statistics

import networkx as nx

# Set the random seed
random.seed(9001)

__author__ = "Vander Meersche Yann"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Yann Vander Meersche"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Vander Meersche Yann"
__email__ = "yann-vm@hotmail.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Reads a Fastq file.

    Arguments
    ---------
    fastq_file: Path to the Fastq file

    Yield
    -----
    Generator yielding Fasta reads
    """
    with open(fastq_file, "r") as f_in:
        for _ in f_in:
            yield next(f_in).strip()
            next(f_in)
            next(f_in)


def cut_kmer(fasta_seq, k):
    """Reads a Fastq file.

    Arguments
    ---------
    fasta_seq: Fasta sequence

    Yield
    -----
    Generator yielding Fasta sequence k-mers
    """
    for i in range(len(fasta_seq) - k+1):
        yield fasta_seq[i:i+k]


def build_kmer_dict(fastq_file, k):
    """Reads a Fastq file.

    Arguments
    ---------
    fasta_seq: Fasta sequence
    k: k-mer size

    Returns
    -------
    kmer_dict: k-mer dictionary
    """
    kmer_dict = {}
    for fasta_seq in read_fastq(fastq_file):
        for kmer in cut_kmer(fasta_seq, k):
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 0
            kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    """Build a NetworkX graph from a k-mer dictionary.

    Arguments
    ---------
    kmer_dict: k-mer dictionary

    Returns
    -------
    graph: NetworkX graph
    """
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        graph.add_edge(kmer[:-1], kmer[1:], weight=kmer_dict[kmer])
    return graph


def get_starting_nodes(graph):
    """Get the starting nodes from the graph.

    Arguments
    ---------
    graph: NetworkX graph

    Returns
    -------
    start_nodes: Starting nodes list
    """
    start_nodes = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            start_nodes.append(node)
    return start_nodes


def get_sink_nodes(graph):
    """Get the sinking nodes from the graph.

    Arguments
    ---------
    graph: NetworkX graph

    Returns
    -------
    sink_nodes: Sinking nodes list
    """
    sink_nodes = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(graph, start_nodes, sink_nodes):
    """Get all the contigs between the staring nodes and the sinking nodes.

    Arguments
    ---------
    graph: NetworkX graph
    start_nodes: Start nodes list
    sink_nodes: Sink nodes list

    Returns
    -------
    contigs_list: Contig list
    """
    contigs_list = []

    for start in start_nodes:
        for sink in sink_nodes:
            for path in nx.all_simple_paths(graph, source=start, target=sink):
                contig = ""
                for kmer in path:
                    if not contig:
                        contig += kmer
                    else:
                        contig += kmer[-1]
                contigs_list.append((contig, len(contig)))

    return contigs_list


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contig_list, contig_filename):
    """Save all the contigs of a list into a Fasta format.

    Arguments
    ---------
    contig_list: Contig list
    contig_filename: Contif file name

    Returns
    -------
    Fasta file in the working directory
    """
    with open(contig_filename, "w") as f_out:
        for i, contig in enumerate(contig_list):
            f_out.write(f">contig_{i} len={contig[1]}\n")
            f_out.write(fill(contig[0]) + "\n")


def std(value_list):
    """Calculate the standard deviation of a list.

    Arguments
    ---------
    value_list: List of values

    Returns
    -------
    Standard deviation of the list
    """
    return statistics.stdev(value_list)


def path_average_weight(graph, path):
    """Calculate the average weight of a path.

    Arguments
    ---------
    graph: NetworkX graph
    path: Path

    Returns
    -------
    Path average weight
    """
    weight = 0
    for i in range(len(path) - 1):
        weight += graph.edges[path[i], path[i+1]]["weight"]

    return weight / (len(path) - 1)


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove all the path of a graph from a list with or witout removing
    the entry or sinking nodes.

    Arguments
    ---------
    graph: NetworkX graph
    path_list: Path list
    delete_entry_node: Entry node deletion boolean
    delete_sink_node: Sink node deletion boolean

    Returns
    -------
    graph: NetworkX graph
    """
    for path in path_list:
        graph.remove_nodes_from(path[1:-1])
        if delete_entry_node:
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[-1])

    return graph


def select_best_path(graph, path_list, path_length, path_weight,
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path all the path from a list with or witout removing
    the entry or sinking nodes.

    Arguments
    ---------
    graph: NetworkX graph
    path_list: Path list
    path_length: Path length
    path_weight: Path weight
    delete_entry_node: Entry node deletion boolean
    delete_sink_node: Sink node deletion boolean

    Returns
    -------
    graph: NetworkX graph
    """
    # Weights
    path_list_wmax = []
    weight_wmax = []
    length_wmax = []

    for i, path in enumerate(path_list):
        if path_weight[i] == max(path_weight):
            path_list_wmax.append(path)
            weight_wmax.append(path_weight[i])
            length_wmax.append(path_length[i])

    # Length
    path_list_lmax = []
    weight_lmax = []
    length_lmax = []

    for i, path in enumerate(path_list_wmax):
        if length_wmax[i] == max(length_wmax):
            path_list_lmax.append(path)
            weight_lmax.append(weight_wmax[i])
            length_lmax.append(length_wmax[i])

    # Random
    while len(path_list_lmax) > 1:
        path_list_lmax.pop(random.randint(0, len(path_list_lmax)))

    wrong_paths = []

    for path in path_list:
        if path not in path_list_lmax:
            wrong_paths.append(path)

    graph = remove_paths(graph, wrong_paths, delete_entry_node, delete_sink_node)

    return graph


def solve_bubble(graph, a_node, d_node):
    """Select the best path between two nodes of a bubble and remove the others.

    Arguments
    ---------
    graph: NetworkX graph
    a_node: Ancestor node
    b_node: Descendant node

    Returns
    -------
    graph: NetworkX graph
    """
    path_list = []
    path_length = []
    path_weight = []

    for path in nx.all_simple_paths(graph, source=a_node, target=d_node):
        path_list.append(path)
        path_length.append(len(path))
        path_weight.append(path_average_weight(graph, path))

    graph = select_best_path(graph, path_list, path_length, path_weight,
                         delete_entry_node=False, delete_sink_node=False)

    return graph


def simplify_bubbles(graph):
    """Search for bubbles in the graph and simplify them.

    Arguments
    ---------
    graph: NetworkX graph

    Returns
    -------
    graph: NetworkX graph
    """
    wrong_nodes = []

    for d_node in graph.nodes:
        pred_list = list(graph.predecessors(d_node))
        if len(pred_list) > 1:
            a_node = nx.lowest_common_ancestor(graph, pred_list[0], pred_list[1])
            wrong_nodes.append([a_node, d_node])

    for nodes_ad in wrong_nodes:
        graph = solve_bubble(graph, nodes_ad[0], nodes_ad[1])

    return graph


def solve_entry_tips(graph, start_nodes):
    """Remove the non optimal entry tips.

    Arguments
    ---------
    graph: NetworkX graph
    start_nodes: List of starting nodes

    Returns
    -------
    graph: NetworkX graph
    """
    mult_ancestors = []

    for node in start_nodes:
        for next_node in nx.descendants(graph, node):
            if len(graph.pred[next_node]) >= 2 and next_node not in mult_ancestors:
                mult_ancestors.append(next_node)

    paths_list = []
    path_length = []
    path_weight = []

    for node in start_nodes:
        for node2 in mult_ancestors:
            for path in nx.all_simple_paths(graph, node, node2):
                paths_list.append(path)
                path_length.append(len(path))
                path_weight.append(path_average_weight(graph, path))

        graph = select_best_path(graph, paths_list, path_length, path_weight,
                             delete_entry_node=True, delete_sink_node=False)

    return graph


def solve_out_tips(graph, sink_nodes):
    """Remove the non optimal sinking tips.

    Arguments
    ---------
    graph: NetworkX graph
    sink_nodes: List of sinking nodes

    Returns
    -------
    graph: NetworkX graph
    """
    mult_descendants = []

    for node in sink_nodes:
        for next_node in nx.ancestors(graph, node):
            if len(graph.succ[next_node]) >= 2 and next_node not in mult_descendants:
                mult_descendants.append(next_node)

    paths_list = []
    path_length = []
    path_weight = []

    for node in sink_nodes:
        for node2 in mult_descendants:
            for path in nx.all_simple_paths(graph, node2, node):
                paths_list.append(path)
                path_length.append(len(path))
                path_weight.append(path_average_weight(graph, path))

        graph = select_best_path(graph, paths_list, path_length, path_weight,
                             delete_entry_node=False, delete_sink_node=True)

    return graph


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Build the k-mer dictionary
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)

    # Build the associated graph
    graph = build_graph(kmer_dict)

    # Get the starting and sinking node of this graph.
    start_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)

    # Remove all of it's bubbles
    graph = simplify_bubbles(graph)

    # Remove the non optimal entry and out tips
    graph = solve_entry_tips(graph, start_nodes)
    graph = solve_out_tips(graph, sink_nodes)

    # Get the contigs and save it in a Fasta format
    contig_list = get_contigs(graph, start_nodes, sink_nodes)
    save_contigs(contig_list, args.output_file)


if __name__ == '__main__':
    main()

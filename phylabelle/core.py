#!/usr/bin/env python
# coding: utf-8
import warnings
from collections import deque
from functools import partial

from networkx import Graph
from networkx.algorithms.components.connected import connected_component_subgraphs
from networkx.algorithms.matching import max_weight_matching
from sortedcontainers import SortedListWithKey

from phylabelle.fileio import read_tree
from phylabelle.orm import Assembly

# need this sorted list for PhyloTree.evaluate
# this is meant to store tuples of (distance, node)
SortedNodeList = partial(SortedListWithKey, key=lambda x: x[0])


class NoResultsException(Exception):
    def __init__(self):
        pass


def merge_lists(lists, threshold, dist_shift):
    """
    :param iterable lists:
    :param threshold:
    :param dist_shift:
    :return: merged sorted node lists
    :rtype: SortedNodeList
    """
    try:
        i, l_ = lists.popleft()
    except IndexError:
        return SortedNodeList([])

    for i, list_ in lists:
        l_.update(list_)

    delete = False

    for i, val in enumerate(l_):
        val[0] += dist_shift
        if val[0] > threshold:
            delete = True
            break

    if delete:
        del l_[i:]

    return l_


def mean(ls):
    """
    calculates the arithmetic mean from a list on numbers
    :param list ls: list of numbers
    :return arithmetic mean of ls:
    """
    if len(ls) == 0:
        return float('nan')
    else:
        return float(sum(ls)) / len(ls)


def same_species_plain(id1, id2, species):
    """
    return true if id1 belongs to the same species as id2
    :param str id1: assembly id1
    :param str id2: assembly id2
    :param dict species: dictionary of species taxon ids
    :return bool eq: species is the same
    """
    if species[id1] is None:
        warnings.warn('Species taxon id is None. Test for same species can not be done.')
    return species[id1] == species[id2]


def cmp_edges(edge1, edge2):
    """
    return edge with shorter distance
    :param tuple edge1: (accession, {'distance': value})
    :param tuple edge2: (accession, {'distance': value})
    :return tuple: (accession, {'distance' : value})
    """
    if edge1[1]['distance'] < edge2[1]['distance']:
        return edge1
    return edge2


class PairGraph(Graph):
    """
    Inherits from networkx's Graph class. PairGraph-objects are used as an
    intermediate storage structure for found pairs, preserving the network-like
    character of the information.
    """
    def __init__(self, max_size=None, graph=None, mode='all', max_dist=None):
        super(PairGraph, self).__init__()

        self.add_pair = None

        self.set_mode(mode)
        self.max_size = max_size
        self.max_dist = max_dist

        if graph:
            self.edge = graph.edge
            self.node = graph.node

    def resize(self, size=None):
        """
        remove edges (and their corresponding nodes, in case they become
        disconnected from the other components), until self.size() == size.
        if no specific size parameter is provided, self.max_size is used.
        Edges with highest distances are removed first.
        :param int size: size of graph after finishing
        """

        #: do nothing if no specific size is given
        if size is None and self.max_size is None:
            return
        elif size is None:
            size = self.max_size

        n_edges = self.number_of_edges()

        if n_edges <= size:
            return

        #: find out how many edges have to be removed in order to
        #: achieve size as final size after the operation
        size_diff = n_edges - size

        #: store removal candidates in sorted list, sorted by distance
        candidates = SortedListWithKey(key=lambda x: x[2])

        #: iterate over all edges and store candidates with highest
        #: distance value
        for node1, node2 in self.edges():
            dist = self.edge[node1][node2]['distance']
            if len(candidates) < size_diff:
                candidates.add((node1, node2, dist))
            elif candidates[0][2] < dist:
                # remove candidate with lowest distance
                del candidates[0]
                # and add new candidate instead
                candidates.add((node1, node2, dist))

        # remove candidates
        self.remove_edges_from(candidates)

        # get a set of nodes from list of tuples
        nodes = set([node for tup in candidates for node in tup[:2]])

        # now remove nodes
        for node in nodes:
            if len(self.edge[node]) == 0:
                self.remove_node(node)

    def _add_same_spec(self, node1, node2, dist):
        """
        Only adds edge if both nodes belong to the same species
        :param node1:
        :param node2:
        :param dist: distance value
        """
        if self.max_dist is not None and self.max_dist < dist:
            return

        if same_species(node1, node2):
            self.add_edge(node1, node2, distance=dist)

        if self.max_size is not None:
            self.resize()

    def _add_diff_spec(self, node1, node2, dist):
        """
        Only adds edge, if noth nodes belong to different species
        :param node1:
        :param node2:
        :param dist: distance
        """
        if self.max_dist is not None and self.max_dist < dist:
            return

        if not same_species(node1, node2):
            self.add_edge(node1, node2, distance=dist)

        if self.max_size is not None:
            self.resize()

    def _add_all(self, node1, node2, dist):
        """
        Adds edges, regardless of species
        :param node1: node1
        :param node2: node2
        :param dist: distance
        :return:
        """
        if self.max_dist is not None and self.max_dist < dist:
            return

        self.add_edge(node1, node2, distance=dist)
        if self.max_size is not None:
            self.resize()

    def set_mode(self, mode):
        """
        Sets adding mode for a PairGraph object.
        Modes can be all, inter, or intra, where inter stands for interspecies,
        i.e. only pairs consisting of organisms from different species are
        considered valid, or intra respectively.
        """
        assert mode in ['all', 'inter', 'intra']
        if mode == 'all':
            self.add_pair = self._add_all
        elif mode == 'inter':
            self.add_pair = self._add_diff_spec
        elif mode == 'intra':
            self.add_pair = self._add_same_spec

    def get_connected_components(self):
        """
        iterate over connected components
        :return: iterator of PairGraphs
        """
        for subgraph in connected_component_subgraphs(self):
            yield subgraph

    def average_distance(self):
        """
        :return: arithmethic mean of distances
        """
        return mean([partner['distance'] for node in self.edge for partner in node.values()])

    def evaluate(self):
        """
        print some information about connected components
        :return:
        """
        n_pairs = []
        n_nodes = []
        average_distances = []

        for component in self.get_connected_components():
            print component.number_of_edges(), component.number_of_nodes(), \
                component.average_distance()
            n_pairs.append(component.size())
            n_nodes.append(len(component.nodes()))
            average_distances.append(component.average_distance())

        return mean(n_pairs), mean(n_nodes), mean(average_distances)

    def get_closest_partner(self, node):
        """
        :param node: node
        :return: closest partner
        """
        try:
            return reduce(cmp_edges, self.edge[node].items())
        except TypeError:
            return None

    def get_partners(self, query_node):
        """
        return partners in same format as in get_closest and get_minimum_matching
        :return:
        """
        for node2 in self.edge[query_node]:
            dist = self.edge[query_node][node2]['distance']
            yield tuple(sorted([query_node, node2]) + [dist])

    def get_closest(self):
        """
        get minimum distance edges
        """
        for node1, node2, data in self.edges_iter(data=True):
            yield tuple(sorted([node1, node2]) + [data['distance']])

    def get_minimum_matching(self):
        """
        calculate a minimum weight matching with maximum cardinality.
        Therefore invert weights by calculating 1/w for ever weight w, and then
        use networkx.algorithms.matching.max_weight_matching
        """

        subgraphs = self.get_connected_components()

        for sg in subgraphs:
            for node1, node2 in sg.edges_iter():
                dict_ = sg[node1][node2]
                dict_['weight'] = 1/dict_['distance']

            seen = []

            for edge1, edge2 in max_weight_matching(sg).items():
                if edge1 in seen or edge2 in seen:
                    continue
                dist = sg.edge[edge1][edge2]['distance']
                seen.append(edge1)
                seen.append(edge2)
                yield tuple(sorted([edge1, edge2]) + [dist])


class PhyloTree(object):
    def __init__(self, tree, session, **kwargs):
        """
        :param str tree: path to tree file
        :param session: database-session
        :param float threshold:   a numerical value indicating below which
                                  distance a pair should be considered as such
        """
        assert type(tree) == str

        self.root = read_tree(tree)

        assert len(self.root.children) == 2, 'Trying to work with unrooted Tree'

        self.species_index = Assembly.get_species_index(session)

        global same_species
        same_species = partial(same_species_plain, species=self.species_index)

        self._rename_leaves()
        self.assemblies = self.load_assemblies(session=session)
        self.labels = {k: v.label for k, v in self.assemblies.iteritems()}
        self._assign_labels()
        if 'threshold' in kwargs:
            self.results = PairGraph(max_dist=kwargs['threshold'])
        else:
            self.results = PairGraph()

    def _assign_labels(self):
        """
        Assign labels to as attributes leaves
        """
        count = 0

        for leaf in self.root.iter_leaves():
            try:
                leaf.label = self.labels[leaf.name]
                count += 1
            except KeyError:
                pass

    def _format_results(self, results):
        """
        provide output as sorted list of tuples (node1, node2, distance)
        :param results:
        :return:
        """
        return SortedListWithKey(((self.assemblies[node1], self.assemblies[node2], dist)
                                  for node1, node2, dist in results),
                                 key=lambda x: x[2])

    def load_assemblies(self, session):
        leaves = [x.name for x in self.root.iter_leaves()]
        asm_index = {x.accession: x for x in Assembly.iter_all(session, subset=leaves)}
        return asm_index

    def get_closest(self):
        """
        calls self.results.get_closest() and adds Assembly objects to plain accessions from graph
        :return list results: [(assembly1, assembly2, distance), ...]
        """
        if not len(self.results.edge) > 0:
            raise NoResultsException()

        return self._format_results(self.results.get_closest())

    def get_minimum_matching(self):
        return self._format_results(self.results.get_minimum_matching())

    def _rename_leaves(self):
        for leaf in self.root.iter_leaves():
            leaf.name = '_'.join(leaf.name.split('_')[:2])

    def evaluate_all_pairs(self, distance,
                           mode='all'):

        # self.pair_graph = PairGraph()
        self.results = PairGraph(mode=mode, max_dist=distance)
        self.evaluate(self.root, distance)

    def get_connected_components(self):
        return connected_component_subgraphs(self.results)

    def get_distance(self, pair):
        """
        get distance for a pair object
        """
        l1 = self.root.get_leaves_by_name(pair.np_acc)[0]
        l2 = self.root.get_leaves_by_name(pair.p_acc)[0]
        return self.root.get_distance(l1, l2)

    def leaves_unique(self, node=None):
        """
        Sanity check function to check if the leaves are ambiguous or not
        :param node: root node
        :return bool:
        """
        if node is None:
            node = self.root
        leaves = [x.name for x in node.iter_leaves()]
        s = len(leaves)
        leaves = set(leaves)
        t = len(leaves)
        return s == t

    def find_leaf(self, leaf_name, node=None):
        """
        find leaf by name
        :param node:
        :param leaf_name:
        :return:
        """
        if node is None:
            node = self.root

        for leaf in node.iter_leaves():
            if leaf_name in leaf.name:
                return leaf

    def find_closest_partner(self, name, mode='all', n_results=1):
        # clear buffer
        self.results.clear()
        self.results.max_size = n_results
        self.results.resize()
        self.results.set_mode(mode)

        # find leaf
        l = self.find_leaf(name)

        # init algorithm
        inf = float('inf')

        node = l.up
        evals = {}

        while self.results.number_of_edges() < n_results:
            evals[node] = self.evaluate(node, inf, evaluated=evals, query=l.name)
            node = node.up
            if not node:
                break

        return self._format_results(self.results.get_partners(l.name))

    def evaluate(self, node, threshold, evaluated=None, query=None):
        """
        :param node:
        :param threshold:
        :param evaluated:
        :param query:
        :return:
        """
        if not evaluated:
            evaluated = {}

        try:
            if node.label:
                return SortedNodeList([[node.dist, node]]), SortedNodeList()
            else:
                return SortedNodeList(), SortedNodeList([[node.dist, node]])
        except:
            pass

        true_ = deque([])
        false_ = deque([])

        if evaluated:
            for i, child in enumerate(node.children):
                if child in evaluated:
                    t_, f_ = evaluated[child]
                else:
                    t_, f_ = self.evaluate(child, threshold, query=query)

                true_.append((i, t_))
                false_.append((i, f_))
        else:
            for i, child in enumerate(node.children):
                t_, f_ = self.evaluate(child, threshold, query=query)
                true_.append((i, t_))
                false_.append((i, f_))

        for i, t_ in true_:
            for j, f_ in false_:
                if i == j:
                    continue
                for t_dist, t_node in t_:
                    for f_dist, f_node in f_:
                        d = t_dist + f_dist
                        if d < threshold:
                            if query is not None and query not in [f_node.name, t_node.name]:
                                continue

                        self.results.add_pair(f_node.name, t_node.name, d)

        t_ = merge_lists(true_, threshold, node.dist)
        f_ = merge_lists(false_, threshold, node.dist)

        return t_, f_

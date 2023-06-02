from moleculekit.molecule import Molecule
from moleculekit.vmdviewer import viewer
import numpy as np
from scipy.spatial.distance import cdist
from collections import Counter
import sys, getopt
from graph_tool.all import *
import math
from typing import Tuple, Dict

import logging
import moleculekit

logging.getLogger(moleculekit.molecule.__name__).setLevel(logging.WARNING)
logger = logging.getLogger('hydrophobic_clusters')
from typing import NamedTuple


class Cluster(NamedTuple):
    area: np.array
    residues: list
    chains: list
    contacts: int
    ratio_contacts_residue: float
    ratio_area_residue: float


# Paths and variables
sel = "protein and chain A and not backbone and noh and resname ILE VAL LEU"
_ATOMIC_RADII = {'C': 1.88}
water_radius = 1.4
sphere_radius_carbon = _ATOMIC_RADII['C'] + water_radius
sphere_points = 610
sphere_area_const = 4.0 * math.pi * (sphere_radius_carbon ** 2) / sphere_points


class Atom:
    """
    Class to handle ILV atoms
    index: the vmd atom index
    mol: an htmd.htmd.molecule object
    """
    radius = _ATOMIC_RADII['C']  # All distances in Angstrom

    def __init__(self, index: int, mol: Molecule) -> None:
        self.index = index
        self.coords = mol.coords[index][:, 0]
        self.resid = mol.resid[index]
        self.point_coords = generate_sphere_points(self.coords, sphere_points, self.radius)
        self.neighbor_indices = self.get_neighbors(mol)

    def get_neighbors(self, mol: Molecule) -> np.array:
        """
        Provides all indices of atoms within 6.56 A of this atom.
        6.56 is the upper bound of a possible neighbor 1.88 (C) + 1.4 + 1.4 + 1.88 (C).
        """
        neighbor_indices = mol.get("index", sel=f"protein\
        and noh and not resid '{self.resid}' and within 6.56 of index '{self.index}'")
        return neighbor_indices


def generate_sphere_points(coords: np.array, n: int = 610, radius: float = 1.88) -> np.array:
    """
    :param coords: The coordinates of the center of the atom
    :param n: number of points to sample
    :param radius: the radius of the atom
    :return: a nx3 vector with the point coordinates
    """
    total_radius = radius + water_radius
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y * y)
        phi = k * inc
        points.append([math.cos(phi) * r, y, math.sin(phi) * r])
    vec = np.asarray(points)
    vec *= total_radius
    vec += coords
    return vec


def retrieve_neighbor_positions(atom: Atom, mol: Molecule) -> Tuple[np.array, Dict[int, int]]:
    """
    :param atom: an Atom object
    :param mol: a htmd.htmd.molecule object
    :return: A tuple object with the positions of the neighboring atoms. A dictionary indexing column positions to resid positions
    """
    positions = mol.coords[atom.neighbor_indices][:, :, 0]
    position_index_to_resid = {index: mol.resid[neighbor_indice] for index, neighbor_indice in
                               enumerate(atom.neighbor_indices)}
    return positions, position_index_to_resid


def retrieve_indices(matrix: np.array, coords: np.array, neighborpositions: np.array, radius: float = 1.88) -> np.array:
    """
    Computes if each of the n sphere points are penetrating neighboring spheres
    :param matrix: n x m Distance matrix where n is the number of sphere points and m the number of neighbors
    :param coords: the coordinates of the atom
    :param neighborpositions: Coordinates of the neighbors
    :param radius: radius of the atom
    :return: The atoms that are in closest distance with each n sphere points.
    """

    ## When a row contains several atoms in contact,
    #  select the one with the closest center to the center of atom A
    sphere_radius = water_radius + radius
    #When a point is within the sphere of other atoms, the atom that is closest is chosen.
    # We need to compute the distances between all atom-neighbor pairs.
    dist_center_atoms = cdist(np.reshape(coords, (1, 3)), neighborpositions)
    ranking = np.argsort(dist_center_atoms)
    valid = matrix <= sphere_radius
    idx2 = []
    for row in valid:
        if row.any():
            # 1. ordering the row with Falses and Trues according to the ranking order
            # 2. Getting the first element that is true, which will be the first in the
            # ranking list, thus the index of the closest atom.
            idx2.append(ranking[0][np.where(np.isin(ranking, np.where(row)))[1][0]])
    return idx2


def fill_matrices(atom: Atom, mol: Molecule,
                  resid_matrix: np.array, indices: np.array, atom_to_residposition) -> Tuple[np.array, np.array]:
    """
    :param atom: An Atom class
    :param mol: an htmd.htmd.molecule object
    :param atom_matrix: the index x index area matrix
    :param resid_matrix: the ILVresid x ILVresid area matrix
    :param indices: the indices that belong to ILV sidechain heavy atoms
    :param resids: the resids that belong to ILV sidechain heavy atoms
    :return: Updated atom_matrix and resid_matrix
    """
    neighbor_positions, position_index_to_resid = retrieve_neighbor_positions(atom, mol)
    # Compute distances between all sphere points and the neighbors.
    # THe shape will be 610 (rows) x nr. of neighbors (columns).
    distances = cdist(atom.point_coords, neighbor_positions)

    #a list with the atoms in closest distance to each of the n (610) points
    column_indices = retrieve_indices(distances, atom.coords, neighbor_positions)
    colpos_occurrences = Counter(column_indices)

    for colpos, occurrences in colpos_occurrences.items():
        if atom.neighbor_indices[colpos] in indices:
            area = sphere_area_const * occurrences
            index_i = atom_to_residposition[atom.index]
            index_j = atom_to_residposition[atom.neighbor_indices[colpos]]
            resid_matrix[index_i, index_j] += area
    return resid_matrix


def create_graph(resid_matrix: np.array, resid_list: np.array, chain_list: np.array, cutoff_area: float = 10.0) -> Graph:
    """
    :param resid_matrix: the ILVresid x ILVresid area matrix
    :param resid_list: the index x index area matrix
    :return: A Graph object where each component is a ILV cluster
    """
    g = Graph()
    g.vp.resid = g.new_vertex_property("int")
    g.vp.chain = g.new_vertex_property("string")
    g.ep.area = g.new_edge_property("float")

    # 1. Create all vertices
    for v in range(len(resid_matrix)):
        v1 = g.add_vertex()
        g.vp.resid[v1] = resid_list[v]
        g.vp.chain[v1] = chain_list[v]

    # 2. Create edges and fill its values with areas
    for row_index, row in enumerate(resid_matrix):
        v1 = g.vertex(row_index)
        for column_index, area in enumerate(row):
            v2 = g.vertex(column_index)
            if area > cutoff_area:
                ae = g.add_edge(v1, v2)
                g.ep.area[ae] = area
    return g

def add_clusters(g: Graph, components: PropertyArray):
    """
    :param mol:
    :param g:
    :param components:
    :return: Molecule representations
            A list with Cluster objects
    """
    clusters = []

    for cluster_index in range(max(components) + 1):
        cluster = [i for i, x in enumerate(components) if x == cluster_index]
        if len(cluster) < 2: continue
        vfilt = g.new_vertex_property('bool')
        for i in cluster: vfilt[i] = True
        sub = GraphView(g, vfilt)
        area = np.sum([g.ep.area[edge] for edge in sub.edges()])

        resid_cluster = [g.vp.resid[i] for i in cluster]
        chain_cluster = [g.vp.chain[i] for i in cluster]

        clusters.append(Cluster(
            area = area,
            residues=resid_cluster,
            chains=chain_cluster,
            contacts=sub.num_edges(),
            ratio_contacts_residue=sub.num_edges() / len(cluster),
            ratio_area_residue=area / sub.num_edges()
        ))

    return clusters
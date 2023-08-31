'''
Created on Jun 8, 2017

@author: Matthew Manning
'''

import collections
import itertools
import math
import math as m

import NanoParticle
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


# def herons_formula(acoords, bcoords, ccoords):
#    sides = [0]*3
#    for i, (coords1, coords2) in enumerate(itertools.combinations([acoords, bcoords, ccoords], 2)):
#        sides[i] = [leaptotalbuilder.calculate_distance(coords1, coords2)
#    S = sum(sides) / 2.
#    area = np.sqrt(S * (S - A) * (S - B) * (S -C))
#    return area
# from nanorodbatchrefactor import __surfthick__, diag_add, capangle

def join_tuples_flatly(tuple_of_tuples):
    return tuple(j for i in tuple_of_tuples for j in (i if isinstance(i, tuple) else (i,)))


def get_nested_depth(sequence):
    from collections import Sequence
    from itertools import chain, count
    sequence = iter(sequence)
    try:
        for level in count():
            sequence = chain([next(sequence)], sequence)
            sequence = chain.from_iterable(s for s in sequence if isinstance(s, Sequence))
    except StopIteration:
        return level


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


# def surf_area_points(coordslist):
#    for
def rotation_matrix_from_vectors(current_vector, target_vector):
    if all(current_vector == target_vector):
        return np.identity(3)
    # Orient ligand along its long principal axis.
    v = np.cross(current_vector, target_vector)
    # sin_angle = np.linalg.norm(v)
    cos_angle = np.dot(current_vector, target_vector)
    v_skew = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    raw_rotation_matrix = np.identity(3) + v_skew + np.matmul(v_skew, v_skew) * (1 / (1 + cos_angle))
    # Normalize matrix
    rotation_matrix = raw_rotation_matrix / np.linalg.norm(raw_rotation_matrix)
    # rotation_matrix = raw_rotation_matrix / row_sums
    return rotation_matrix


def lattice_builder_from_ranges(lattice_params, min_max_tuple, lattice_type='sc'):
    cells_tuple = []
    # for cell_num, param in zip(unit_cells, lattice_params):
    #    cells_tuple.append([param * coord for coord in list(range(cell_num))])
    coords_tuple = []
    for dimension, param in zip(min_max_tuple, lattice_params):
        # print(dimension, param)
        cells_tuple.append(np.arange(dimension[0], dimension[1], step=param))
    # print(cells_tuple)
    for x, y, z in itertools.product(cells_tuple[0], cells_tuple[1], cells_tuple[2]):
        coords_tuple.append([x, y, z])
    return coords_tuple


def convert_from_simple_cubic(coords_tuple, lattice_type, lattice_param):
    expanded_coords_tuple = []
    if lattice_type == 'bcc':
        fractional_coords_tuple = np.array(
            [[0., 0., 0.], [0.5 * lattice_param[0], 0.5 * lattice_param[0], 0.5 * lattice_param[0]]])
    elif lattice_type == 'fcc':
        fractional_coords_tuple = np.array([[0., 0., 0.], [.5, .5, 0.], [.5, 0., .5], [0., .5, .5]])
    for coords, fractional_coords in itertools.product(coords_tuple, fractional_coords_tuple):
        expanded_coords_tuple.append([coordinate + move for coordinate, move in zip(coords, fractional_coords)])
    return expanded_coords_tuple


def get_taxicab_distance(coords1, coords2):
    distance = abs(coords1[0] - coords2[0]) + abs(coords1[1] - coords2[1]) + abs(coords1[2] - coords2[2])
    return distance


def calculate_distance(coords1, coords2):
    # x1, y1, z1 = float(coords1[0]), float(coords1[1]), float(coords1[2])
    # x2, y2, z2 = float(coords2[0]), float(coords2[1]), float(coords2[2])
    x1, y1, z1 = coords1
    x2, y2, z2 = coords2
    distance = np.sqrt(float((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2))
    if type(coords1[0]) is int: distance = int(distance)
    return distance


def calculate_spherical_contour_length(coords1, coords2, center=[0., 0., 0.]):
    # float_vector1 = [float(coords1[0]), float(coords1[1]), float(coords1[2])] - center
    # float_vector2 = [float(coords2[0]), float(coords2[1]), float(coords2[2])] - center
    radius1 = calculate_distance(coords1, center)
    radius2 = calculate_distance(coords2, center)
    average_radius = (radius1 + radius2) / 2
    if 2 * (radius1 - radius2) / (radius1 + radius2) > 0.100: print(
        "Mismatch in calculated radii: Check that particle is sphere centered on {}".format(center))
    radians_between_vectors = angle_between(coords1, coords2)
    contour_length = average_radius * radians_between_vectors / m.pi()
    # contour_length = average_radius * radians_between_vectors
    return contour_length


def contour_length_pairs(coords_list):
    contour_pairs = np.zeros([len(coords_list), len(coords_list)])
    # print dist_pairs
    # print(list(itertools.combinations(enumerate(coords_list), 2))[0])
    for (ind1, coords1), (ind2, coords2) in itertools.combinations(enumerate(coords_list), 2):
        if ind1 < ind2:
            contour_pairs[ind1, ind2] = calculate_spherical_contour_length(coords1, coords2)
        else:
            contour_pairs[ind2, ind1] = calculate_spherical_contour_length(coords1, coords2)
    return contour_pairs


def distance_pairs(coords_tuple, same_index_distance=0.000001):
    coords_matrix = np.array(coords_tuple)
    repeat_rows = np.broadcast_to(coords_matrix, (len(coords_tuple), len(coords_tuple), 3))
    repeat_columns = np.swapaxes(repeat_rows, 0, 1)
    distance_pairs_array = np.sqrt(np.sum(np.square(repeat_rows - repeat_columns), axis=2).astype(float))

    '''
    dist_pairs_array = np.empty([len(coords_tuple), len(coords_tuple)], dtype=ligand_type(coords_tuple[0][0]))
    for ind1, coords1 in enumerate(coords_tuple):
        for ind2, coords2 in enumerate(coords_tuple):
            if ind1 == ind2:
                dist_pairs_array[ind1, ind2] = same_index_distance
            elif ind1 < ind2:
                dist_pairs_array[ind1, ind2] = calculate_distance(coords1, coords2)
                dist_pairs_array[ind2, ind1] = calculate_distance(coords1, coords2)
            else:
                continue
    '''
    return distance_pairs_array


def plot_rad_dist_func_from_dist_array(coords_list):
    sns.set()
    distance_array = distance_pairs(coords_list)
    ax = sns.distplot(np.flatnonzero(distance_array))
    plt.show()
    return ax


def plot_points_on_sphere_wireframe(coords_tuple, radius):
    # Great circles for sphere wireframe.
    phi = np.linspace(0, np.pi, 20)
    theta = np.linspace(0, 2 * np.pi, 40)
    x = radius * np.outer(np.sin(theta), np.cos(phi))
    y = radius * np.outer(np.sin(theta), np.sin(phi))
    z = radius * np.outer(np.cos(theta), np.ones_like(phi))
    # Convert nested list of coords into numpy array.
    coords_array = np.array(coords_tuple)
    # Plot points.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')  # subplot_kw={'projection': '3d', 'aspect': 'equal'})
    ax.plot_wireframe(x, y, z, color='k', rstride=1, cstride=1)

    ax.scatter(coords_array[:, 0], coords_array[:, 1], coords_array[:, 2], s=5000, c='k', zorder=10)
    ax.set_aspect("equal")
    plt.show()
    return fig, ax


def plot_nanoparticle_with_ligands(surface_coords_list, uncharged_index_list, charged_index_list, radius):
    atom_size = radius
    surface_coords_array = np.array(surface_coords_list) / 1000.
    if len(uncharged_index_list) > 0:
        uncharged_coords_array = np.take(surface_coords_array, np.array(uncharged_index_list), axis=0)
    if len(charged_index_list) > 0:
        charged_coords_array = np.take(surface_coords_array, np.array(charged_index_list), axis=0)
    if len(uncharged_index_list) == 0 and len(charged_index_list) == 0:
        uncharged_coords_array = list(range(len(surface_coords_list)))
    unbound_coords_array = np.delete(surface_coords_array, np.array(uncharged_index_list + charged_index_list), axis=0)

    # Plot points.
    # fig, ax = plt.subplots(1, 1, subplot_kw={'projection': '3d', 'aspect': 'equal'})
    fig = plt.figure()
    # ax = Axes3D(fig)
    ax = fig.add_subplot(111, projection='3d')
    if unbound_coords_array.size > 0: ax.scatter(unbound_coords_array[:, 0], unbound_coords_array[:, 1],
                                                 unbound_coords_array[:, 2], s=atom_size, c='y', marker='.')
    if len(uncharged_index_list) > 0: ax.scatter(uncharged_coords_array[:, 0], uncharged_coords_array[:, 1],
                                                 uncharged_coords_array[:, 2], s=atom_size, c='bk', marker='.')
    if len(charged_index_list) > 0: ax.scatter(charged_coords_array[:, 0], charged_coords_array[:, 1],
                                               charged_coords_array[:, 2], s=atom_size, c='b', marker='.')
    ax.set_aspect("equal")
    plt.tight_layout()
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()
    return fig, ax


def get_number_of_bonds_from_pairs_tuple(bond_pairs_tuple, number_of_atoms):
    number_of_bonds_list = [0] * number_of_atoms
    for index1, index2 in bond_pairs_tuple:
        number_of_bonds_list[index1] += 1
        number_of_bonds_list[index2] += 1
    return number_of_bonds_list


def get_nearest_neighbors(coords_tuple=[], rough_distance_matrix=[], nearest_neighbor_distance=0, rounding_decimal=2,
                          tolerance=1.05):
    # print(coords_tuple)
    if rough_distance_matrix == [] and coords_tuple != []:
        rough_distance_matrix = distance_pairs(coords_tuple)
    elif rough_distance_matrix == [] and coords_tuple == []:
        print(
            'Either Numpy array of pairwise distances(distance_matrix) or tuple of Cartesian coordinates(coords_tuple) required.')
        raise IOError
    distance_matrix = rough_distance_matrix.round(rounding_decimal)
    nearest_neighbor_distance = np.nanmin(distance_matrix.ravel()[np.flatnonzero(distance_matrix)])
    nearest_neighbor_indices_tuple = []
    for np_row in distance_matrix:
        nearest_neighbor_indices_tuple.append(np.where(np_row < tolerance * nearest_neighbor_distance)[0].tolist())
        # print(np_row)
    return nearest_neighbor_indices_tuple


def slow_nearest_neighbors(coords_tuple, rounding_decimal=2, tolerance=1.05):
    rough_distance_matrix = distance_pairs(coords_tuple)
    distance_matrix = rough_distance_matrix.round(rounding_decimal)
    nearest_neighbor_distance = np.nanmin(distance_matrix.ravel()[np.flatnonzero(distance_matrix)])
    nearest_neighbor_indices_tuple = []
    for np_row in distance_matrix:
        nearest_neighbor_indices_tuple.append(np.where(np_row < tolerance * nearest_neighbor_distance)[0].tolist())
        # print(np_row)
    return nearest_neighbor_indices_tuple


def get_num_of_nearest_neighbors(nearest_neighbor_indices=[], coords_tuple=[], distance_pairs_np_array=[]):
    if len(nearest_neighbor_indices) <= 0:
        nearest_neighbor_indices = get_nearest_neighbors(coords_tuple=coords_tuple,
                                                         rough_distance_matrix=distance_pairs_np_array)
    num_of_nearest_neighbors = []
    for row_of_indices in nearest_neighbor_indices:
        num_of_nearest_neighbors.append(len(row_of_indices))
    return num_of_nearest_neighbors


def check_if_surface_atom(min_bulk_coord_num=10, nearest_neighbor_indices=[], coords_tuple=[],
                          distance_pairs_np_array=[]):
    num_of_nearest_neighbors = get_num_of_nearest_neighbors(nearest_neighbor_indices=nearest_neighbor_indices,
                                                            coords_tuple=coords_tuple,
                                                            distance_pairs_np_array=distance_pairs_np_array)
    surface_bool_list = [coord_num < min_bulk_coord_num for coord_num in num_of_nearest_neighbors]
    return surface_bool_list


def get_surface_atom_indices(min_bulk_coord_num=10, nearest_neighbor_indices=[], coords_tuple=[],
                             distance_pairs_np_array=[]):
    surface_bool_list = check_if_surface_atom(min_bulk_coord_num=min_bulk_coord_num,
                                              nearest_neighbor_indices=nearest_neighbor_indices,
                                              coords_tuple=coords_tuple,
                                              distance_pairs_np_array=distance_pairs_np_array)
    surface_index_list = [ind for ind, boo in enumerate(surface_bool_list) if boo]  # surface_bool_list.index(True)
    return surface_index_list


def simple_get_surface_indices(flat_num_nearest_neighbors, min_bulk_coord_num=10):
    surface_indices_list = []
    for flat_index, num_neighbors in enumerate(flat_num_nearest_neighbors):
        if num_neighbors < min_bulk_coord_num:
            surface_indices_list.append(flat_index)
    return surface_indices_list


def pick_from_neighbors(central_atom_index, neighbor_index_list, num_bonds_list, bonded_index_pairs_tuple,
                        number_of_bonds=1, max_bonds_central=6, max_bonds_neighbor=7):
    ### Sorts neighbor_index_list by number of bonds and bonds atom with least number of bonds first.
    sort_index_list = np.argsort([num_bonds_list[neighbor_index] for neighbor_index in neighbor_index_list])
    # print(sort_index_list)
    for neighbor_iter, sort_index in zip(list(range(number_of_bonds)), sort_index_list):
        neighbor_index = neighbor_index_list[sort_index]
        if num_bonds_list[central_atom_index] >= max_bonds_central:
            break
        if num_bonds_list[neighbor_index] < max_bonds_neighbor:
            bonded_index_pairs_tuple.append(sorted([central_atom_index, neighbor_index]))
            num_bonds_list[neighbor_index] += 1
            num_bonds_list[central_atom_index] += 1
    return num_bonds_list, bonded_index_pairs_tuple


def bond_atoms_universal(central_index_list, available_index_list=[], bond_within_central=False,
                         distance_pairs_matrix=[], nearest_neighbor_index_tuple=[], coords_tuple=[],
                         max_central_bonds=7, max_neighbor_bonds=7, bond_num_list=[], bonded_index_pairs_tuple=[]):
    if len(bond_num_list) == 0:
        if len(coords_tuple) > 0:
            bond_num_list = [0] * len(coords_tuple)
        else:
            bond_num_list = [0] * (len(central_index_list) + len(available_index_list))
    if len(nearest_neighbor_index_tuple) == 0:
        nearest_neighbor_index_tuple = get_nearest_neighbors(coords_tuple=coords_tuple,
                                                             rough_distance_matrix=distance_pairs_matrix)
        # distance_pairs_matrix = distance_pairs(coords_tuple)
    for row_num, atom_row in enumerate(nearest_neighbor_index_tuple):
        nearest_neighbor_index_tuple[row_num] = [index for index in atom_row if index in available_index_list]
    for central_index in central_index_list:
        central_num_bonds = bond_num_list[central_index]
        if central_num_bonds >= max_central_bonds:
            continue
        if not bond_within_central:
            neighbor_index_list = [index for index in nearest_neighbor_index_tuple[central_index] if
                                   index not in central_index_list]
        else:
            neighbor_index_list = nearest_neighbor_index_tuple[central_index]
        # print(nearest_neighbor_index_tuple)
        neighbor_index_list = [index for index in neighbor_index_list if not index == central_index and [index,
                                                                                                         central_index] not in bonded_index_pairs_tuple and [
                                   central_index, index] not in bonded_index_pairs_tuple]
        new_bond_num_list, new_bonded_index_pairs_tuple = pick_from_neighbors(central_index, neighbor_index_list,
                                                                              bond_num_list,
                                                                              bonded_index_pairs_tuple,
                                                                              max_bonds_neighbor=max_neighbor_bonds)
    return new_bonded_index_pairs_tuple, new_bond_num_list


def bond_by_groups(index_list_tuple, distance_pairs_matrix=[], nearest_neighbor_index_tuple=[], coords_tuple=[],
                   max_central_bonds=7, max_neighbor_bonds=7, bond_num_list=[], bonded_index_pairs_tuple=[]):
    for group_index, central_index_list in enumerate(index_list_tuple):
        available_index_list = join_tuples_flatly(index_list_tuple[group_index + 1:])
        new_bonded_index_pairs_tuple, bond_num_list = bond_atoms_universal(central_index_list,
                                                                           available_index_list=available_index_list,
                                                                           distance_pairs_matrix=distance_pairs_matrix,
                                                                           nearest_neighbor_index_tuple=nearest_neighbor_index_tuple,
                                                                           max_central_bonds=max_central_bonds,
                                                                           max_neighbor_bonds=max_neighbor_bonds,
                                                                           bond_num_list=bond_num_list,
                                                                           bonded_index_pairs_tuple=bonded_index_pairs_tuple)
    return new_bonded_index_pairs_tuple, bond_num_list


def atom_bonder_with_surface(coords_tuple, bound_index_list, nearest_neighbor_distance, max_num_bonds=7,
                             neighbor_dist_tolerance=1.1):
    # Generate list of indices and shuffle, then loop through indices, leaving them in sequential order.
    # If atoms are not in spatially sequential order, Amber will throw a "fixatomorder/setMolecules" error.
    bond_num_list = []
    bonded_index_pairs_tuple = []
    nearest_neighbor_distance = nearest_neighbor_distance * neighbor_dist_tolerance
    unbound_index_list = [index for index in list(range(len(coords_tuple))) if index not in bound_index_list]
    # if distance_pairs_matrix.size == 0:
    distance_pairs_matrix = distance_pairs(coords_tuple)
    nearest_neighbor_index_tuple = get_nearest_neighbors(rough_distance_matrix=distance_pairs_matrix,
                                                         rounding_decimal=0)
    bonded_index_pairs_tuple, bond_num_list = bond_atoms_universal(bound_index_list,
                                                                   available_index_list=unbound_index_list,
                                                                   nearest_neighbor_index_tuple=nearest_neighbor_index_tuple,
                                                                   max_central_bonds=6,
                                                                   max_neighbor_bonds=6,
                                                                   bond_num_list=[0] * len(coords_tuple))

    number_of_bonds = get_number_of_bonds_from_pairs_tuple(bonded_index_pairs_tuple, len(coords_tuple))
    for count1, count2 in zip(bond_num_list, number_of_bonds):
        if count1 == count2:
            continue
            # else:
            # print('Discrepancy in number of bonds.')
    bond_num_list = number_of_bonds
    core_neighbor_index_tuple = []
    for row in nearest_neighbor_index_tuple:
        core_neighbor_index_tuple.append([index for index in row if index not in bound_index_list])
    bonded_index_pairs_tuple, bond_num_list = bond_atoms_universal(unbound_index_list,
                                                                   available_index_list=unbound_index_list,
                                                                   bond_within_central=True,
                                                                   nearest_neighbor_index_tuple=core_neighbor_index_tuple,
                                                                   bonded_index_pairs_tuple=bonded_index_pairs_tuple,
                                                                   bond_num_list=bond_num_list)

    num_bonds_list_ligand_bound = [bond_num_list[ligand_index] for ligand_index in bound_index_list]
    print('Number of total bonds: {}'.format(len(bonded_index_pairs_tuple)))
    print("Distribution of coordination numbers for ligand-bound atoms: {}".format(
        collections.Counter(num_bonds_list_ligand_bound)))
    print("Distribution of coordination numbers: {}".format(collections.Counter(bond_num_list)))

    return bonded_index_pairs_tuple.copy(), bond_num_list.copy()
    # Bond all ligand sites, first to surface atoms, then to available core atoms.


'''
    for bound_index in bound_index_list:
        bound_coords = coords_tuple[bound_index]
        neighbor_index_list = []
        if num_bonds_list[bound_index] >= max_num_bonds - 1: continue
        for check_index in unbound_index_list:
            neighbor_coords = coords_tuple[check_index]
            if num_bonds_list[check_index] >= max_num_bonds:
                continue
            pair_distance = calculate_distance(bound_coords, neighbor_coords)
            if pair_distance - nearest_neighbor_distance < 0.001:
                neighbor_index_list.append(check_index)
        num_bonds_list, bonded_index_pairs_tuple = pick_from_neighbors(bound_index, neighbor_index_list, num_bonds_list, bonded_index_pairs_tuple, number_of_bonds=6, max_bonds_central=6, max_bonds_neighbor=7)

    for atom_double_index, unbound_index in enumerate(unbound_index_list):
        if num_bonds_list[unbound_index] >= max_num_bonds: continue
        unbound_coords = coords_tuple[unbound_index]
        neighbor_index_list = []
        #for check_index in unbound_index_list[atom_double_index + 1:]:
        for check_index in unbound_index_list[atom_double_index + 1:]:
            if num_bonds_list[check_index] >= max_num_bonds or unbound_index == check_index:
                continue
            neighbor_coords = coords_tuple[check_index]
            pair_distance = calculate_distance(unbound_coords, neighbor_coords)
            if pair_distance - nearest_neighbor_distance < 0.001:
                neighbor_index_list.append(check_index)
        # Continue if neighbor_index_list is empty.
        if neighbor_index_list == []: continue
        ### Sorts neighbor_index_list by number of bonds and bonds atom with least number of bonds first.
        num_bonds_list, bonded_index_pairs_tuple = pick_from_neighbors(unbound_index, neighbor_index_list, num_bonds_list, bonded_index_pairs_tuple, number_of_bonds=7, max_bonds_central=7, max_bonds_neighbor=7)

        sort_list = np.argsort(neighbor_bond_nums)
        for sort_index in sort_list:
            neighbor_index = neighbor_index_list[sort_index]
            if num_bonds_list[bound_index] > 6: break
            bonded_pair_index_list.append([bound_index, neighbor_index])
            num_bonds_list[neighbor_index] += 1
            num_bonds_list[bound_index] += 1

    for bonded_index_pair in bonded_index_pairs_tuple:
        if [bonded_index_pair[1], bonded_index_pair[0]] in bonded_index_pairs_tuple and [
            bonded_index_pair[0], bonded_index_pair[1]] in bonded_index_pairs_tuple:
            print('Duplicate\n')
    # print(num_bonds_list)
    # print(collections.Counter(num_bonds_list))
    num_bonds_list_ligand_bound = [num_bonds_list[ligand_index] for ligand_index in bound_index_list]
    print("Distribution of coordination numbers for ligand-bound atoms: {}".format(collections.Counter(num_bonds_list_ligand_bound)))
    print("Distribution of coordination numbers: {}".format(collections.Counter(num_bonds_list)))
    print('Number of bonds: {}.'.format(len(bonded_index_pairs_tuple)))
    return bonded_index_pairs_tuple
'''
if __name__ == '__main__':
    pass


# Takes a 3-D numpy array (unit_cell_3d_array) consisting of 2-D numpy arrays, where each row is the x,y,z coords of the atoms in that unit cell
# Returns a 1-D numpy array of nested lists
def get_neighbor_list_addresses(unit_cell_3d_array, center_cell_index_list, neighbor_distance=0,
                                neighbor_tolerance=1.05):
    coords_tuple = []
    index_list = []
    x, y, z = center_cell_index_list
    # print(center_cell_index_list)
    center_2d_array = unit_cell_3d_array[x, y, z]
    if center_2d_array.size == 0:
        return [[]]
    x_range, y_range, z_range = [], [], []
    if x > 0:
        x_range = [x - 1, x]
    else:
        x_range = [x]
    if x + 1 < unit_cell_3d_array.shape[0]:
        x_range.append(x + 1)
    if y > 0:
        y_range = [y - 1, y]
    else:
        y_range = [y]
    if y + 1 < unit_cell_3d_array.shape[1]:
        y_range.append(y + 1)
    if z > 0:
        z_range = [z - 1, z]
    else:
        z_range = [z]
    if z + 1 < unit_cell_3d_array.shape[2]:
        z_range.append(z + 1)
    address_nested = [[] for i in list(range(center_2d_array.shape[0]))]
    if center_2d_array.size == 0:
        return [[]]
    # print('Ranges:')
    # print(x_range, y_range, z_range)

    for x_cell, y_cell, z_cell in itertools.product(x_range, y_range, z_range):
        neighbor_2d_array = unit_cell_3d_array[x_cell, y_cell, z_cell]
        if neighbor_2d_array.size == 0:
            # print(neighbor_2d_array.shape)
            continue
        # print(x)
        # neighbor_2d_array = np.array([[0,0,0], [10,4,3]])
        # center_2d_array = np.array([[10,0,0], [0, 20, 0], [8,8,3]])
        neighbor_broadcast_T = np.broadcast_to(neighbor_2d_array,
                                               (center_2d_array.shape[0], neighbor_2d_array.shape[0], 3))
        neighbor_broadcast = np.swapaxes(neighbor_broadcast_T, 0, 1)

        center_broadcast = np.broadcast_to(center_2d_array, (neighbor_2d_array.shape[0], center_2d_array.shape[0], 3))
        '''
        print(center_broadcast.shape)
        print(neighbor_broadcast.shape)
        '''
        # print(center_broadcast)
        # print('\nNeighbor broadcast: {}\n'.format(neighbor_broadcast))
        # print(neighbor_2d_array)
        # center_broadcast = center_2d_array
        # neighbor_broadcast = neighbor_2d_array
        # print('\n{}\n'.format(center_broadcast-neighbor_broadcast))
        distance_2d_matrix = np.sqrt(np.sum(np.square(center_broadcast - neighbor_broadcast), axis=2)).T
        # print(neighbor_2d_array.shape)
        # print(center_2d_array.shape)
        # print(distance_2d_matrix)
        if neighbor_distance == 0:
            neighbor_distance = np.min(distance_2d_matrix[distance_2d_matrix > 0])

        for center_index in list(range(center_2d_array.shape[0])):
            if neighbor_2d_array.shape[0] > 1:
                for target_index in list(range(neighbor_2d_array.shape[0])):
                    dist = distance_2d_matrix[center_index, target_index]
                    if dist > neighbor_distance / 2 and dist < neighbor_tolerance * neighbor_distance:
                        address_nested[center_index].append([x_cell, y_cell, z_cell, target_index])
            else:
                dist = distance_2d_matrix[center_index]
                if dist > neighbor_distance / 2 and dist < neighbor_tolerance * neighbor_distance:
                    address_nested[center_index].append([x_cell, y_cell, z_cell, target_index])
                    # neighbor_indices = np.where(nearest_mask[center_index, :])
                    # print(neighbor_indices)
    # print(address_nested)
    address_array = np.array(address_nested)
    return address_array


'''
       else:
            neighbor_indices = np.where(nearest_mask)
            if neighbor_indices:
                #print(neighbor_indices)
                for target_index in neighbor_indices:
                    address_nested[0].append([x_cell, y_cell, z_cell, target_index])
'''


def get_unit_cell_coords(unit_cell_length, crystal_struct="fcc", origin_adj=0):
    # surface_thickness_integer = int(np.round((unit_cell_length / np.sqrt(2))))
    # print(surface_thickness_integer)
    # cell_length_x, cell_length_y, cell_length_z = unit_cell_length, unit_cell_length, unit_cell_length
    fcc_fractional_coords_array = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    unit_cell_float = fcc_fractional_coords_array * unit_cell_length
    unit_array = unit_cell_float - origin_adj * unit_cell_length
    return unit_array


def sphere_check(coords, sphere_radius, surface_thickness=0):
    # if abs(coords[0]) > sphere_radius or abs(coords[1]) > sphere_radius or abs(coords[2]) > sphere_radius:
    #    return False
    distance = np.sum(np.square(coords))
    if distance - sphere_radius ** 2 < -0.001:
        if distance - (sphere_radius - surface_thickness) ** 2 > 0.001:
            return "surface"
        else:
            return "core"


def nanoparticle_generator(lattice, geometry_check, *args, element='Au'):
    surface_thickness_integer = int(np.round((lattice.unit_cell_lengths[0] / 2)))
    for x_cell, y_cell, z_cell in itertools.product(*lattice.unit_cell_indices):
        new_unit_cell = lattice.create_unit_cell(np.array([x_cell, y_cell, z_cell]))
        for x_cell_coords, y_cell_coords, z_cell_coords in lattice.unit_cell_coords:
            x_coords = x_cell_coords + x_cell * lattice.unit_cell_lengths[0]
            y_coords = y_cell_coords + y_cell * lattice.unit_cell_lengths[1]
            z_coords = z_cell_coords + z_cell * lattice.unit_cell_lengths[2]
            new_coords = [x_coords, y_coords, z_coords]
            check = geometry_check(new_coords, *args)  # , surface_thickness=surface_thickness_integer)
            if check:
                # Create a new atom in SuperLattice class if it's part of the nanoparticle.
                new_atom = lattice.create_atom(element, new_coords, new_unit_cell)
                if check == "surface":
                    lattice.surface_atoms.append(new_atom)
    return lattice


def nanoparticle_builder(radius_ang, element='Au', unit_cell_length_int=4070):
    origin_offset = 0
    radius_integer = int(np.round(1000 * radius_ang))  # - cell_length_x / (2 * np.sqrt(2))))
    radius_in_unit_cells = int(np.round((math.ceil(float(radius_integer) / float(unit_cell_length_int)))) + 1)
    unit_cell_range = list(range(-radius_in_unit_cells, radius_in_unit_cells + 1))
    unit_array = get_unit_cell_coords(unit_cell_length_int, origin_adj=origin_offset)
    vdw_radius_integer = int(np.round(unit_cell_length_int / (2. * np.sqrt(2)), decimals=0))
    lattice = NanoParticle.NanoparticleSuperLattice([unit_cell_range, unit_cell_range, unit_cell_range],
                                                    unit_cell_coords=unit_array,
                                                    unit_cell_lengths=[unit_cell_length_int] * 3,
                                                    neighbor_distance=2 * vdw_radius_integer)
    lattice = nanoparticle_generator(lattice, sphere_check, radius_integer)
    lattice.get_all_neighbors()
    # lattice.get_surface_atoms()
    print('There are {} surface atoms.'.format(len(lattice.surface_atoms)))
    print('There are {} total atoms.'.format(len(lattice.atoms)))
    lattice.assign_atom_names()
    # Choose ligand_sites.

    # Bond all gold atoms.
    lattice.bond_all_atoms()
    return lattice


def make_fcc_111(surface_dimensions, atomic_radius):
    coords_list = []
    num_atoms = [10, 6, 2]
    lattice_dim = [list(range(total)) for total in num_atoms]
    unit_cell = [[2 * atomic_radius, 0, 0], [atomic_radius, 3 ** 0.5 * atomic_radius, 0],
                 [0, 2 * atomic_radius * 3 ** 0.5, 0]]
    z_operator = np.array([atomic_radius, 3 ** 0.5 * atomic_radius / 3, -atomic_radius * 2 * 2 ** 0.5 / 3 ** 0.5])
    # print(lattice_dim)
    for xi, yi, zi in itertools.product(*lattice_dim):
        # print(xi,yi,zi)
        for xp, yp, zp in unit_cell:
            coords = np.array(
                [xp + xi * 2 * atomic_radius, yp + yi * 2 * 3 ** 0.5 * atomic_radius, zi]) + zi * z_operator - np.array(
                [0, 0, atomic_radius])
            tuple_coords = coords.tolist()
            # print(xp,yp,zp)
            # print(coords)
            if tuple_coords not in coords_list:
                coords_list.append(tuple_coords)
    return coords_list


def nanorod_xyz_builder(length_int, radius_int, direction=False, cap_angle=math.pi / 4, cap_length_int=-1):
    # Important relations.
    # Side length = Width / 2.414 = Radius / 1.207
    # Or by inradius: 4 = (1 + sqrt(2))/2 * Side length
    # Boundary of cross-section in first quadrant:
    # y = -x + b
    # b = 2*sqrt(2)/(sqrt(2)+1) * radius => Derived by substituting coords for first corner of octagon into slope-intercept form.
    def nanorod_xsection_check(x_coord, y_coord, length_center_int, radius_center_int, body_slope_height_int=0):
        if abs(x_coord) < radius_center_int or abs(y_coord) < radius_center_int:
            return False
        if abs(x_coord) + abs(y_coord) - body_slope_height_int < 0.001:
            return True

    def nanorod_check(coords, length_center_int, radius_center_int, cap_length_int, cap_angle=m.pi / 4.,
                      body_slope_height_int=0):

        def nanorod_body_check(coords, radius_center_int, body_slope_height_int):
            # print(abs(coords[0]) + abs(coords[1]) - body_slope_height_int < 0.001)
            correct_100_plane = vdw_radius_integer
            return abs(coords[0]) < radius_center_int and abs(coords[1]) < radius_center_int and abs(coords[0]) + abs(
                coords[1]) - body_slope_height_int - correct_100_plane < 0.001

        def nanorod_cap_check(coords, half_body_length, radius_center_int, cap_length_int, cap_angle,
                              slope_factor=(1. + m.tan(m.pi / 8.))):
            diagonal_correction = 0
            # slope_factor = 1
            cap_slope_height = int(
                np.round(radius_center_int - np.tan(cap_angle) * (abs(coords[2]) - half_body_length)))
            x_check = abs(coords[0]) - cap_slope_height < 0.001
            y_check = abs(coords[1]) - cap_slope_height < 0.001
            diag_check = abs(coords[0]) + abs(coords[1]) - slope_factor * cap_slope_height + diagonal_correction < 0.001
            # print(coords, x_check, y_check, diag_check)
            return x_check and y_check and diag_check

        if (abs(coords[2]) > length_center_int // 2):
            return False

        slope_factor = (1 + m.tan(m.pi / 8))  # Valid for 45 degree slope.
        half_body_length = length_center_int // 2 - cap_length_int

        # if (body_slope_height_int == 0):
        #    body_slope_height_int = int(np.round(slope_factor * radius_center_int, decimals=0)) // 2
        oct_slope_intercept = radius_int * 2 * np.sqrt(2) / (np.sqrt(2) + 1)
        # body_slope_height_int = int(np.round(slope_factor * radius_center_int, decimals=0)) // 2
        if (abs(coords[2]) - half_body_length < 0.001):
            # print('Part of main nanorod')
            return nanorod_body_check(coords, radius_center_int, oct_slope_intercept)
        else:
            # print(coords[2])
            return nanorod_cap_check(coords, half_body_length, radius_center_int, cap_length_int, cap_angle,
                                     slope_factor=slope_factor)

    unit_cell_length_int = 4070
    element = 'Au'
    surface_thickness_integer = int(np.round((unit_cell_length_int / np.sqrt(2))))
    # print(surface_thickness_integer)
    cell_length_x, cell_length_y, cell_length_z = unit_cell_length_int, unit_cell_length_int, unit_cell_length_int
    cell_lengths = [unit_cell_length_int, unit_cell_length_int, unit_cell_length_int]
    fcc_fractional_coords_array = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    unit_cell_float = fcc_fractional_coords_array * unit_cell_length_int
    unit_cell_origin_fractional_adjust = 0.5
    unit_array = unit_cell_float - unit_cell_origin_fractional_adjust * unit_cell_length_int
    ### Set radius_center_integer of cluster. Subtract VdW radius_center_integer from radius_center_integer.
    vdw_radius_integer = int(np.round(cell_length_x / (2. * np.sqrt(2)), decimals=0))
    radius_center_integer = int(np.round(radius_int - vdw_radius_integer, decimals=0))
    length_center_integer = int(np.round(length_int - vdw_radius_integer, decimals=0))
    if cap_length_int == -1:
        cap_length_int = radius_center_integer // 2
        print('Cap length in milliAngstrom is {}'.format(cap_length_int))
    # Set how many times to extend lattice in the x_fractional,y_fractional,z_fractional direction;
    radius_in_unit_cells = int(np.ceil((radius_center_integer / cell_length_x)))
    length_in_unit_cells = int(2 * np.ceil(length_center_integer / (2 * cell_length_z)))
    if length_in_unit_cells // 2 * 2 < length_in_unit_cells:
        length_in_unit_cells += 1
    # x_unit_cells = range(-radius_in_unit_cells, radius_in_unit_cells)
    # y_unit_cells = range(-radius_in_unit_cells, radius_in_unit_cells)
    # z_unit_cells = range(-length_in_unit_cells, length_in_unit_cells)
    x_cells_range = list(range(-radius_in_unit_cells - 1, radius_in_unit_cells + 2))
    y_cells_range = list(range(-radius_in_unit_cells - 1, radius_in_unit_cells + 2))
    z_cells_range = list(range(-length_in_unit_cells // 2, 1 + length_in_unit_cells // 2))
    # Set dimensions of nanoparticle

    # ## Set dimensions of FCC unit cell
    # if direction == "110":
    #    maxx, maxy, maxz = 4.07 * np.sqrt(2), 4.07 * np.sqrt(2), 4.07
    #    unit_array = np.array(
    #        [[0., 0., 0.], [1. / np.sqrt(2.), 0., 0.], [0., 1. / np.sqrt(2.), 0.], [0., 2.035, 2.035]])
    lattice = NanoParticle.NanoparticleSuperLattice([x_cells_range, y_cells_range, z_cells_range],
                                                    unit_cell_coords=unit_array, unit_cell_lengths=cell_lengths,
                                                    neighbor_distance=2 * vdw_radius_integer)
    lattice = nanoparticle_generator(lattice, nanorod_check, length_center_integer, radius_center_integer,
                                     cap_length_int)

    """
    lattice = nanoparticle_generator(lattice, nanorod_check())
    for x_cell in x_cells_range:
        #x_cell_translation = x_cell - radius_in_unit_cells
        for y_cell in y_cells_range:
            #y_cell_translation = y_cell - radius_in_unit_cells
            for z_cell in z_cells_range:
                #z_cell_translation = z_cell - length_in_unit_cells // 2
                # Create new unit cell.
                new_unit_cell = lattice.create_unit_cell(np.array([x_cell, y_cell, z_cell]))
                for x_cell_coords, y_cell_coords, z_cell_coords in unit_array:
                    x_coords = x_cell_coords + x_cell * cell_length_x
                    y_coords = y_cell_coords + y_cell * cell_length_y
                    z_coords = z_cell_coords + z_cell * cell_length_z
                    new_coords = [x_coords, y_coords, z_coords]
                    if nanorod_check(new_coords, length_center_integer, radius_center_integer, cap_length_int):
                        # Create a new atom in SuperLattice class if it's part of the nanoparticle.
                        lattice.create_atom(element, new_coords, new_unit_cell)
    """
    print('There are {} atoms in the nanorod.\n\n'.format(len(lattice.atoms)))
    # Determine which atoms are on surface.
    lattice.get_all_neighbors()
    lattice.get_surface_atoms()
    print('There are {} surface atoms.'.format(len(lattice.surface_atoms)))
    lattice.assign_atom_names()
    # Choose ligand_sites.

    # Bond all gold atoms.
    # lattice.bond_all_atoms()
    return lattice

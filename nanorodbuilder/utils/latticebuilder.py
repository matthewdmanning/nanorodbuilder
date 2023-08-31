import math

import geometry_tools
import script_writing_tools


def read_xyz_file(file_name, include_hydrogens=True):
    coords_tuple = []
    atom_type_list = []
    xyz_file = open(file_name, 'r')
    min_max_tuple = [[0, 0, 0], [0, 0, 0]]
    for row_num, row in enumerate(xyz_file):
        if row_num < 2:
            continue
        atom_type, xs, ys, zs = row.split()
        x, y, z = float(xs), float(ys), float(zs)
        if row_num == 2:
            min_max_tuple = [[x, x], [y, y], [z, z]]
        if not include_hydrogens and atom_type in ['N', 'O', 'H']:
            continue
        coords_tuple.append([x, y, z])
        for num, coord in enumerate([x, y, z]):
            if coord - min_max_tuple[num][0] < 0.0001:
                min_max_tuple[num][0] = math.floor(coord)
            if coord - min_max_tuple[num][1] > 0.0001:
                min_max_tuple[num][1] = math.ceil(coord)
        atom_type_list.append(atom_type)
    print(min_max_tuple)
    return coords_tuple, atom_type_list, min_max_tuple


pdb_coords_tuple, atom_type_list, min_max_tuple = read_xyz_file('/home/mdmannin/Downloads/1kx5_remainder_trimmed.xyz',
                                                                include_hydrogens=False)
distance_cutoff = 1.
lattice_params = [5., 5., 5.]
pdb_size = [math.ceil(mm[1] - mm[0]) for mm in min_max_tuple]
supercell_size = [math.ceil(length / param) for length, param in zip(pdb_size, lattice_params)]
sc_lattice = geometry_tools.lattice_builder_from_ranges(lattice_params, min_max_tuple)
print(sc_lattice)
# print(bcc_lattice)

cg_coords_tuple = []
for coords in sc_lattice:
    match_found = False
    for pdb_coords in pdb_coords_tuple:
        # print(pdb_coords)
        if match_found:
            continue
        if abs(coords[0] - pdb_coords[0]) > distance_cutoff or abs(coords[1] - pdb_coords[1]) > distance_cutoff or abs(
                coords[2] - pdb_coords[2]) > distance_cutoff:
            continue
        # distance = geometry_tools.get_taxicab_distance(coords, pdb_coords)
        # if distance - distance_cutoff < 0.0001:
        if geometry_tools.calculate_distance(coords, pdb_coords) - distance_cutoff < 0.0001:
            cg_coords_tuple.append(coords)
            match_found = True
            print(coords)
print(cg_coords_tuple)
bcc_lattice = geometry_tools.convert_from_simple_cubic(cg_coords_tuple, 'bcc', lattice_params)

script_writing_tools.xyz_writer('/home/mdmannin/Desktop/cg_1kx5_remainder_3.xyz', cg_coords_tuple, element_list='ZZ')

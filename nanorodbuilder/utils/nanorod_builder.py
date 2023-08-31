#!/usr/bin/env python
import itertools
import logging
import math
import os
import random

import NanoParticle
import System
import numpy as np
from numpy import int

import script_writing_object
from nanorodbuilder.utils import geometry_tools

nearest_neighbor_distance = float(4070 / np.sqrt(2) * 1.1)
print(nearest_neighbor_distance)
solvent_buffer = int(10)
force_charge_ratio = False
# ## Outline
# 1. Fetch dimensions and files (gold mol2 files (surf and core) and ligand mol2).
# 2. Check for overlap and delete non-ligand gold.
# 3. Open shell file and import parameters and FFs. (Load frcmod files from RESP output)
# 3. Throw ligands, bind to gold, and 'impose bond length'.
# 4. Bind gold atoms.
#    Save mol2? or prmtop? file and use antechamber to rename/combine gold residues?
# 5. Solvate/salt and build topology

# Merging nanorodbatch.py
# Leap Commands after confirming part of nanorod.
# system = createUnit system?
# AU = createResidue
# createAtom name ligand_type(Au) charge
# ??? set ATOM element Au
# set ATOM position { x y z }

filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/'
parmpath = '{}ParmFiles/'.format(filepath)
# print '\n\n INPUT FILES COMING FROM RENAMED MOL2. MAKE SURE THEY HAVE BEEN CREATED.\n\n'

print('\n\n Nanorod Builder')


# ===========================================================================
def overlap_check(atoms, goldid, protected, goldsulflig, residuelist, namelist, residuedelete):
    atomsdelete = []
    ausurvivors = []
    atomsize = atoms.shape[0]
    for counter, atom in enumerate(atoms):
        if counter in atomsdelete or residuelist[int(atom[0])] in residuedelete: continue
        for i in range(counter + 1, atomsize - 1):
            checkatom = atoms[i]
            if i in atomsdelete or residuelist[int(checkatom[0])] in residuedelete: continue
            if residuelist[int(atom[0])] == residuelist[int(checkatom[0])] and int(atom[0]) not in goldid: continue
            if math.sqrt((atom[1] - checkatom[1]) ** 2 + (atom[2] - checkatom[2]) ** 2 + (
                    atom[3] - checkatom[3]) ** 2) - 1.0 < 0.001 or (
                    atom[1] == checkatom[1] and atom[2] == checkatom[2] and atom[3] == checkatom[3]):
                if int(atom[0]) in goldid and int(checkatom[0]) not in goldid:
                    # print 'Overlapping ligand containing Atom #{}.'.format(atom[0])
                    residuedelete.append(residuelist[int(checkatom[0])])
                    print('Overlapping ligand found. Deleting ligand residue: nanorod.{}.\n'.format(
                        residuelist[int(checkatom[0])]))
                elif int(atom[0]) not in goldid and int(checkatom[0]) in goldid:
                    residuedelete.append(residuelist[int(atom[0])])
                    print('Overlapping ligand found. Deleting ligand residue: nanorod.{}.\n'.format(
                        residuelist[int(atom[0])]))
                elif int(atom[0]) in goldid and int(checkatom[0]) in goldid:
                    # ##Delete overlapping atom if not bound to ligand.
                    if int(checkatom[0]) not in protected:
                        atomsdelete.append(i)
                        ausurvivors.append(counter)
                        # print atom
                        # print checkatom
                        # print 'Gold overlap deleted ({})'.format(namelist[i])
                    # ##Atom being check is overlapping with a ligand-bound atom.
                    elif int(atom[0]) not in protected:
                        atomsdelete.append(counter)
                        ausurvivors.append(i)
                        # print 'Gold overlap deleted ({})'.format(namelist[counter])
                        break
                    # ## Both atoms are bound to ligands. -> Not good. Choose one gold/ligand and delete.
                    else:
                        row = goldsulflig[:, 0].index(int(checkatom[0]))
                        atomsdelete.append(i)
                        ausurvivors.append(counter)
                        residuedelete.append(goldsulflig[row, 2])
    return 1, atomsdelete, residuedelete, ausurvivors


def nano_cylinder_surf(length, radius):
    cylinder_area = 2. * np.pi * float(radius) * float(length)
    cap_area = 4. * np.pi * (float(radius)) ** 2
    total_area = cylinder_area + cap_area
    return total_area


def calc_ligand_num(length_int, radius_int, in_density_nm2, ratio, lig_length=12):
    # Correct percentage form of ratio.
    if ratio - 1.01 > 0.0001:
        ratio = ratio / 100.
    # Takes input as milli-Angstrom (ie. basis vector is 4070)
    # Density taken directly from file as 1/nm**2.
    # in_density_nm2 = in_density#_nm2 * 10 ** (-8)
    length_nm = length_int / 10000.
    radius_nm = radius_int / 10000.
    total_area = nano_cylinder_surf(length_nm, radius_nm)
    ligand_num = int(total_area * float(in_density_nm2))  # Necessary to convert from sq Angstroms to sq nanometers.
    density_actual = float(ligand_num) / total_area
    # print('Surface density for {} ligands is {}.'.format(ligand_num, density_actual))
    if ratio == 100:
        charged_num = ligand_num
    else:
        charged_num = int(np.round(ligand_num * ratio, decimals=0))
    uncharged_num = ligand_num - charged_num
    ### Change for different ligands.
    ligand_surface_area = nano_cylinder_surf(length_nm, (radius_nm + lig_length))
    charge_density = charged_num / ligand_surface_area
    return charged_num, ligand_num


def calc_end_coords(length, radius, coords_list, lig_radius):
    length, radius, lig_radius = float(length), float(radius), float(lig_radius)
    end_group_coords_list = []
    for coords in coords_list:
        if abs(coords[2]) - float(length) / 2 - float(radius) < 0.0001:
            newx = coords[0] * (radius + lig_radius) / lig_radius
            newy = coords[1] * (radius + lig_radius) / lig_radius
            end_group_coords_list.append([newx, newy, coords[2]])
        else:
            center = coords[2] - (length / 2 - radius) * abs(coords[2]) / coords[2]
            newz = center * (radius + lig_radius) / lig_radius
            newy = coords[1] * (radius + lig_radius) / lig_radius
            newx = coords[0] * (radius + lig_radius) / lig_radius
            end_group_coords_list.append([newx, newy, newz])
    return end_group_coords_list


def calc_potential_coords(length, radius, coords_list, lig_radius):
    # Set both charges and permitivity to 1 to simplify calcs.
    potentials = np.zeros(len(coords_list))
    end_coords = calc_end_coords(length, radius, lig_radius, coords_list)
    for (ind1, coords1), (ind2, coords2) in itertools.combinations(enumerate(end_coords), 2):
        potentials[ind1, ind2] = 1. / geometry_tools.calculate_distance(coords1, coords2)
    total_pot = sum(potentials)
    return total_pot


def calc_potential_distances(charged_index, dist_pairs):
    potential = 0
    charged_index.sort()
    for i, j in itertools.combinations(charged_index, 2):
        next_pot = abs(1. / dist_pairs[i, j])
        potential = potential + next_pot
    return potential


def hard_core_ligand_distributor(length, radius, density, ratio, bound_array):
    charged_num = int(len(bound_array) * ratio / 100)
    if charged_num == len(bound_array):
        return bound_array, []
    elif charged_num == 0:
        return [], bound_array
    charged_array, uncharged_array = [], []
    dist_pairs = geometry_tools.distance_pairs(bound_array)
    # while len(charged_array) < charged_num:


def ligand_distributor(length, radius, density, ratio, bound_array):
    # Random placement of charged vs uncharged ligands. Selected by lowest electrostatic potential of corona.
    # Shuffling is only done on list of indices, not on the original coordinates or any of its pointers.
    mutate_num = 10000
    charged_num = int(len(bound_array) * ratio / 100)
    if charged_num == len(bound_array):
        return bound_array, []
    elif charged_num == 0:
        return [], bound_array
    charged_array, uncharged_array = [], []
    dist_pairs = geometry_tools.distance_pairs(bound_array)

    # Generate random initial placement of charged ligands.
    rand_index = list(range(len(bound_array)))
    random.shuffle(rand_index)
    charged_index = rand_index[:charged_num]
    uncharged_index = rand_index[charged_num:]
    # Loop that mutates ligand sites randomly, calculates the electrostatic potential of corona, and accepts or rejects the new placement.
    for i in range(mutate_num):
        random.shuffle(rand_index)
        charged_copy = rand_index[:charged_num]
        uncharged_copy = rand_index[charged_num:]
        new_pot = calc_potential_distances(charged_copy, dist_pairs)
        old_pot = calc_potential_distances(charged_index, dist_pairs)
        if new_pot - old_pot < 0.000001:
            # print("Lower configuration found.")
            charged_index = charged_copy[:]
            uncharged_index = uncharged_copy[:]
            print(new_pot)
    for ind in charged_index: charged_array.append(bound_array[ind])
    for ind in uncharged_index: uncharged_array.append(bound_array[ind])

    return charged_array, uncharged_array


def parameters():
    lengthlist, radiuslist, densitylist, ligand1_list, ligand2_list, sulfindex1_list, sulfindex2list, ratiolist = [], [], [], [], [], [], [], []
    batch_file_name = '../data/batch_input.txt'
    batch_file = open(batch_file_name, 'r')
    ligand_name_info_dict = {}
    for line in batch_file:
        if line.find("LENGTH=") != -1:
            lengthbite = line[line.index('=') + 1:].split(',')
            for i in lengthbite: lengthlist.append(int(i.rstrip()))
        elif line.find("RADIUS=") != -1:
            radiusbite = line[line.index('=') + 1:].split(',')
            for i in radiusbite: radiuslist.append(int(i.rstrip()))
        elif line.find("DENSITY=") != -1:
            densitybite = line[line.index('=') + 1:].split(',')
            for i in densitybite: densitylist.append(int(i.rstrip()))
        elif line.find("RATIO=") != -1:
            ratiobite = line[line.index('=') + 1:].split(',')
            for i in ratiobite: ratiolist.append(100 * float(i.rstrip()))
        elif line.find("LIGAND1=") != -1:
            ligand1bite = line[line.index('=') + 1:].split(',')
            for i in ligand1bite:
                ligand_name = i.rstrip()
                ligand1_list.append(ligand_name)
                if ligand_name not in ligand_name_info_dict.keys():
                    new_ligand = NanoParticle.LigandInfo(ligand_name)
                    ligand_name_info_dict[ligand_name] = new_ligand
        elif line.find("LIGAND2=") != -1:
            ligand2bite = line[line.index('=') + 1:].split(',')
            for i in ligand2bite:
                ligand_name = i.rstrip()
                ligand2_list.append(ligand_name)
                if ligand_name not in ligand_name_info_dict.keys():
                    new_ligand = NanoParticle.LigandInfo(ligand_name)
                    ligand_name_info_dict[ligand_name] = new_ligand
    batch_file.close()
    # Get ligand info.
    ligand_info_filename = "../data/ligand_info.txt"
    with open(ligand_info_filename, 'r') as batch_file:
        for line in batch_file:
            for ligand_name, ligand_info in ligand_name_info_dict.items():
                if line.find(ligand_name) != -1:
                    ligand_info_string = line[line.index(ligand_name) + len(ligand_name) + 1:].split()
                    sulfur_index = ligand_info_string[0]
                    ligand_info.attach_index = sulfur_index
    return lengthlist, radiuslist, densitylist, ligand1_list, ligand2_list, ligand_name_info_dict, ratiolist


def master_file(in_names, leap_names):
    # ## Create master file for tLeap execution.
    grep_list = ['max', 'improper', 'missing', 'WARNING -A 1']
    filepath = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/'
    parmpath = '{}ParmFiles/'.format(filepath)
    shellpath = '{}LeapScripts/'.format(filepath)
    master_shell_file = open('{}masterscript.sh'.format(shellpath), 'w+')
    master_shell_file.write("#/bin/bash\n\n")
    for in_file, in_leap in zip(in_names, leap_names):
        master_shell_file.write('rm {}.log\n'.format(in_file))
        master_shell_file.write('tleap -f {} -log {}.log\n'.format(in_leap, in_file))
        master_shell_file.write('sleep 5\n')
    for in_file, search_string in itertools.product(in_names, grep_list):
        master_shell_file.write('grep {} {}.log\n'.format(search_string, in_file))
    master_shell_file.close()
    print('Master shell script created: {}masterscript.sh.\n'.format(shellpath))
    return


def main():
    log_path = os.getcwd()
    logging.basicConfig(filename='{}leaptotalbuilder.log'.format(log_path), filemode='w', level=logging.INFO)
    # ## Main Code
    script_path = '{}LeapScripts/'.format(filepath)
    ### Get nanorod parameters from batch_input.txt file.
    in_names, leap_names = [], []
    systems_list = []
    for nanorod_ang_length, nanorod_ang_radius in itertools.product(length_list, radius_list):
        parmpath = '{}ParmFiles/'.format(filepath)
        # Get lists of coordinates of gold atoms.
        nanorod_int_length = int(np.round(1000 * nanorod_ang_length, decimals=0))
        nanorod_int_radius = int(np.round(1000 * nanorod_ang_radius, decimals=0))
        print('Nanorod length as milliAngstrom: {}'.format(nanorod_int_length))
        print('Nanorod width as milliAngstrom: {}'.format(2 * nanorod_int_radius))
        # Create instance of NanoparticleSuperLattice.
        nanorod_lattice = geometry_tools.nanorod_xyz_builder(nanorod_int_length, nanorod_int_radius)
        print(len(nanorod_lattice.atoms))
        # nanorod_lattice.bond_all_atoms()
        script_writing_object.xyz_writer(nanorod_lattice,
                                         '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/NanorodXYZFiles/{}x{}-nanorod.xyz'.format(
                                             nanorod_ang_length, nanorod_ang_radius))
        # geometry_tools.plot_nanoparticle_with_ligands([atom.coords for atom in nanorod_lattice.atoms], [], [], nearest_neighbor_distance // 2)
        continue
        print(ratio_list)
        for ligand_density in density_list:
            blank, total_ligand_num = calc_ligand_num(nanorod_int_length, nanorod_int_radius, ligand_density, 1.0)
            generic_nanorod = NanoParticle.GenericParticle(nanorod_lattice)
            generic_nanorod.choose_unspecified_sites(total_ligand_num)
            for ligand_ratio in ratio_list:
                charged_num = int(np.round(total_ligand_num * ligand_ratio / 100.))
                uncharged_num = total_ligand_num - charged_num
                # Create placeholder ligand types.
                if ligand_ratio == 100:
                    num_ligand_types = 1  # Change to variable based on input file. Not a priority right now.
                    ligand_type_charges = [1]
                    ligand_num_list = [total_ligand_num]
                else:
                    num_ligand_types = 2
                    ligand_type_charges = [1, 0]
                    ligand_num_list = [charged_num, uncharged_num]
                # Create new particle for assigning ligand types.
                generic_nanorod.clear_placeholder_sites()
                # print('Ligand numbers: {}'.format(ligand_num_list))
                generic_nanorod.choose_placeholder_sites(ligand_num_list=ligand_num_list)
                print(ligand1_list, ligand2_list)
                for ligand1_name, ligand2_name in zip(ligand1_list, ligand2_list):
                    if ligand_ratio < 100:
                        combined_ligands = '{}-{}'.format(ligand1_name, ligand2_name)
                        ligand_names_list = [ligand1_name, ligand2_name]
                    else:
                        combined_ligands = ligand1_name
                        ligand_names_list = [ligand1_name]
                    # Create Nanoparticle object that specifies ligand types.
                    complete_nanorod = NanoParticle.Nanoparticle(generic_nanorod)
                    ligand_info_list = [ligand_name_info_dict[name] for name in ligand_names_list]
                    ligand_info_num_dict = {ligand_name_info_dict[ligand1_name]: charged_num,
                                            ligand_name_info_dict[ligand2_name]: uncharged_num}
                    complete_nanorod.create_ligand_sites(ligand_info_num_dict=ligand_info_num_dict)
                    mol2_filename = (
                        '{}x{}-{}r-{}d-{}'.format(nanorod_ang_length, nanorod_ang_radius, int(ligand_ratio),
                                                  ligand_density, combined_ligands))
                    leap_filename = '{}{}_leap.in'.format(script_path, mol2_filename)
                    new_system = System.LeapSystem(leap_filname=leap_filename, amber_unit_name=mol2_filename,
                                                   system_name=mol2_filename)

                    systems_list.append(new_system)
                    print('Saving script to {}.'.format(leap_filename))
                    script_writing_object.write_leap_input(new_system, complete_nanorod, prmtop_name=mol2_filename)

    script_writing_object.master_file(systems_list, script_path='/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/')


if __name__ == "__main__":
    main()

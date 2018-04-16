import collections
import copy
import itertools
import math
import random
import string
import geometry_tools

import numpy as np

file_path = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/'
ligand_file_path = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/'
default_parm_path = "{}ParmFiles/".format(file_path)


def get_next_unused_letter(letter, list_of_used_letters, case_sensitive=True, random_letter=False, case=False):
    unused_letters = [chara for chara in string.ascii_letters if chara not in list_of_used_letters]
    if case == 'upper':
        available_letters = [chara for chara in unused_letters if chara.isupper()]
    elif case == 'lower':
        available_letters = [chara for chara in unused_letters if chara.islower()]
    if random_letter:
        return random.choice(available_letters)
    else:
        while ord(available_letters[0]) < ord(letter):
            available_letters.append(available_letters.pop(0))
        return available_letters[0]
        # Scraps if other case of a used letter shouldn't be considered.
        # if letter.swapcase() not in list_of_used_letters: return letter.swapcase()
        # else:
        #    if not random_letter: letter = chr(ord(letter) + letter_assignment)
        #    else: letter = chr()


def generate_atom_names_from_elements():
    print('This function is not ready.')
    raise Exception
    max_char_in_name = 4
    if len(atom_letter_prefix_list) == 0:
        atom_letter_prefix_list = [symbol[0] for symbol in element_symbol_list]
        prefix_counter = collections.Counter(copy.deepcopy(atom_letter_prefix_list))
        for prefix, prefix_count in list(prefix_counter):
            max_prefix_count = 10 ** (max_char_in_name - len(prefix))
            if prefix_count > max_prefix_count:
                indices = [index for index, unique_prefix in atom_letter_prefix_list if unique_prefix is prefix]
                new_prefix = prefix.swapcase()
                while prefix_count > max_prefix_count:
                    if prefix.swapcase() not in list(prefix_counter):
                        new_prefix = prefix.swapcase()
                for index in indices[max_prefix_count:min(2 * max_prefix_count, len(indices))]:
                    atom_letter_prefix_list[index] = atom_letter_prefix_list[index].swapcase()
                new_prefix = get_next_unused_letter(new_prefix, used_prefixed)


def check_file_exists(file_path):
    return True
    import os
    return os.path.isfile(file_path)


def check_amber_file_exists(file_name):
    return True
    import os
    try:
        file_path = os.environ['AMBERHOME'] + file_name
        return os.path.isfile(file_path)
    except:
        print('$AMBERHOME environment variable not set. Try launching this script from the command line.')
        return False


def check_pdbseq_validity(seq, amber_recognized_residues=None):
    if not amber_recognized_residues:
        return
    residue_string = seq.strip('{}')
    for residue in residue_string.split():
        if residue in amber_recognized_residues:
            continue
        else:
            print('Warning: Unrecognized residue {} in loadpdbUsingSeq command. Check sequence.'.format(residue))
            return False
    return True


def load_force_fields(leap_script, force_field_list, frcmod_list):
    for forcefield in force_field_list:
        force_field_name = 'leaprc.{}'.format(forcefield)
        if check_amber_file_exists(force_field_name):
            leap_script.write('source leaprc.{}\n'.format(forcefield))
        else:
            print('Force field file, {}, cannot be found in $AMBERHOME directory. Check spelling or try launching this '
                  'script from the command line.'.format(force_field_name))
    for frcmod in frcmod_list:
        # if check_amber_file_exists(frcmod):
        leap_script.write('loadamberparams {}\n'.format(frcmod))
        # else:
        #    print(
        #        '{} cannot be found in $AMBERHOME directory. Check spelling or try launching this script from the command'
        #        ' line.'.format(frcmod))
    return leap_script


def load_molecule_auto(leap_script, molecule_unit_alias, file_name, unit=False, seq=False, alignaxes=False):
    check_amber_file_exists(file_name)
    if file_name.endswith('.mol2'):
        load_command = 'loadmol2'
    elif seq:
        load_command = 'loadPdbUsingSeq'
        seq = check_pdbseq_validity(seq)
        file_name = ' '.join(file_name)
    elif file_name.endswith('.pdb') and not seq:
        load_command = 'loadpdb'
    else:
        print('File ligand_type not recognized. Check ending for ".pdb" or ".mol2".')
        return leap_script
    leap_script.write('{} = {} {}\n'.format(molecule_unit_alias, load_command, file_name))
    if alignaxes:
        leap_script.write('alignaxes {}\n'.format(molecule_unit_alias))
    return leap_script


def load_ligand_template(leap_system, ligand_type):
    load_molecule_auto(leap_system.input_file, ligand_type.amber_template_alias, ligand_type.file_path)


def create_atom(leap_script, atom, nanoparticle=None, leap_residue_number=1, coords_mult=1000.):
    atom_leap_name = '{}_{}'.format(atom.element, atom.amber_name)
    leap_script.write('{} = createAtom {} {} {}\n'.format(atom.amber_name, atom_leap_name, atom.element, atom.charge))
    leap_script.write('set {} element {}\n'.format(atom.amber_name, atom.element))
    # if type(atom.coords[0]) is int or type(atom.coords[0]) is np.int:
    #    leap_script.write('set {0} position {{ {1[0]} {1[1]} {1[2]} }}\n'.format(atom.amber_name, atom.coords / 1000.))
    # else:
    leap_script.write('set {0} position {{ {1[0]} {1[1]} {1[2]} }}\n'.format(atom.amber_name, atom.coords / coords_mult))
    if nanoparticle:
        leap_script.write('add {}.{} {} \n'.format(nanoparticle.amber_unit_name, leap_residue_number, atom.amber_name))
    return leap_script


def create_unit(leap_script, unit_alias, unit_name):
    leap_script.write('{} = createUnit {}\n'.format(unit_alias, unit_name))
    return leap_script


def create_residue(leap_script, residue_alias, residue_name, unit_alias=None):
    leap_script.write('{} = createResidue {}\n'.format(residue_alias, residue_name))
    if unit_alias:
        leap_script.write('add {} {}\n'.format(unit_alias, residue_alias))
    return leap_script


def auto_atom_creator(leap_script, coords_list, leap_unit_alias, leap_residue_number, element_symbol_list,
                      atom_leap_name_list=[], atom_charge_list=None, atom_letter_prefix_list=[], coords_multiplier=1):
    # If atom names and prefixes for atom names are not provided, assign atom names based on first character of element.
    # If number runs over max space(1000 for a###), use other case and then increment through letters not being used.
    import sys
    if not atom_charge_list:
        atom_charge_list = [0] * len(coords_list)
    if len(atom_leap_name_list) == 0:
        if len(coords_list) < np.sum([10 ** (4 - len(prefix)) for prefix in atom_letter_prefix_list]):
            atom_leap_name_list = ['{}{}'.format(atom_letter_prefix, atom_number) for atom_number, atom_letter_prefix in
                                   enumerate(atom_letter_prefix_list)]
        else:
            print(
                'More than 1000 atoms passed without a list of prefixes for Leap names. Pass a tuple to {} containing '
                'strings. Example: pass ["g", "G"] for Au atoms.'.format(sys._getframe().f_code.co_name))
    for atom_py_index, (atom_coords, atom_leap_name, atom_element_symbol, atom_charge) in enumerate(
            zip(coords_list, atom_leap_name_list, element_symbol_list, atom_charge_list)):
        coords_angstrom = [xyz / coords_multiplier for xyz in atom_coords]
        atom_alias = '{}{}'.format(atom_element_symbol, atom_py_index)
        atom_leap_name = '{}{}'.format(atom_letter_prefix, atom_letter_number)
        atom_letter_number += 1
        if atom_letter_number == 10 ** (4 - len(atom_letter_prefix)):
            atom_letter_prefix = atom_letter_prefix_list.pop()
            atom_letter_number = 0
    return leap_script

    # Bond NP atoms to each other.


def create_nanoparticle_core(leap_system, nanoparticle):
    # Create containers in Leap.
    # try:
    #    nano_unit_alias = nanoparticle.amber_unit_name
    #except:
    nano_unit_alias = nanoparticle.amber_unit_name
    nano_unit = '{}_unit'.format(nano_unit_alias)
    nano_residue = '{}_res'.format(nano_unit_alias)
    nano_residue_alias = '{}_res_alias'.format(nano_unit_alias)
    create_unit(leap_system.input_file, unit_alias=nano_unit_alias, unit_name=nano_unit)
    create_residue(leap_system.input_file, residue_alias=nano_residue_alias, residue_name=nano_residue,
                   unit_alias=nano_unit_alias)
    # Create all atoms.
    for atom in nanoparticle.atoms:
        create_atom(leap_system.input_file, atom, nanoparticle=nanoparticle)
    return

def copy_ligand(leap_system, ligand_object, nanoparticle=None):
    if ligand_object not in leap_system.residue_objects:
        leap_system.create_residue(ligand_object, amber_alias='ligand{}'.format(len(leap_system.residue_objects)))
    leap_system.input_file.write(
        '{}_unit = copy {}\n'.format(leap_system.residue_objects[ligand_object], ligand_object.amber_template_alias))
    leap_system.input_file.write('{0} = {0}_unit.1\n'.format(leap_system.residue_objects[ligand_object]))
    if nanoparticle:
        leap_system.input_file.write('remove {0}_unit {0}\n'.format(leap_system.residue_objects[ligand_object]))
        leap_system.input_file.write(
            'add {} {}\n'.format(nanoparticle.amber_unit_name, leap_system.residue_objects[ligand_object]))
    return


def move_ligands(leap_system, nanoparticle, coords_multiplier=1000., throw=4., rotate_lig=False, ligand_vector=np.array([0,0,-1.])):
    for site, ligand in nanoparticle.site_ligand_dict.items():
        # if type(ligand.site.coords[0]) is int or type(ligand.site.coords[0]) is np.int:
        coords = site.coords / coords_multiplier * throw
        if rotate_lig:
            rotation_vector = geometry_tools.rotation_matrix_from_vectors(ligand_vector, coords)
            #print(rotation_vector)
            #print(np.linalg.norm(np.array(rotation_vector)))
            #leap_system.input_file.write('translate {0} {{ 0 0 30 }}\n'.format(ligand.amber_alias))
            leap_system.input_file.write('transform {0} {{ {{ {1[0][0]} {1[0][1]} {1[0][2]} }} {{ {1[1][0]} {1[1][1]} {1[1][2]} }} {{ {1[2][0]} {1[2][1]} {1[2][2]} }} }}\n'.format(ligand.amber_alias, rotation_vector))
            #leap_system.input_file.write('translate {0} {{ 0 0 -30 }}\n'.format(ligand.amber_alias))
        leap_system.input_file.write('translate {0} {{ {1[0]} {1[1]} {1[2]} }}\n'.format(ligand.amber_alias, coords))
    return

# Loads all ligand templates, makes copies, moves to nanoparticle unit, and bonds ligands to surface atoms.
def add_all_ligands(leap_system, nanoparticle):
    for ligand_type in nanoparticle.ligand_types:
        load_ligand_template(leap_system, ligand_type)
    for ind, (site, ligand) in enumerate(nanoparticle.site_ligand_dict.items()):
        if ligand.amber_alias is None:
            setattr(ligand, 'amber_alias', 'ligand{}'.format(len(leap_system.residue_objects)))
        copy_ligand(leap_system, ligand, nanoparticle=nanoparticle)
    move_ligands(leap_system, nanoparticle, throw=4.)
    return

def bond_atoms_by_aliases(leap_script, leap_alias_pairs_tuple, bond_order=1):
    for alias_list in leap_alias_pairs_tuple:
        if bond_order == 1:
            leap_script.write('bond {0[0]} {0[1]}\n'.format(alias_list))
        else:
            leap_script.write('bond {0[0]} {0[1]} {1}\n'.format(alias_list, bond_order))
    return leap_script


def bond_atoms_by_number(leap_script, leap_atom_numbers_list, unit_number_list=[1, 1], residue_number_list=[1, 1],
                         bond_order=1):
    for alias_list in leap_atom_numbers_list:
        if bond_order == 1:
            leap_script.write(
                'bond {0[0]}.{1[0]}.{2[0]} {0[0]}.{1[0]}.{2[1]}\n'.format(unit_number_list, residue_number_list,
                                                                          leap_atom_numbers_list))
        else:
            leap_script.write(
                'bond {0[0]}.{1[0]}.{2[0]} {0[0]}.{1[0]}.{2[1]} {3}\n'.format(unit_number_list, residue_number_list,
                                                                              leap_atom_numbers_list, bond_order))
    return leap_script


def bond_all_nanoparticle_atoms(leap_system, nanoparticle):
    # Bond all ligands.
    # print('Bonding all atoms')
    for site_atom, ligand in nanoparticle.site_ligand_dict.items():
        if [site_atom, ligand] not in leap_system.bond_pairs_list and [ligand,
                                                                       site_atom] not in leap_system.bond_pairs_list:
            leap_system.bond_atoms(site_atom, ligand)
            # print('Ligand-sulfur: {}.{}'.format(ligand.ligand_type.mol2_name, ligand.attach_index))
            leap_system.input_file.write(
                'bond {0} {1}.{2}\n'.format(site_atom.amber_name, leap_system.residue_objects[ligand],
                                            ligand.attach_index))
    # Bond all core atoms, checking for duplicates.
    for atom1 in nanoparticle.atoms:
        for atom2 in atom1.bonded_atoms:
            if [atom1, atom2] in leap_system.bond_pairs_list or [atom2, atom1] in leap_system.bond_pairs_list:
                continue
            else:
                leap_system.bond_atoms(atom1, atom2)
                leap_system.input_file.write('bond {0} {1}\n'.format(atom1.amber_name, atom2.amber_name))
    return


def add_ions(leap_system, ions_num_dict=None, neutralize=True):
    # Neutralize system.
    if neutralize:
        leap_system.input_file.write('addions {} Cl- 0\n'.format(leap_system.amber_unit_name))
        leap_system.input_file.write('addions {} Na+ 0\n'.format(leap_system.amber_unit_name))
    if ions_num_dict:
        for ion in ions_num_dict:
            leap_system.input_file.write('addions {} {} {}\n'.format(leap_system.amber_unit_name, ion, ions_num_dict[ion]))
    return


def save_amber_parm(leap_script, unit_alias, system_name=None, parm_path=''):
    if not system_name:
        system_name = unit_alias
    parm_file_name = parm_path + system_name
    leap_script.write('saveamberparm {0} {1}.prmtop {1}.rst7\n'.format(unit_alias, parm_file_name))
    return leap_script


def assign_atom_names(num_atoms, prefix_list):
    atom_names_list = []
    max_names = [10 ** (4 - len(prefix)) for prefix in prefix_list]
    if num_atoms <= max_names[0]:
        atom_names_list = [prefix_list[0] + format(index, '03') for index in list(range(num_atoms))]
        return atom_names_list
    cumulative_names = np.cumsum(max_names)
    full_prefix_bool_list = cumulative_names < num_atoms
    for prefix_num, prefix in enumerate(itertools.compress(prefix_list, full_prefix_bool_list)):
        #print(prefix)
        atom_names_list.extend([prefix + format(index, '03') for index in list(range(max_names[prefix_num]))])
    partial_prefix_index = np.size(np.nonzero(cumulative_names < num_atoms))
    num_in_partial = num_atoms - max_names[partial_prefix_index]
    prefix = prefix_list[prefix_num + 1]
    atom_names_list.extend([prefix + format(index, '03') for index in list(range(num_in_partial))])
    return atom_names_list


def remove_file_extension(file_name, file_extensions=[]):
    for extension in file_extensions:
        if not extension[0] == '.':
            extension = '.' + extension
        length = len(extension)
        if file_name[-length:] is extension:
            new_name = file_name[:-length]
            return new_name


def master_file(systems_list, enviro='#!/usr/bin/env bash', script_path='', warning_grep=False, logs=True):
    # ## Create master file for tLeap execution.
    grep_list = ['max', 'improper', 'missing', 'WARNING -A 1']
    script_path = '{}LeapScripts/'.format(script_path)
    master_shell_file = open('{}master_script.sh'.format(script_path), 'w')
    master_shell_file.write('rm leap.log\n')
    master_shell_file.write('rm error.log\n')
    master_shell_file.write('{}\n\n'.format(enviro))
    for system in systems_list:
        master_shell_file.write('tleap -f {} > leap.log\n'.format(
            system.leap_filename))  # | grep -v "maximum number of bonds" | grep -v "triangular"\n'.format(system.leap_filename))
        master_shell_file.write('echo {} >> error.log\n'.format(system.system_name))
        master_shell_file.write('grep "type String" leap.log >> error.log\n')  # .format(system.system_name)))
        #'tleap -f {} | grep -v "There is a bond of" | grep -v "triangular" \n'.format(system.leap_filename))
        master_shell_file.write('mv leap.log {}.log \n'.format(system.system_name))
        master_shell_file.write('mv {} LeapIn \n'.format(system.leap_filename))
        master_shell_file.write('mv {}.log logs \n\n'.format(system.system_name))
    if warning_grep:
        for system, search_string in itertools.product(systems_list, grep_list):
            master_shell_file.write(
                'grep -v "There is a bond of" {1}.log | grep -v "maximum number of bonds" | grep -v "triangular" | grep {0}.log \n'.format(
                    search_string, system.leap_filename))
            master_shell_file.write('sleep 2\n')
    #for system in systems_list:
    master_shell_file.close()
    print('Master shell script created: {}master_script.sh. \n'.format(script_path))
    return


# Create .xyz file from NanoParticle instance. Coordinates should be in milli-Angstrom.
def xyz_writer(nanoparticle, file_name):
    with open(file_name, 'w') as xyz_file:
        xyz_file.write('{}\n\n'.format(len(nanoparticle.atoms)))
        for atom in nanoparticle.atoms:
            # Convert to Angstroms.
            print(atom.coords)
            coords = atom.coords / 1000.
            xyz_file.write('{} {:.4f} {:.4f} {:.4f} \n'.format(atom.element, coords[0], coords[1], coords[2]))
    return


# Main function for writing tLeap input files for nanoparticle construction.
def write_leap_input(system, nanoparticle, force_field_list=None, prmtop_name=None, prmtop_path=default_parm_path):
    print('Writing tLeap input')
    if force_field_list is None:
        force_field_list = ['gaff2', 'water.tip3p', 'protein.ff14SB']
    frcmod_names = system.get_frcmods()
    for ligand_type in nanoparticle.ligand_types:
        frcmod_names.append('{}frcmod/{}.frcmod'.format(ligand_file_path, ligand_type.mol2_name))
    #print(frcmod_names)
    system.open_for_writing()
    load_force_fields(system.input_file, force_field_list=force_field_list, frcmod_list=frcmod_names)
    # Create nanoparticle core.
    create_nanoparticle_core(leap_system=system, nanoparticle=nanoparticle.generic_parent.parent_lattice)
    # Place all ligands.
    add_all_ligands(leap_system=system, nanoparticle=nanoparticle)
    bond_all_nanoparticle_atoms(leap_system=system, nanoparticle=nanoparticle)
    add_ions(leap_system=system)
    save_amber_parm(leap_script=system.input_file,
                    unit_alias=nanoparticle.generic_parent.parent_lattice.amber_unit_name, system_name=prmtop_name,
                    parm_path=prmtop_path)
    system.write('quit\n')
    system.close_input_file()
    return


'''
def leap_script_writer(coords_list, bond_pairs_list, site_coords_index_tuple,
                       nanocluster_name, script_path, ligand_names_list, sulfur_index_list):
    print(site_coords_index_tuple)

    leap_script = open(script_path, 'w')
    leap_script.write('verbosity 0\n')
    # Load parameters.
    forcefield_list = ['gaff2', 'water.tip3p']
    frcmod_list = ['AuNP.frcmod']  # /home/mdmannin/amber16/dat/leap/parm/AuNP.frcmod
    for ligand_name, site_coords_index_list in zip(ligand_names_list, site_coords_index_tuple):
        if len(site_coords_index_list) == 0:
            continue
        print('Atoms bound to {}: {}.'.format(ligand_name, len(site_coords_index_list)))
        frcmod_list.append('{}frcmod/{}.frcmod'.format(ligand_file_path, ligand_name))
    print(script_path)
    load_force_fields(leap_script, forcefield_list, frcmod_list)
    # ## Create gold residue and atoms and move to coordinates in list.
    core_unit_alias = 'nano_unit'
    gold_residue_alias = 'gold_res'
    create_unit(leap_script, core_unit_alias, 'goldunit')
    create_residue(leap_script, gold_residue_alias, 'goldcore', unit_alias=core_unit_alias)
    # Create atoms in NP and bond to each other.
    leap_atom_names = assign_atom_names(len(coords_list), ['g', 'a', 'i', 'u'])
    site_alias_tuple = []
    for ligand_index_list in site_coords_index_tuple:
        site_alias_tuple.append([leap_atom_names[ind] for ind in ligand_index_list])
    for atom_alias, atom_coords_int in zip(leap_atom_names, coords_list):
        atom_leap_name = atom_alias + '_name'
        atom_coords_angstrom = [x / 1000. for x in atom_coords_int]
        create_atom(leap_script, atom_alias, atom_leap_name, 'Au', 0, atom_coords_angstrom,
                    core_unit_alias, 1)
    bond_alias_tuple = [[leap_atom_names[index[0]], leap_atom_names[index[1]]] for index in bond_pairs_list]

    # leap_script = bond_atoms_by_aliases(leap_script, bond_alias_tuple)
    # print(site_alias_tuple)
    leap_script = bond_atoms_by_number(leap_script, bond_pairs_list)
    ligand_files_list = ['{}.mol2'.format(ligand_name) for ligand_name in ligand_names_list]
    leap_ligand_placer(leap_script, core_unit_alias, coords_list, site_coords_index_tuple,
                       site_alias_tuple, ligand_file_path, ligand_files_list, sulfur_index_list,
                       ligand_names_list=ligand_names_list, previous_residue_number=1,
                       ligand_coords_multiplier=1000, core_residue_number=1)

    leap_script.write('setBox nano_unit vdw 5\n')
    leap_script.write('addionsrand nano_unit Cl- 0\n')
    # leap_script.write('solvateBox nanorod TIP3PBOX {}\n'.format(int(solvent_buffer)))
    save_amber_parm(leap_script, unit_alias=core_unit_alias, system_name=nanocluster_name,
                    parm_path=parm_path)
    # leap_script.write('check nanorod\n')
    leap_script.write('quit\n')
    leap_script.close()
    print('Tleap input file ready!.\n')

    return
'''

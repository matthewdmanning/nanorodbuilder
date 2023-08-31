import collections
import copy
import itertools
import random
import string

import numpy as np


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


def load_forcefields(leap_script, forcefield_list, frcmod_list):
    for forcefield in forcefield_list:
        force_field_name = 'leaprc.{}'.format(forcefield)
        if check_amber_file_exists(force_field_name):
            leap_script.write('source leaprc.{}\n'.format(forcefield))
        else:
            print('Force field file, {}, cannot be found in $AMBERHOME directory. Check spelling or try launching this '
                  'script from the command line.'.format(force_field_name))
    for frcmod in frcmod_list:
        if check_amber_file_exists(frcmod):
            leap_script.write('loadamberparams {}\n'.format(frcmod))
        else:
            print(
                '{} cannot be found in $AMBERHOME directory. Check spelling or try launching this script from the command'
                ' line.'.format(frcmod))
    return leap_script


def check_pdbseq_validity(seq):
    residue_string = seq.strip('{}')
    for residue in residue_string.split():
        if residue in amber_recognized_residues:
            continue
        else:
            print('Warning: Unrecognized residue {} in loadpdbUsingSeq command. Check sequence.'.format(residue))
            return False
    return seq


def create_unit(leap_script, unit_alias, unit_name):
    leap_script.write('{} = createUnit {}\n'.format(unit_alias, unit_name))
    return leap_script


def create_residue(leap_script, residue_alias, residue_name, unit_alias=None):
    leap_script.write('{} = createResidue {}\n'.format(residue_alias, residue_name))
    if unit_alias:
        leap_script.write('add {} {}\n'.format(unit_alias, residue_alias))
    return leap_script


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


def bond_atoms(leap_script, leap_unit_name_tuple, leap_residue_number_tuple, leap_atom_name_tuple):
    for leap_unit_name_list, leap_residue_number_list, leap_atom_name_list in bond_pairs_list:
        leap_script.write(
            'bond {0[0]}.{1[0]}.{2[0]} {0[1]}.{1[1]}.{2[1]}\n'.format(leap_unit_name_list, leap_residue_number_list,
                                                                      leap_atom_name_list))
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


def create_atom(leap_script, atom_alias, atom_leap_name, atom_element_symbol, atom_charge=0, atom_coords_angstrom=None,
                leap_unit_alias=None, leap_residue_number=None):
    leap_script.write('{} = createAtom {} {} {}\n'.format(atom_alias, atom_leap_name, atom_element_symbol, atom_charge))
    leap_script.write('set {} element {}\n'.format(atom_alias, atom_element_symbol))
    if atom_coords_angstrom:
        leap_script.write('set {0} position {{ {1[0]} {1[1]} {1[2]} }}\n'.format(atom_alias, atom_coords_angstrom))
    if leap_unit_alias and leap_residue_number:
        leap_script.write('add {}.{} {} \n'.format(leap_unit_alias, leap_residue_number, atom_alias))
    return leap_script


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
        print(prefix)
        atom_names_list.extend([prefix + format(index, '03') for index in list(range(max_names[prefix_num]))])
    partial_prefix_index = np.size(np.nonzero(cumulative_names < num_atoms))
    num_in_partial = num_atoms - max_names[partial_prefix_index]
    prefix = prefix_list[prefix_num + 1]
    atom_names_list.extend([prefix + format(index, '03') for index in list(range(num_in_partial))])
    return atom_names_list


def auto_atom_creator(leap_script, coords_list, leap_unit_alias, leap_residue_number, element_symbol_list,
                      atom_leap_name_list=[], atom_charge_list=None, atom_letter_prefix_list=[],
                      coords_multiplier=1):
    # If atom names and prefixes for atom names are not provided, assign atom names based on first character of element.
    # If number runs over max space(1000 for a###), use other case and then increment through letters not being used.
    import sys
    if not atom_charge_list:
        atom_charge_list = [0] * len(coords_list)
    if len(atom_leap_name_list) == 0:
        if len(coords_list) < np.sum([10 ** (4 - len(prefix)) for prefix in atom_letter_prefix_list]):
            atom_leap_name_list = ['{}{}'.format(atom_letter_prefix, atom_number) for
                                   atom_number, atom_letter_prefix in enumerate(atom_letter_prefix_list)]
        else:
            print(
                'More than 1000 atoms passed without a list of prefixes for Leap names. Pass a tuple to {} containing '
                'strings. Example: pass ["g", "G"] for Au atoms.'.format(sys._getframe().f_code.co_name))
    for atom_py_index, (atom_coords, atom_leap_name, atom_element_symbol, atom_charge) in enumerate(
            zip(coords_list, atom_leap_name_list, element_symbol_list, atom_charge_list)):
        coords_angstrom = [xyz / coords_multiplier for xyz in atom_coords]
        atom_alias = '{}{}'.format(atom_element_symbol, atom_py_index)
        # atom_leap_name = '{}{}'.format(atom_letter_prefix, atom_letter_number)
        atom_letter_number += 1
        if atom_letter_number == 10 ** (4 - len(atom_letter_prefix)):
            atom_letter_prefix = atom_letter_prefix_list.pop()
            atom_letter_number = 0
    return leap_script

    # Bond NP atoms to each other.


def remove_file_extension(file_name, file_extensions=[]):
    for extension in file_extensions:
        if not extension[0] == '.':
            extension = '.' + extension
        length = len(extension)
        if file_name[-length:] is extension:
            new_name = file_name[:-length]
            return new_name


def attach_ligand(leap_script, ligand_residue_alias, ligand_residue_number, ligand_attach_index, ligand_coords_angstrom,
                  parent_unit_alias, site_alias, ligand_unit_alias=None):
    if not ligand_unit_alias:
        ligand_unit_alias = '{}_unit'.format(ligand_residue_alias)
    leap_script.write('translate {0} {{ {1[0]} {1[1]} {1[2]} }}\n'.format(ligand_unit_alias, ligand_coords_angstrom))
    leap_script.write('{} = {}.1\n'.format(ligand_residue_alias, ligand_unit_alias))
    leap_script.write('remove {} {}\n'.format(ligand_unit_alias, ligand_residue_alias))
    leap_script.write('add {} {}\n'.format(parent_unit_alias, ligand_residue_alias))
    leap_script.write(
        'bond {0}.{1}.{2} {3}\n'.format(parent_unit_alias, ligand_residue_number, ligand_attach_index, site_alias))
    return leap_script


def leap_ligand_placer(leap_script, parent_unit_alias, coords_list, site_coords_index_tuple, site_leap_alias_tuple,
                       ligand_path, ligand_files_list, ligand_attach_index_list, ligand_names_list=None, throw=4.0,
                       previous_residue_number=1, ligand_coords_multiplier=1, core_residue_number=1):
    # Check for single ligand ligand_type?
    if len(ligand_names_list) == 0:
        ligand_names_list = [remove_file_extension(ligand_name, ['.mol2', '.pdb']) for ligand_name in ligand_files_list]
    # Loop through ligand types, creating and bonding ligands.
    total_ligand_num = 1
    for ligand_type_num, (
            ligand_name, ligand_file, ligand_attach_index, site_coords_index_list, site_leap_alias_list) in enumerate(
        zip(ligand_names_list, ligand_files_list, ligand_attach_index_list, site_coords_index_tuple,
            site_leap_alias_tuple)):
        if not ligand_name:
            ligand_name = remove_file_extension(ligand_file, ['.mol2', '.pdb'])
        if len(ligand_attach_index) == 0:
            continue
        ligand_file_name = '{}{}'.format(ligand_path, ligand_file)
        # if not check_amber_file_exists(ligand_file_name):
        #     print('Ligand file not found. Please check {}.'.format(ligand_file_name))
        ligand_master_copy_name = 'ligand{}_master'.format(ligand_type_num)
        load_molecule_auto(leap_script, ligand_master_copy_name, ligand_file_name, unit=None, alignaxes=True)
        # leap_script.write('set {}_unit.1 name "{}"\n'.format(ligand_master_copy_name, ligand_name))
        for ligand_type_number, (site_alias, site_coords_index) in enumerate(
                zip(site_leap_alias_list, site_coords_index_list)):
            site_coords = coords_list[site_coords_index]
            ligand_total_residue_num = total_ligand_num + 1
            ligand_residue_alias = '{}{}'.format(ligand_name, ligand_type_number)
            ligand_unit_alias = '{}_unit'.format(ligand_residue_alias)
            ligand_coords_angstrom = [xyz * throw / ligand_coords_multiplier for xyz in site_coords]
            leap_script.write('{} = copy {}\n'.format(ligand_unit_alias, ligand_master_copy_name))
            attach_ligand(leap_script, ligand_residue_alias, ligand_total_residue_num, ligand_attach_index,
                          ligand_coords_angstrom, parent_unit_alias, site_alias, ligand_unit_alias=ligand_unit_alias)
            total_ligand_num += 1
    return leap_script


def master_file(input_files_list, script_names_list, enviro='#/bin/bash', script_path='', warning_grep=True,
                logs=True):
    # ## Create master file for tLeap execution.
    grep_list = ['max', 'improper', 'missing', 'WARNING -A 1']
    script_path = '{}LeapScripts/'.format(script_path)
    master_shell_file = open('{}master_script.sh'.format(script_path), 'w')
    master_shell_file.write(enviro + "\n\n")
    for input_file, script_name in zip(input_files_list, script_names_list):
        master_shell_file.write('rm leap.log\n')
        master_shell_file.write(
            'tleap -f {} | grep -v "There is a bond of" | grep -v "triangular" \n'.format(script_name))
        master_shell_file.write('mv leap.log {}.log\n'.format(input_file))
    for input_file, search_string in itertools.product(input_files_list, grep_list):
        master_shell_file.write(
            'grep -v "There is a bond of" {1}.log | grep -v "maximum number of bonds" | grep -v "triangular" | grep {0} \n'.format(
                search_string, input_file))
        master_shell_file.write('sleep 2\n')
    for input_file in input_files_list:
        master_shell_file.write('mv {}.in LeapIn\n'.format(input_file))
        master_shell_file.write('mv {}.log logs\n'.format(input_file))
    master_shell_file.close()
    print('Master shell script created: {}master_script.sh.\n'.format(script_path))
    return master_shell_file


def xyz_writer(file_name, coords_tuple, element_list=['Au']):
    if len(element_list) == 1 or type(element_list) is str:
        element_list = element_list * len(coords_tuple)
        print(len(coords_tuple))
    elif len(element_list) < len(coords_tuple):
        print(
            'Number of element names does not match the number of coordinates. Repeating given element list: {}'.format(
                element_list))
    elif len(element_list) > len(coords_tuple):
        print('Too many element names given for number of coordinates. Using the first {} names'.format(
            len(coords_tuple)))
    xyz_file = open(file_name, 'w')
    xyz_file.write('{}\n\n'.format(len(coords_tuple)))
    for row, element_symbol in zip(coords_tuple, element_list):
        xyz_file.write('{} {:.4f} {:.4f} {:.4f} \n'.format(element_symbol, row[0], row[1], row[2]))
    return


def leap_script_writer(coords_list, bond_pairs_list, site_coords_index_tuple,
                       nanocluster_name, script_path, ligand_names_list, sulfur_index_list):
    print(site_coords_index_tuple)
    file_path = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/'
    ligand_file_path = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/'
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
    load_forcefields(leap_script, forcefield_list, frcmod_list)
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

    parm_path = "{}ParmFiles/".format(file_path)
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

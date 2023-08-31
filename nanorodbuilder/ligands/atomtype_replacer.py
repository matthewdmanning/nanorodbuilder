#!/usr/bin/env python
import os
import subprocess

# #execfile('/home/mdmannin/amber16/amber.sh')
# #subprocess.call(['source', 'home/mdmannin/amber16/amber.sh'])

# batch_file_name = 'batch_input.txt'
# batch_file = open(batch_file_name, 'r')
main_dir_with_backslash = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/'
input_mol2_file_dir = 'home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/'
ac_output_file_dir = 'home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/Renamed/'
renamed_mol2_file_dir = 'home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/Renamed/'


def rewriter(input_mol2_file_dir, ligand_name, renamed_mol2_file_dir='', ac_output_dir=''):
    if ac_output_dir == '':
        ac_output_dir = '{}Renamed/'.format(input_mol2_file_dir)
    if renamed_mol2_file_dir == '':
        renamed_mol2_file_dir = '{}Renamed/'.format(input_mol2_file_dir)
    input_mol2_file_name = '/{}{}.mol2'.format(input_mol2_file_dir, ligand_name)
    input_mol2_file = open(input_mol2_file_name, 'r')
    ac_output_file_name = '/{}{}.ac'.format(ac_output_dir, ligand_name)
    ac_output_file = open(ac_output_file_name, 'r')
    renamed_mol2_file_name = '/{}{}.mol2'.format(renamed_mol2_file_dir, ligand_name)
    renamed_mol2_file = open(renamed_mol2_file_name, 'w')
    amber_atom_types_list = []
    for line in ac_output_file:
        # print line[:4] == 'ATOM'
        if line[:4] == 'ATOM':
            new_atom_type = line.strip()[-2:]
            if new_atom_type == ' s': new_atom_type = 'ss'
            amber_atom_types_list.append(new_atom_type)
            # print line[-3:]
    renamed_atom_counter = 0
    input_mol2_file.seek(0)
    atoms_trigger = 0
    DUnum = 0


print('The size of the renamed list is {}.\n'.format(len(amber_atom_types_list)))
print(amber_atom_types_list)
for line in input_mol2_file:
    if line[:13] == '@<TRIPOS>ATOM':
        atoms_trigger = 1
        renamed_mol2_file.write(line)
        continue
    if line[:13] == '@<TRIPOS>BOND':
        atoms_trigger = 0
    if atoms_trigger == 1:
        print(line.split())
        old_atom_type = line.split()[5]
        if amber_atom_types_list[renamed_atom_counter].find('DU') == 1:
            newline = line.replace('DU', 'hc')
            DUnum = DUnum + 1
        else:
            # newline = line.replace(old_atom_type, amber_atom_types_list[renamed_atom_counter])
            newline = line[:10] + line[10:].replace(old_atom_type, amber_atom_types_list[renamed_atom_counter])
            print(newline)
        renamed_mol2_file.write(newline)
        renamed_atom_counter = renamed_atom_counter + 1
    else:
        renamed_mol2_file.write(line)
    # Leave off all info related to substructure,formal charges, etc.
    # if line[:21] == '@<TRIPOS>SUBSTRUCTURE': break
# print '{} atoms reassigned as from DU to hc.'.format(DUnum)
renamed_mol2_file.close()
ac_output_file.close()
input_mol2_file.close()
return


def get_file_names_from_batch_input(batch_file):
    lengthlist, radiuslist, densitylist, ligandlist, ratiolist = [], [], [], [], []
    for line in batch_file:
        line = line.strip()
        if line.find("LENGTH=") != -1:
            lengthbite = line[line.index('=') + 1:].split(',')
            for i in lengthbite: lengthlist.append(int(i))
        elif line.find("RADIUS=") != -1:
            radiusbite = line[line.index('=') + 1:].split(',')
            for i in radiusbite: radiuslist.append(int(i))
        elif line.find("DENSITY=") != -1:
            densitybite = line[line.index('=') + 1:].split(',')
            for i in densitybite: densitylist.append(int(i))
        elif line.find("LIGAND1=") != -1:
            ligandbite = line[line.index('=') + 1:].split(',')
            for i in ligandbite: ligandlist.append(i)
        elif line.find("RATIO=") != -1:
            ratiobite = line[line.index('=') + 1:].split(',')
            for i in ratiobite: ratiolist.append(int(100 * float(i)))
    return lengthlist, radiuslist, densitylist, ligandlist, ratiolist
    # ===========================================================================
    # elif line.find("LIGANDS2=") == 1:
    #     ligand2list = line[line.index('=')+1,:].split(delimiter=',')
    # elif line.find("SULF1INDEX=") == 1:
    #     sulf1indexlist = line[line.index('=')+1,:].split(delimiter=',')
    # elif line.find("SULF2INDEX=") == 1:
    #     sulf2indexlist = line[line.index('=')+1,:].split(delimiter=',')
    # ===========================================================================


# for length, radius, density, ligand, ratio in itertools.product(lengthlist, radiuslist, densitylist, ligandlist, ratiolist):
#    mol2_name = '{}x{}-{}ratio{}dense{}'.format(length, radius, ratio, density, ligand)

def write_antechamber_retype_script(mol2_input_dir, ligand_name, ac_output_file_path='', script_file_path='',
                                    ac_script_file_name=''):
    if ac_output_file_path == '':
        ac_output_file_path = '{}Renamed/'.format(mol2_input_dir)
    if ac_script_file_name == '':
        ac_script_file_name = '{}renamingscript.sh'.format(script_file_path)
    ac_script_file = open(ac_script_file_name, 'w+')
    ac_script_file.write('#!/bin/sh\n')

    mol2_file_name = '{}{}.mol2'.format(mol2_input_dir, ligand_name)
    ac_output_file_name = '{}{}_rn.ac'.format(ac_output_file_path, ligand_name)
    subprocess.call(
        ['home/mdmannin/amber16/bin/antechamber', '-i ' + mol2_file_name, '-f mol2', '-o ' + ac_output_file_name,
         '-p gaff2'])
    # ac_script_file.write('atomtype -i ' + mol2_file_name + ' -f mol2 -o ' + ac_output_file_name + ' -p gaff2 \n')
    ac_script_file.close()
    return ac_script_file_name


def run_script_from_python(script_file_dir, script_name='', script_file_name=''):
    import stat
    if script_file_name == '':
        if script_name == '':
            script_name = 'renaming_script.sh'
        script_file_name = '{}{}'.format(script_file_dir, script_name)
    st = os.stat(script_file_name)
    os.chmod(script_file_name, st.st_mode | stat.S_IEXEC)
    rename = subprocess.Popen(['source {}'.format(script_file_name)])
    rename.wait()
    return


# for ligand_name in [f.replace('.mol2', '') for f in os.listdir((input_mol2_file_dir)) if f.endswith('.mol2')]:# and f.replace('.mol2', '.ac') not in os.listdir(renamed_mol2_file_dir)]:
#    antechamber_script_file_name = write_antechamber_retype_script(input_mol2_file_dir, ligand_name, script_file_path=input_mol2_file_dir)
# run_script_from_python(input_mol2_file_dir, script_file_name=antechamber_script_file_name)
for ligand_name in [f.replace('.mol2', '') for f in os.listdir(main_dir_with_backslash) if f.endswith('.mol2')]:
    rewriter(input_mol2_file_dir, ligand_name)

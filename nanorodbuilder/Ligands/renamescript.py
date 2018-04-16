#!/usr/bin/env python
import glob

infilepath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/'
outfilepath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/Renamed/'
renamedfilepath = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/Renamed/'

###Replace atom names in mol2 file.
glob_pattern = '3C12'
ligandlist = glob.glob('{}{}.mol*'.format(infilepath, glob_pattern))
for ligand_file in ligandlist:
#	try:
ligand_name = ligand_file.replace(infilepath, '').replace('.mol2', '')
red_filename = '{}{}.mol2'.format(infilepath, ligand_name)
ante_filename = '{}{}.ac'.format(outfilepath, ligand_name)
red_file = open(red_filename, 'r')
ante_file = open(ante_filename, 'r')
renamed_filename = '{}{}_rn.mol2'.format(renamedfilepath, ligand_name)
renamed_file = open(renamed_filename, 'w')
print('Renamed file is {}'.format(renamed_filename))
amber_name_list = []
for line in ante_file:
    # print line[:4] == 'ATOM'
    if line[:4] == 'ATOM':
        ac_line_split = line.split()
        # Amber atom type is last column in antechamber output files.
        new_atom_name = ac_line_split[-1]
        if new_atom_name == 's':
            new_atom_name = 'ss'
        amber_name_list.append(new_atom_name)
        # print line[-3:]
atom_index = 0
red_file.seek(0)
trigger = 0
print('The size of the renamed list is {}.\n'.format(len(amber_name_list)))
for line in red_file:
    if '@<TRIPOS>ATOM' in line:
        trigger = 1
        renamed_file.write(line)
        continue
    if '@<TRIPOS>BOND' in line:
        trigger = 0
    if trigger == 1:
        # if amber_name_list[atom_index] == 's':
        #    newline = line[:43]+'ss'+line[46:]
        #    print(newline)
        # else:
        newline = line[:43] + amber_name_list[atom_index] + line[46:]
        renamed_file.write(newline)
        atom_index += 1
    else:
        renamed_file.write(line)
        # Leave off all info related to substructure,formal charges, etc.
        # if line[:21] == '@<TRIPOS>SUBSTRUCTURE': break
renamed_file.close()
ante_file.close()
red_file.close()
#    except IOError:
#        print 'File not found.'
#        continue

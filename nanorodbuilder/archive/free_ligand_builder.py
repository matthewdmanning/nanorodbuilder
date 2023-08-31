import random

from nanorodbuilder.utils import ligandsysbuilder, script_writing_tools

leap_script_name = "leap.in"
nano_file_name = '19A-cluster-120_DEA.mol2'
ligand_mol2_path = "/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/"
ligand_frcmod_path = "/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/frcmod/"
force_field_list = ['gaff2', 'water.tip3p']
lig_name_list = ['DEAF', 'DEAUF']
lig_num_list = [50, 50]
min_buffer_Ang = 25
max_throw_range = 10
lig_file_name_list = ["{}{}.mol2".format(ligand_mol2_path, ligand_name) for ligand_name in lig_name_list]
lig_frcmod_name_list = ["{}{}.frcmod".format(ligand_frcmod_path, lig_name_list) for ligand_name in lig_name_list]
lig_name_mol2_dict = {key: value for (key, value) in zip(lig_name_list, lig_file_name_list)}

leap_script = open(leap_script_name, 'w')
script_writing_tools.load_forcefields(leap_script, force_field_list, frcmod_list=lig_frcmod_name_list)
# Load ligands.
total_ligand_num = sum(lig_num_list)
ligand_iter_list = [[name] * number for name, number in zip(lig_name_list, lig_num_list)][0]
random.shuffle(ligand_iter_list)
print(ligand_iter_list)
ligand_alias_list = []
for ligand_count, ligand_name in enumerate(ligand_iter_list):
    ligand_alias = "{}{}".format(ligand_name, ligand_count)
    ligand_alias_list.append(ligand_alias)
    mol2_file_name = lig_name_mol2_dict[ligand_name]
    translate_Ang = random.randrange(1, max_throw_range) + min_buffer_Ang
    leap_script = ligandsysbuilder.ligand_mover(translate_Ang, leap_script, ligand_alias, mol2_file_name)
nano_alias = 'nano'
leap_script.write('{} = loadmol2 {}\n'.format(nano_alias, nano_file_name))
ligand_alias_string = ' '.join(ligand_alias_list)
print(ligand_alias_string)
leap_script.write('system = combine {{ {} {} }}\n'.format(ligand_alias_string, nano_alias))
leap_script.close()

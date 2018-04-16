import NanoParticle
import System
# Custom packages
import geometry_tools
import script_writing_object

script_path = "/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/"


def parameters():
    radius_list, charged_num_list, ligand_num_list, ligand1_list, ligand2_list = [], [], [], [], []
    params_dict = {"RADIUS": radius_list, "CHARGED_LIGANDS": charged_num_list, "TOTAL_LIGANDS": ligand_num_list,
                   "LIGAND1": ligand1_list, "LIGAND2": ligand2_list}
    ligand_name_info_dict = {}
    batch_filename = 'cluster_input.txt'
    with open(batch_filename, 'r') as batch_file:
        for line in batch_file:
            for keyword, param_list in params_dict.items():
                if keyword in line:
                    line_bite = line[line.index('=') + 1:].split(',')
                    for value in line_bite:
                        param_list.append(value.rstrip())
    radius_list = [int(num) for num in radius_list]
    charged_num_list = [int(num) for num in charged_num_list]
    ligand_num_list = [int(num) for num in ligand_num_list]
    for ligand_name in ligand1_list + ligand2_list:
        if ligand_name not in ligand_name_info_dict.keys():
            new_ligand = NanoParticle.LigandInfo(ligand_name)
            ligand_name_info_dict[ligand_name] = new_ligand
    ligand_info_filename = "ligand_info.txt"
    with open(ligand_info_filename, 'r') as batch_file:
        for line in batch_file:
            for ligand_name, ligand_info in ligand_name_info_dict.items():
                if line.find(ligand_name) != -1:
                    ligand_info_string = line[line.index(ligand_name) + len(ligand_name) + 1:].split()
                    sulfur_index = ligand_info_string[0]
                    ligand_info.attach_index = sulfur_index
    return radius_list, ligand_num_list, charged_num_list, ligand1_list, ligand2_list, ligand_name_info_dict

# Output two lists: one with the distances of the nearest neighbors and another with the number of atoms at that distance.

def main():
    # radius_list, charged_num_list, ligand_num_list, ligand1_list, ligand2_list, sulf1_index_list, sulf2_index_list = initializer()
    radius_list, ligand_num_list, charged_num_list, ligand1_list, ligand2_list, ligand_name_info_dict = parameters()
    nanocluster_names_list = []
    systems_list = []
    for radius_angstroms in radius_list:
        # Get lists of coordinates of gold atoms.
        xyz_dir = '/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/NanorodXYZFiles/'
        xyz_file_name = '{}{}A-cluster.xyz'.format(xyz_dir, radius_angstroms)
        ### coords_list must be a 2-deep tuple with len(coords_list[i]) == 3 shell_index_list and core_index_list should be lists
        lattice = geometry_tools.nanoparticle_builder(radius_angstroms, element='Au', unit_cell_length_int=4070)
        lattice.get_surface_atoms()
        script_writing_object.xyz_writer(lattice, xyz_file_name)
        for charged_num, ligand_num in zip(charged_num_list, ligand_num_list):
            if ligand_num > len(lattice.surface_atoms):
                print("Number of ligands exceeds number of surface atoms. Script not written.")
                continue
            uncharged_num = ligand_num - charged_num
            ligand_num_list = [charged_num, uncharged_num]
            generic_np = NanoParticle.GenericParticle(lattice)
            generic_np.choose_unspecified_sites(ligand_num)
            generic_np.clear_placeholder_sites()
            # print(ligand_num_list)
            generic_np.choose_placeholder_sites(ligand_num_list=ligand_num_list)
            for ligand1_name, ligand2_name in zip(ligand1_list, ligand2_list):
                if uncharged_num > 0 and charged_num > 0:
                    nanocluster_name = '{0}A-{1}_{2}-{3}_{4}'.format(radius_angstroms, charged_num,
                                                                             ligand1_name, uncharged_num, ligand2_name)
                    ligand_names_list = [ligand1_name, ligand2_name]
                elif uncharged_num == 0:
                    nanocluster_name = '{0}A-{1}_{2}'.format(radius_angstroms, charged_num, ligand1_name)
                    ligand_names_list = [ligand1_name]
                elif charged_num == 0:
                    nanocluster_name = '{0}A-{1}_{2}'.format(radius_angstroms, uncharged_num, ligand2_name)
                    ligand_names_list = [ligand2_name]
                complete_nanocluster = NanoParticle.Nanoparticle(generic_np)
                # ligand_info_list = [ligand_name_info_dict[name] for name in ligand_names_list]
                ligand_info_num_dict = {ligand_name_info_dict[ligand1_name]: charged_num,
                                        ligand_name_info_dict[ligand2_name]: uncharged_num}
                complete_nanocluster.create_ligand_sites(ligand_info_num_dict=ligand_info_num_dict)
                leap_filename = '{}LeapScripts/{}_leap.in'.format(script_path, nanocluster_name)
                new_system = System.LeapSystem(leap_filname=leap_filename, amber_unit_name=nanocluster_name,
                                               system_name=nanocluster_name)
                systems_list.append(new_system)
                print('Saving script to {}.'.format(leap_filename))
                script_writing_object.write_leap_input(new_system, complete_nanocluster, prmtop_name=nanocluster_name)

            script_writing_object.master_file(systems_list,
                                              script_path='/home/mdmannin/Desktop/Nanoparticles/RoughNanorods/')


if __name__ == "__main__":
    main()

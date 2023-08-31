mport
os


# Calculate angle b/w Au CoM-sulfur-functional/head group (eg. n4)
def ligand_np_angle(trajin_file_object, traj_name, ligand_residue_num_list, au_core_residue_num, head_group_atom_name,
                    sulfur_atom_name='ss'):
    # Define prefixes for dataset names.
    au_com_name = 'au_com'
    au_com_sulfur_name = 'au_com_sulfur'
    sulfur_head_name = 'sulfur_head'
    au_sulfur_head_angle_name = 'au_sulfur_head_angle'
    trajin_file_object.write('vector au_com center out au_com.{}.dat\n'.format(traj_name))
    for ligand_residue_num in ligand_residue_num_list:
        # Calculate CoM from Au core to sulfur. This will store vectors from all ligands in the same data file.
        trajin_file_object.write(
            'vector {0}{1} mask out {0}.dat :{2} :{1}@{3}\n'.format(au_com_sulfur_name, ligand_residue_num,
                                                                    au_core_residue_num, sulfur_atom_name))
        # Calculate sulfur-head vector. This will store vectors from all ligands in the same data file.
        trajin_file_object.write(
            'vector {0}{1} mask out {0}.dat :{1}@{2} :{1}@{3}\n'.format(sulfur_head_name, ligand_residue_num,
                                                                        sulfur_atom_name, head_group_atom_name))
        trajin_file_object.write(
            'vectormath vec1 {0}{1} vec2 {2}{1} out {3}.dat name {3}{1} dotangle\n'.format(au_com_sulfur_name,
                                                                                           ligand_residue_num,
                                                                                           sulfur_head_name,
                                                                                           au_sulfur_head_angle_name))
        # return trajin_file_object


# Calculates RMSD of headgroup with symmetry correction based on sulfur atoms. This eliminates noise from rotational/translational motion of NP.
def head_group_mobility(trajin_file_object, traj_name, ligand_residue_num_list, head_group_atom_name, reference_name,
                        sulfur_atom_name='ss'):
    trajin_file_object.write('symmrmsd @%{} remap ref {}\n'.format(head_group_atom_name, reference_name))
    trajin_file_object.write('atomicfluct out rmsf.{0}.{1}.dat @%{0} bymask\n'.format(head_group_atom_name, traj_name))
    # return trajin_file_object


# Define frame to be used as reference. Default is restart file from second equilibration.
def define_reference_structure(trajin_file_object, system_name, reference_file_name='', ref_name='eq2_ref'):
    if len(reference_file_name) == 0:
        reference_file_name = '{}.equil2.rst7'.format(system_name)
    if ref_name[0] is not '[' and ref_name[-1] is not ']':
        ref_name = '[{}]'.format(ref_name)
    trajin_file_object.write('reference {} [{}]\n'.format(reference_file_name, ref_name))
    return ref_name


# Get inputs from command line.
def get_names(files_list):
    head_group_atom_name = input('Atom name for head group. Do not include @%.   ')
    au_core_residue_num = int(input('Residue number for Au core.   '))
    ligand_residue_start = int(input('First ligand residue number.   '))
    ligand_residue_stop = int(input('Last ligand residue number.   '))
    ligand_residue_num_list = list(range(ligand_residue_start, ligand_residue_stop + 1))
    trajectory_prefix = input('Prefix of trajectory...   ')
    for file_name in files_list:
        if trajectory_prefix in file_name and '.nc' in file_name:
            print('{}\n'.format(file_name))
    trajectory_start = int(input('First trajectory...   '))
    trajectory_stop = int(input('Last trajectory...   '))
    trajectory_range = list(range(trajectory_start, trajectory_stop + 1))

    return head_group_atom_name, au_core_residue_num, ligand_residue_num_list, trajectory_prefix, trajectory_range


# Grab subdirectories for use as system names.
def get_immediate_subdirectories():
    print(next(os.walk('.'))[1])
    sub_dir_list = next(os.walk('.'))[1]
    print(sub_dir_list)
    return [sub_dir for sub_dir in sub_dir_list if 'ins' not in sub_dir]


# Load .prmtop and .nc trajectories for analysis.
def load_parm_traj(trajin_file_object, parm_name, trajectories_name_list):
    trajin_file_object.write('parm {}\n'.format(parm_name))
    for trajectory_name in trajectories_name_list:
        trajin_file_object.write('trajin {}\n'.format(trajectory_name))
        # return trajin_file_object


def run_and_quit(trajin_file_object):
    trajin_file_object.write('run\n')
    trajin_file_object.write('quit\n')


def main():
    systems_list = get_immediate_subdirectories()
    print('List of systems to be analyzed: {}\n'.format('\n'.join(systems_list)))
    for system_name in systems_list:
        files_list = next(os.walk('./{}'.format(system_name)))[2]
        print('\n\nAnalyzing system: {}\n'.format(system_name))
        if 'y' not in input('Enter y to analyze this system...'):
            print('Not running this system...\n\n')
            continue
        head_group_atom_name, au_core_residue_num, ligand_residue_num_list, trajectory_prefix, trajectory_range = get_names(
            files_list)
        short_traj_name = 'md{}-{}'.format(trajectory_range[0], trajectory_range[-1])
        with open('{}/corona_traj.in'.format(system_name), 'w') as trajin_file_object:
            if len(trajectory_prefix) == 0:
                parm_name = '{}.prmtop'.format(system_name)
            else:
                parm_name = '{}.{}.prmtop'.format(trajectory_prefix, system_name)
            trajectory_names_list = ['{}.md{}.nc'.format(parm_name, trajectory_num) for trajectory_num in
                                     trajectory_range]
            load_parm_traj(trajin_file_object, parm_name, trajectory_names_list)
            ref_name = define_reference_structure(trajin_file_object, system_name)
            head_group_mobility(trajin_file_object, short_traj_name, ligand_residue_num_list, head_group_atom_name,
                                ref_name, sulfur_atom_name='ss')
            ligand_np_angle(trajin_file_object, short_traj_name, ligand_residue_num_list, au_core_residue_num,
                            head_group_atom_name, sulfur_atom_name='ss')
            run_and_quit(trajin_file_object)


if __name__ == "__main__":
    main()

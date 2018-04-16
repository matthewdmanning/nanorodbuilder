import geometry_tools
import script_writing_tools
import script_writing_object

lattice = geometry_tools.make_fcc_111([500,500,2],1.668)
file_name = '/home/mdmannin/Desktop/genetic_ligands/au_surf.xyz'
script_name = '/home/mdmannin/Desktop/genetic_ligands/au_surf.in'
leap_script = open(script_name, 'w')
#script_writing_tools.xyz_writer(file_name, lattice)
#script_writing_tools.xyz_writer(file_name, lattice)
print(len(lattice))
charge_list = [0] * len(lattice)
name_list = ['Au{}'.format(index+1) for index in list(range(len(lattice)))]
symbol_list = ['Au'] * len(lattice)
#script_writing_tools.auto_atom_creator(leap_script, lattice, 'gold', 0, symbol_list, atom_leap_name_list=name_list, atom_charge_list=name_list)
for index, coords in enumerate(lattice):
    script_writing_tools.create_atom(leap_script, 'g{}'.format(index), 'Au{}'.format(index), 'Au', atom_coords_angstrom=coords, leap_unit_alias='gold', leap_residue_number=1)

leap_script.close()


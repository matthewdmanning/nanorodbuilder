import collections
import itertools
import random
from abc import ABCMeta

import numpy as np

import geometry_tools


class NanoParticle:
    """
    Defines parameters and atoms of a single-element nanoparticle, for use in creating Leap scripts, mol2, pdb, etc.
    Subclasses are used for specific geometries.
    *** Note: all lengths are given in milli-Angstrom to avoid using float coordinates.

    Attributes:
        element *May need to include support for bimetallic, etc. NPs.
        fractional_center: fractional coordinates of NP center
        crystal_structure: 'fcc', 'bcc', 'amorphous', etc
        basis_vectors: Numpy 1-D array of the three lattice vectors
        basis_angles: Numpy 1-D array of the three angles of unit cell
        atoms: Tuple of (x,y,z) coordinates of atoms
        surface_indices: Index of atoms not fully coordinated
        core_indices: Index of (x,y,z) coordinates of atoms not on surface.
        ligand_sites: Index [atoms] of atoms bound to ligands
        ligand_types: List of names of ligands (should correspond to .mol2 file names.
        ligand_frcmods: List of frcmod locations for each ligand. *
            May want to save with each class to make outputs independent of run location.

    """

    __metaclass__ = ABCMeta

    def get_crystal_structure(self):
        fcc_elements = ['Au', 'Ag', 'Rh', 'Pd', 'Pt', 'Ni', 'Cu', 'Al', 'Pb', 'Ca', 'Sr', 'Th', 'Yb']
        if self.element in fcc_elements:
            crystal_structure = 'fcc'
        else:
            crystal_structure = 'none'
        return crystal_structure

    def get_cell_length(self):
        cell_length = {'Au': np.array([4070, 4070, 4070])}
        return cell_length[self.element]

    def __init__(self, element, fractional_center=np.array([0, 0, 0]), crystal_structure='canon', atoms=None,
                 surface_indices=None, core_indices=None, ligand_sites=None, ligand_types=None, ligand_frcmods=None):
        self.element = element
        if crystal_structure == 'canon':
            self.crystal_structure = self.get_crystal_structure()
        else:
            self.crystal_structure = crystal_structure
        self.atoms = atoms
        self.surface_indices = surface_indices
        self.core_indices = core_indices
        self.ligand_sites = ligand_sites
        self.ligand_types = ligand_types
        self.ligand_frcmods = ligand_frcmods

    def add_atom(self, coords, location):
        if location == 'shell':
            self.atoms.append(coords)
            self.surface_indices.append(len(self.atoms) - 1)
        elif location == 'core':
            self.atoms.append(coords)
            self.core_indices.append(len(self.atoms) - 1)
        else:
            self.atoms.append(coords)

    '''        
    def find_surface_atoms(self):
        return
    '''


# class NanoCluster(NanoParticle):
#    def


default_ligand_path = '/home/mdmannin/Desktop/Nanoparticles/Ligands/Charged/'
default_frcmod_path = '{}frcmod/'.format(default_ligand_path)


class CrystalSystem:
    def __init__(self, bravais_lattice=None):
        self.bravais_lattice = bravais_lattice
        if self.bravais_lattice:
            self.cell_lengths = get_cell_lengths(self.bravais_lattice)
            self.cell_angles = get_cell_angles(self.bravais_lattice)
            self.fractional_coords = get_fractional(self.bravais_lattice)


class PlaceholderType:
    def __init__(self, charge=None):
        self.charge = charge


class LigandInfo:
    def __init__(self, ligand_name, charge=None, length=None, attach_index=None):
        self.ligand_name = ligand_name
        self.charge = charge
        self.length = length
        self.attach_index = attach_index


class LigandType:
    # mol2_name should indicate the name of the ligand, without path or .mol2 extensions.
    def __init__(self, ligand_info=None, mol2_name=None, frcmod_name=None, amber_template_alias=None, charge=None,
                 length=None, attach_index=None, file_path=None):
        self.amber_template_alias = False
        self.file_path = file_path
        self.frcmod_name = frcmod_name
        if ligand_info:
            self.mol2_name = ligand_info.ligand_name
            self.charge = ligand_info.charge
            self.length = ligand_info.length
            self.attach_index = ligand_info.attach_index
            if not amber_template_alias:
                self.amber_template_alias = self.mol2_name
        if mol2_name and not self.mol2_name:
            self.mol2_name = mol2_name
        if not self.frcmod_name:
            self.frcmod_name = '{}{}.frcmod'.format(default_frcmod_path, self.mol2_name)
        if self.mol2_name and not self.amber_template_alias:
            self.amber_template_alias = self.mol2_name
        if length and not self.length:
            self.length = length
        if charge and not self.charge:
            self.charge = charge
        if attach_index and not self.attach_index:
            self.attach_index = attach_index  # Specifies index (starting with 1) of linking atom (eg. sulfur for thiolate)
        if not file_path and not self.file_path:
            self.file_path = '{}{}.mol2'.format(default_ligand_path, self.mol2_name)
        elif file_path and not self.file_path:
            self.file_path = file_path
        if amber_template_alias:
            self.amber_template_alias = 'template_{}'.format(amber_template_alias)
        elif self.mol2_name:
            self.amber_template_alias = self.mol2_name
            # print('Ligand-sulfur: {} {}'.format(self.mol2_name, self.attach_index))


# From (https://stackoverflow.com/questions/1081253/inheriting-from-instance-in-python/1081925#1081925)
_marker = object()


class Ligand:
    _inherited = ['mol2_name', 'frcmod_name', 'amber_template_alias', 'charge', 'length', 'attach_index']
    def __init__(self, parent, amber_alias=None, coords=np.array([0, 0, 0]), residue_num=np.NaN, site=None):

        self.amber_alias = amber_alias
        self.coords = coords
        self.residue_num = residue_num
        self.site = None  # Atom to which Ligand is bound.
        self.ligand_type = parent
        if self.site and self.site.ligand is not self:
            self.site.ligand = self
            # print(self)

    def impose_new_type(self, ligand_type):
        self.ligand_type = ligand_type
        return

    def __getattr__(self, name, default=_marker):
        if name in self._inherited:
            try:
                # print(getattr(self.ligand_type, name))
                return getattr(self.ligand_type, name)
            except AttributeError:
                if default is _marker:
                    raise
                return default
        if name not in self.__dict__:
            raise AttributeError(name)
        return self.__dict__[name]


# def nanoparticle_from_xyz(xyz_filename, unit_cell_indices, unit_cell_coords=None, unit_cell_lengths=None, neighbor_distance=np.NaN):


class NanoparticleSuperLattice:
    __ligand_attributes__ = ['mol2_name', 'frcmod_name', 'amber_template_alias', 'charge', 'length', 'attach_index']
    __max_num_bonds__ = 7

    def __init__(self, unit_cell_indices, unit_cell_coords=None, unit_cell_lengths=None, neighbor_distance=np.NaN,
                 amber_unit_name='nano'):
        # Must be 1x3 Numpy array.
        self.unit_cell_indices = unit_cell_indices
        self.neighbor_distance = neighbor_distance  # Nearest neighbor distance, based on crystal lattice
        self.amber_unit_name = amber_unit_name
        self.unit_cells = []
        self.unit_cell_coords = unit_cell_coords
        self.unit_cell_lengths = unit_cell_lengths
        self.atoms = []
        self.surface_atoms = []
        self.bulk_atoms = []
        self.min_bulk_coordination = 11
        self.bonded_pairs = []
        self.frcmod_name = ['AuNP.frcmod']
        # Particle Geometry
        # self.ligand_sites = []

    def add_unit_cells(self, cells_list):
        if len(self.unit_cells) > 0:
            self.unit_cells.extend(cells_list)
        elif geometry_tools.get_nested_depth(cells_list) == 2:
            for cell in cells_list:
                self.unit_cells.append(cell)
        elif geometry_tools.get_nested_depth(cells_list) == 1:
            self.unit_cells.append(cells_list)
        else:
            print('Nested list depth is not 1 or 2. Unit cells not added.')
        return

    def get_neighbor_cells(self, unit_cell):
        center_address = unit_cell.address
        neighbor_cell_list = []
        neighbor_address_list = []
        for ind in list(range(3)):
            address_range = [[num] for num in center_address]
            x = center_address[ind]
            address_range[ind] = [x - 1, x, x + 1]
            neighbor_address_list.extend(
                [address for address in itertools.product(address_range[0], address_range[1], address_range[2])])
        address_range = [[x - 1, x, x + 1] for x in center_address]
        neighbor_address_list = [address for address in
                                 itertools.product(address_range[0], address_range[1], address_range[2]) if
                                 address is not center_address]
        neighbor_set = set(neighbor_address_list)
        # print(neighbor_set)
        # print(len(neighbor_set))
        for unit_cell in self.unit_cells:
            if any((np.array(unit_cell.address) == x).all() for x in neighbor_address_list):
                neighbor_cell_list.append(unit_cell)
        unit_cell.neighbor_cells = neighbor_cell_list
        return unit_cell.neighbor_cells

    def create_unit_cell(self, address):
        new_cell = UnitCell(address, super_lattice=self)
        self.unit_cells.append(new_cell)
        return new_cell

    def create_atom(self, element, coords, unit_cell):
        np_index = len(self.atoms)
        new_atom = unit_cell.create_atom(element, coords, np_index=np_index, super_lattice=self)
        self.atoms.append(new_atom)
        return new_atom

    def get_all_neighbors(self):
        for atom in self.atoms:
            neighbors = atom.get_nearest_neighbors(self.neighbor_distance)
        return

    def get_surface_atoms(self):
        print(self.min_bulk_coordination)
        for atom in self.atoms:
            if atom.coordination_num == 0:
                atom.get_nearest_neighbors(self.neighbor_distance)
            if atom.coordination_num < self.min_bulk_coordination:
                self.surface_atoms.append(atom)
                atom.is_surface = True
            else:
                self.bulk_atoms.append(atom)
                atom.is_surface = False
        return self.surface_atoms, self.bulk_atoms

    def bond_pair(self, atom1, atom2):
        atom1.bind_atom(atom2)
        atom2.bind_atom(atom1)
        self.bonded_pairs.append([atom1, atom2])
        return

    def bond_to_neighbors(self, atom, neighbors_list=None, max_atom_bonds=7, max_neighbor_bonds=7):
        atom.nearest_neighbors.sort(key=lambda x: x.num_bonds, reverse=False)
        if not neighbors_list:
            neighbors_list = atom.nearest_neighbors
        for neighbor in neighbors_list:
            if atom.num_bonds >= max_atom_bonds:
                break
            if neighbor.num_bonds < self.__max_num_bonds__ and [neighbor, atom] not in self.bonded_pairs:
                self.bond_pair(atom, neighbor)
        return

    def bond_all_atoms(self):
        # Sort self.surface_atoms by coordination number.
        self.surface_atoms.sort(key=lambda x: x.coordination_num, reverse=False)
        for atom in self.surface_atoms:
            self.bond_to_neighbors(atom, atom.nearest_neighbors, max_atom_bonds=6, max_neighbor_bonds=6)
        for atom in [core for core in self.atoms if core not in self.surface_atoms]:
            self.bond_to_neighbors(atom, [neighbor for neighbor in atom.nearest_neighbors if
                                          neighbor not in self.surface_atoms], max_atom_bonds=7, max_neighbor_bonds=7)
        return

    def assign_atom_names(self):
        for index, atom in enumerate(self.atoms):
            if atom.element == 'Au':
                atom.amber_name = 'g{}'.format(index)
            else:
                atom.amber_name = '{}{}'.format(atom.element, index)
        return

        # def assign_ligand_names(self):
        #    for ind, ligand in enumerate(self.ligands):
        #        setattr(ligand, 'amber_alias', 'ligand{}'.format(ind))


# The GenericParticle class is a 'subclass' of NanoparticleSuperLattice, in which ligand sites are specified
class GenericParticle:
    lig_tol = 1.05

    def __init__(self, lattice_parent):
        self.parent_lattice = lattice_parent
        self.ligand_sites = []
        self.site_type_dict = {}
        self.surface_sites = []
        self.create_generic_sites()

    def __getattr__(self, item):
        inherit_list = ['dimensions', 'neighbor_distance', 'amber_unit_name', 'unit_cells', 'atoms', 'surface_atoms',
                        'bulk_atoms', 'min_bulk_coordination', 'bonded_pairs', 'frcmod_name']
        if item in inherit_list:
            return getattr(self.parent_lattice, item)
        else:
            return self.__getattribute__(item)

    def calculate_all_ligand_distances(self, potential_sites, ligand_list=None):
        if not ligand_list:
            ligand_list = self.ligand_sites
        for site in potential_sites:
            site.calculate_ligand_distances(ligand_list)

    def create_generic_sites(self):
        for atom in self.parent_lattice.surface_atoms:
            new_site = GenericSite(atom, parent_particle=self)
            self.surface_sites.append(new_site)
        return

    # This method is meant to choose which atoms are bound to ANY ligand or to all to single ligand.
    # Use place_ligands method if multiple ligand types are to be used.

    def choose_ligand_sites(self, num_ligands, potential_sites=None):
        # Clear any existing ligand sites.
        self.ligand_sites = []
        self.site_type_dict = {}
        if not potential_sites:
            potential_sites = self.surface_sites
            #print('No potential site list detected. Using surface atoms instead...')
        # print(len(potential_sites), num_ligands, len(self.site_ligand_dict))
        if len(potential_sites) < num_ligands:
            print('Number of ligands exceeds the number of potential sites. Placing ligands on all atoms.')
            # self.ligand_sites.append([site for site in potential_sites if site not in self.ligand_sites])
            for site in potential_sites:
                self.add_ligand_site(site)
                # if placeholder_type:
                #    self.site_type_dict[site] = placeholder_type
                # atom.bind_ligand(ligand_type)
            return potential_sites
        new_site = random.choice(potential_sites)
        self.add_ligand_site(new_site)
        new_site.is_ligand_site = True
        # print('Ligand sites: {}'.format(len(self.ligand_sites)))
        # for lig_index in list(range(num_ligands)):
        while len(self.ligand_sites) < num_ligands:
            potential_sites = [site for site in potential_sites if site not in self.ligand_sites]
            self.calculate_all_ligand_distances(potential_sites, self.ligand_sites)
            if len([site for site in potential_sites if site not in self.ligand_sites]) <= 0:
                break
            # print('Nanmin: {}'.format([np.nanmin(site.distances_from_ligands) for site in potential_sites if site not in self.ligand_sites]))
            # farthest_nearest = np.nanmax(
            #    [np.nanmin(site.distances_from_ligands[np.nonzero(site.distances_from_ligands)]) for site in potential_sites])
            # least_crowded_sites = [site for site in potential_sites if farthest_nearest - np.nanmin(site.distances_from_ligands[np.nonzero(site.distances_from_ligands)]) < 0.1]
            # least_crowded_sites.sort(key=lambda x: np.sum(np.power(x.distances_from_ligands[x.distances_from_ligands - 2.1 * farthest_nearest < 0.01].astype(np.float), -2)),
            #                         reverse=False)
            # crowding_index = [np.sum(np.power(site.distances_from_ligands[site.distances_from_ligands - 2.1 * farthest_nearest < 0.01].astype(np.float), -2)) for site in least_crowded_sites]
            least_crowded_sites = self.get_least_crowded_sites(potential_sites, self.ligand_sites)
            # print(crowding_index)
            # print([np.nanmin(site.distances_from_ligands) for site in least_crowded_sites])
            new_site = least_crowded_sites[0]
            self.add_ligand_site(new_site)
            # print('Chosen nearest: {}'.format(np.nanmin(new_site.distances_from_ligands[np.nonzero(new_site.distances_from_ligands)])))
            # print(np.nanmin(new_site.distances_from_ligands[np.nonzero(new_site.distances_from_ligands)]) - farthest_nearest)
        return self.ligand_sites

    def choose_unspecified_sites(self, num_total_ligands):
        self.choose_ligand_sites(num_total_ligands)

    def clear_placeholder_sites(self):
        self.site_type_dict = {}

    def choose_placeholder_sites(self, name_num_dict=None, ligand_num_list=None):
        self.clear_placeholder_sites()
        if ligand_num_list:
            name_num_dict = {}
            for num in ligand_num_list:
                new_placeholder = PlaceholderType()
                name_num_dict[new_placeholder] = num
        if not name_num_dict:
            print('Numbers of ligands not specified.')
            return None
        if not self.ligand_sites:
            total_lig_num = sum(self.ligand_sites)
            self.choose_unspecified_sites(total_lig_num)
        potential_sites = self.ligand_sites
        self.place_ligands(name_num_dict)
        return name_num_dict
        # for type, num in name_num_dict.items():
        #    new_site_list = self.choose_ligand_sites(num, placeholder_type=type, potential_sites=potential_sites)
        #    potential_sites = [site for site in potential_sites if site not in new_site_list]


    def add_ligand_site(self, site, ligand_type=None, **kwargs):
        if not site.is_ligand_site:
            site.make_ligand_site(ligand_type=ligand_type)
        self.ligand_sites.append(site)

    # Function for picking ligands (typically based on charged vs uncharged) to equally distribute ligands by type
    def place_ligands(self, type_num_dict):
        """
        :ligand_type ligand_type_list: Class LigandType
        """
        total_ligand_num = sum([type_num_dict[lig] for lig in type_num_dict])
        if 0 >= len(self.site_type_dict) > total_ligand_num:
            self.choose_ligand_sites(total_ligand_num, self.surface_sites)
        # ligand_type_by_charge = sorted(type_num_dict, key=lambda x: abs(operator.itemgetter('charge')), reverse=True)
        for ligand_type, ligand_num in type_num_dict.items():
            # print("Here's a ligand_type with {} ligands: {}\n".format(ligand_num, ligand_type))
            potential_sites = [site for site in self.ligand_sites if site not in self.site_type_dict.keys()]
            if ligand_num == len(potential_sites):
                for site in potential_sites:
                    self.add_ligand_site(site, ligand_type=ligand_type)
                    self.site_type_dict[site] = ligand_type
                return
            for ligand_index in list(range(ligand_num)):
                potential_sites = [site for site in self.ligand_sites if site not in self.site_type_dict.keys()]
                if len(potential_sites) <= 0:
                    print('No more available ligand sites. Ending assignments')
                    return
                assigned_sites = [site for site in self.site_type_dict if self.site_type_dict[site] is ligand_type]
                # print('assigned: {}\n'.format(assigned_sites))
                self.calculate_all_ligand_distances(potential_sites, assigned_sites)
                # nearest_ligand_list = [np.nanmin(site.distances_from_ligands) for site in potential_sites]
                # farthest_nearest = np.nanmax(nearest_ligand_list)
                # least_crowded_sites = [site for site in potential_sites if
                #                       np.nanmin(site.distances_from_ligands) - self.lig_tol * farthest_nearest < 0.001]
                # least_crowded_sites.sort(key=lambda x: np.sum(np.power(x.distances_from_ligands, -1)), reverse=False)
                least_crowded_sites = self.get_least_crowded_sites(potential_sites, assigned_sites, search_radius=2.1)
                next_assigned = least_crowded_sites[0]
                self.add_ligand_site(next_assigned, ligand_type=ligand_type)
                self.site_type_dict[next_assigned] = ligand_type
        return

    def get_least_crowded_sites(self, potential_sites, ligand_sites, search_radius=3.):
        # Returns a list of sites, filtered by farthest nearest ligand, and sorted by crowding function.
        farthest_nearest = np.nanmax(
            [np.nanmin(site.calculate_ligand_distances(ligand_sites)) for site in potential_sites])
        least_crowded_sites = [site for site in potential_sites if
                               farthest_nearest - np.nanmin(site.calculate_ligand_distances(ligand_sites)) < 0.1]
        least_crowded_sites.sort(key=lambda x: np.sum(np.power(x.calculate_ligand_distances(ligand_sites)[
            x.calculate_ligand_distances(ligand_sites) - search_radius * farthest_nearest < 0.01].astype(np.float),
            -2)), reverse=False)
        return least_crowded_sites


'''
    def clear_ligand_sites(self):
        for site in self.ligand_sites:
            site.clear_ligand()
        self.ligand_sites = []
        self.site_ligand_dict ={} # sites : ligands


    def update_ligand_types(self, new_type_list):
        conversion_dict = dict(zip(self.ligand_types, new_type_list))
        for site_atom, old_ligand in self.site_ligand_dict.items():
            # print(site_atom.ligand)
            try:
                # new_ligand = Ligand(conversion_dict[old_ligand.ligand_type], site=site_atom)
                # setattr(site_atom, 'ligand', new_ligand)
                self.add_ligand_site(site_atom, ligand_type=conversion_dict[old_ligand.ligand_type])
                # print(site_atom.ligand.ligand_type)
            except AttributeError or KeyError:
                print('No ligand assigned to {}'.format(site_atom))
                # print(site_atom.ligand)
        # Change attributes of all individual ligands.
        # for ligand, old_type in self.ligands_dict.items():
        #    for attribute in self.__ligand_attributes__:
        #        setattr(ligand, attribute, getattr(conversion_dict[old_type], attribute))
        #    self.ligands_dict[ligand] = conversion_dict[old_type]
        # Change values for all LigandType instances.
        # print(self.ligand_types)
        # self.ligand_types = [conversion_dict[old_type] for old_type in self.ligand_types]
        # new_type_list = [conversion_dict[old_type] for old_type in self.ligand_types]
        self.ligand_types = new_type_list
        print('New ligand_type list: {}\n'.format(self.ligand_types))
        return
'''


class Nanoparticle:
    def __init__(self, generic_parent):
        self.generic_parent = generic_parent
        self.site_ligand_dict = {}
        self.ligand_types = []

    def __getattr__(self, item):
        inherit_list = ['dimensions', 'neighbor_distance', 'amber_unit_name', 'unit_cells', 'atoms', 'surface_atoms',
                        'bulk_atoms', 'min_bulk_coordination', 'bonded_pairs', 'frcmod_name', 'surface_sites',
                        'ligand_sites', 'site_type_dict']
        if item in inherit_list:
            return getattr(self.generic_parent, item)
        else:
            return self.__getattribute__(item)

    def create_ligand_type(self, ligand_info=None, mol2_name=None, frcmod_name=None, amber_template_alias=None,
                           charge=0, length=np.NaN,
                           attach_index=1):
        new_ligand_type = LigandType(ligand_info=ligand_info, mol2_name=mol2_name, frcmod_name=frcmod_name,
                                     amber_template_alias=amber_template_alias, charge=charge, length=length,
                                     attach_index=attach_index)
        self.ligand_types.append(new_ligand_type)
        return new_ligand_type

    def create_type_dict(self, ligand_info_list=None, ligand_name_list=None, **kwargs):
        ligand_type_dict = {}
        generic_type_set = set(val for val in self.generic_parent.site_type_dict.values())
        # print('Generic set: {}'.format(generic_type_set))
        if ligand_info_list:
            for new_type_info, generic_type in zip(ligand_info_list, generic_type_set):
                new_type = self.create_ligand_type(ligand_info=new_type_info, **kwargs)
                ligand_type_dict[generic_type] = new_type
            return ligand_type_dict
        elif ligand_name_list:
            for new_type_name, generic_type in zip(ligand_name_list, generic_type_set):
                new_type = self.create_ligand_type(mol2_name=new_type_name, **kwargs)
                ligand_type_dict[generic_type] = new_type
            return ligand_type_dict

    def create_ligand_sites(self, ligand_info_num_dict=None, ligand_info_list=None, ligand_type_dict=None,
                            ligand_name_list=None, **kwargs):
        if not ligand_type_dict:
            ligand_type_dict = {}
        if ligand_info_num_dict:
            type_list = [type for type in self.generic_parent.site_type_dict.values()]
            ligand_counter = collections.Counter(type_list)
            print(ligand_counter)
            for ligand_info, ligand_num in ligand_info_num_dict.items():
                print(ligand_info, ligand_num)
                for generic_type, generic_num in ligand_counter.items():
                    print(generic_type, generic_num)
                    if ligand_num == generic_num and generic_type not in ligand_type_dict.keys():
                        print('Found a match')
                        new_type = self.create_ligand_type(ligand_info=ligand_info, **kwargs)
                        ligand_type_dict[generic_type] = new_type
                        break
        elif ligand_name_list:
            ligand_type_dict = self.create_type_dict(ligand_name_list=ligand_name_list, **kwargs)
        elif ligand_info_list:
            ligand_type_dict = self.create_type_dict(ligand_info_list=ligand_info_list, **kwargs)
        if ligand_type_dict:
            print(ligand_type_dict.keys())
            for generic_site, generic_type in self.generic_parent.site_type_dict.items():
                new_type = ligand_type_dict[self.generic_parent.site_type_dict[generic_site]]
                new_ligand = Ligand(new_type)
                new_site = LigandSite(generic_site, ligand=new_ligand)
                new_ligand.site = new_site
                self.site_ligand_dict[new_site] = new_ligand
        else:
            print('No ligand types specified. Ligands not placed.')
            return None


class UnitCell:
    def __init__(self, address, super_lattice=None):
        self.address = address
        self.atoms = []
        self.neighbor_cells = []
        self.super_lattice = super_lattice

    def add_atom(self, atom):
        self.atoms.append(atom)
        return

    def create_atom(self, element, coords, np_index=np.NaN, super_lattice=None):
        cell_index = len(self.atoms)
        new_atom = NanoparticleAtom(element, coords, np_index, unit_cell_index=cell_index, unit_cell=self,
                                    super_lattice=super_lattice)
        self.atoms.append(new_atom)
        return new_atom

    def get_neighbor_cells(self):
        if self.neighbor_cells:
            return self.neighbor_cells
        else:
            self.neighbor_cells = self.super_lattice.get_neighbor_cells(self)
            return self.neighbor_cells


class NanoparticleAtom:
    def __init__(self, element, coords_int, amber_name=None, unit_cell_index=np.NaN, unit_cell=None, super_lattice=None,
                 charge=0):
        self.element = element
        self.charge = charge
        self.coords = coords_int
        if type(self.coords) is list:
            self.coords = np.array(self.coords)
        self.super_lattice = super_lattice  # Class NanoparticleSuperLattice
        self.amber_name = amber_name
        self.unit_cell_index = unit_cell_index
        self.unit_cell = unit_cell  # Class UnitCell
        self.nearest_neighbors = []
        self.num_bonds = 0
        self.bonded_atoms = []
        self.coordination_num = 0
        self.is_surface = False

    def assign_superlattice(self, super_lattice):
        self.super_lattice = super_lattice
        return

    def assign_unit_cell(self, unit_cell):
        self.unit_cell = unit_cell
        return

    def get_nearest_neighbors(self, neighbor_distance):
        if self.super_lattice and self.unit_cell:
            neighbor_cells = self.unit_cell.get_neighbor_cells()
            neighbors = [atom for atom in self.unit_cell.atoms if atom is not self]
            # = self.super_lattice.get_neighbor_cells(self.unit_cell.address)
            for unit_cell in neighbor_cells:
                atoms = [atom for atom in unit_cell.atoms if atom is not self]
                for atom in atoms:
                    if geometry_tools.calculate_distance(self.coords, atom.coords) - neighbor_distance < 0.001:
                        neighbors.append(atom)
            self.nearest_neighbors = list(set(neighbors))
            self.coordination_num = len(self.nearest_neighbors)
            # print(neighbors)
            # print('Num neighbors: {}'.format(self.coordination_num))
            return neighbors
        else:
            print('Either self.unit_cell or self.super_lattice is not defined.')
            return []

    def bind_atom(self, atom):
        if atom not in self.bonded_atoms:
            self.bonded_atoms.append(atom)
        self.update_num_bonds()
        return

    def update_num_bonds(self):
        self.num_bonds = len(self.bonded_atoms)
        # if self.is_ligand_site:
        #    self.num_bonds += 1
        # return


'''
    def clear_ligand(self):
        self.is_ligand_site = False
        self.ligand = None
        self.update_num_bonds()
        
   
             
    def bind_ligand(self, ligand):
        self.is_ligand_site = True
        self.ligand = ligand
        self.update_num_bonds()
        return'''

'''
sup = NanoparticleSuperLattice([1,1,10])
atom_list = []
for i in list(range(10)):
    new_cell = sup.create_unit_cell([0,0,i])
    for j in list(range(4)):
        new_atom = sup.create_atom('Au', [j,j,i], new_cell)
print(sup.unit_cells)
print(sup.atoms)
for cell in sup.unit_cells:
    atoms = cell.atoms
    for atom in atoms:
        neighbors = atom.get_nearest_neighbors()
        print(neighbors)
        for neigh in neighbors:
            print(neigh.coords)

'''


# Composition of NanoparticleAtom for a surface atom within a GenericParticle object.
# Created for all surface atoms which may be bound to a ligand.
# Attribute ligand_type is a placeholder for specifying placeholder type #1, 2, etc.
# Realization should be made by creating instance of LigandSite, where a specific Ligand and LigandType are indicated.
class GenericSite:
    def __init__(self, parent_atom, parent_particle=None, charge=None, ligand_type=None):
        self.parent_atom = parent_atom
        self.parent_particle = parent_particle
        self.is_ligand_site = False
        self.distances_from_ligands = []
        self.site_distance_dict = {}
        self.ligand_type = ligand_type

    def __getattr__(self, item):
        inherit_list = ['coords', 'element', 'super_lattice', 'amber_name', 'unit_cell_index', 'unit_cell',
                        'nearest_neighbors', 'num_bonds', 'bonded_atoms', 'coordination_num']
        if item in inherit_list:
            return getattr(self.parent_atom, item)
        else:
            return self.__getattribute__(item)

    def calculate_ligand_distances(self, site_list=None):
        if not site_list:
            site_list = self.parent_particle.ligand_sites
        new_sites_list = [site for site in site_list if site not in self.site_distance_dict and site is not self]
        for site in new_sites_list:
            self.site_distance_dict[site] = geometry_tools.calculate_distance(self.coords, site.coords)
        distance_array = np.array(
            [self.site_distance_dict[site] for site in self.site_distance_dict if site in site_list])
        # Should I exclude atoms not in passed site_list?
        self.distances_from_ligands = distance_array[np.nonzero(distance_array)]
        return distance_array

    def make_ligand_site(self, ligand_type=None):
        self.is_ligand_site = True
        self.ligand_type = ligand_type


# Composition of GenericSite for atom bound to a specific Ligand object, with a specific LigandType.
# Only created for surface atom bound to an unspecified ligands (ie. where parent_site.is_ligand_site == True)
# Contained within a Nanoparticle object.
class LigandSite:
    def __init__(self, parent_site, ligand=None):
        self.parent_site = parent_site
        self.ligand = ligand

    def __getattr__(self, item):
        inherit_list = ['coords', 'element', 'super_lattice', 'amber_name', 'unit_cell_index', 'unit_cell',
                        'nearest_neighbors', 'num_bonds', 'bonded_atoms', 'coordination_num', 'parent_atom',
                        'parent_particle']
        if item in inherit_list:
            return getattr(self.parent_site, item)
        else:
            return self.__getattribute__(item)

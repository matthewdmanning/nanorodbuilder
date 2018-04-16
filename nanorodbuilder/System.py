'''
Creates a class that handles building systems in tLeap via an input script.
'''

### Creates class for tLeap system.
def merge_systems(systems_list, amber_name=None):
    combined_system = LeapSystem(amber_unit_name=amber_name)
    for system in systems_list:
        for residue in system.residue_objects:
            combined_system.residue_objects.append(residue)
        for pair in system.bond_pairs_list:
            if type(pair) is list and len(pair) == 2:
                combined_system.bond_pairs_list.append(pair)


class LeapSystem:
    def __init__(self, system_name=None, leap_filname=None, amber_unit_name=None):
        self.system_name = system_name
        self.leap_filename = leap_filname
        self.amber_unit_name = amber_unit_name
        self.amber_unit_num = None
        self.bond_pairs_list = []
        self.residue_objects = {}  # residue object : amber_alias
        self.input_file = None

    def open_for_writing(self):
        self.input_file = open(self.leap_filename, 'w')
        return self.input_file

    def close_input_file(self):
        self.input_file.close()
        return

    def write(self, in_line):
        self.input_file.write(in_line)
        return

    def check_bond_duplicate(self, atom1, atom2):
        if [atom1, atom2] in self.bond_pairs_list or [atom2, atom1] in self.bond_pairs_list:
            return True
        else:
            return False

    def bond_atoms(self, atom1, atom2, bond_order=1):
        if not self.check_bond_duplicate(atom1, atom2):
            self.bond_pairs_list.append([atom1, atom2])
            return True
        else:
            return False

    def create_residue(self, residue, amber_alias=None):
        if not amber_alias and hasattr(residue, 'amber_alias'):
            amber_alias = getattr(residue, 'amber_alias')
        self.residue_objects[residue] = amber_alias
        residue_num = len(self.residue_objects)
        return residue_num

    def get_frcmods(self):
        frcmod_list = []
        for residue in self.residue_objects:
            try:
                frcmod_name = residue.frcmod_name
                if frcmod_name not in frcmod_list:
                    frcmod_list.append(frcmod_name)
            except:
                continue
        if 'AuNP.frcmod' not in frcmod_list:
            frcmod_list.append('AuNP.frcmod')
        #frcmod_list.append('VACMIN_ONLY_AuNP.frcmod')
        # print(frcmod_list)
        return frcmod_list

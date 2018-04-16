class Ligand():
    '''
    Class for construction of ligands for use with NanoParticle class.

    Attributes:
        Name: should match mol2 name *Not all atoms should be in same order as mol2.
        frcmod: for use with AMBER
        Attach index: index of open shell attaching atom (eg. sulfur in thiolate ligand)
        Charge: indicates cationic, anionic, or neutral
        element_list: lists elements of each atom
        atom_type_list: list of each AMBER atom ligand_type, should match frcmod file
        ligand_length: approximate length of ligand (twice hydrodynamic radius [can measure in DS from geom. optimized?])

    '''

    def __init__(self, name):
        self.name = name

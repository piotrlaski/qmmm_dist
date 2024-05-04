import re
import os

input1 = r'C:\Users\piotr\Documents\working_dirs_lapek\oniom_dist\example_files\cam_6311_lan_S1_qmmm.log'

class Atom:
    def __init__(self, x:float, y:float, z:float, name:str=None,  at_id:int=None) -> None:
        self.xyz = (float(x),float(y),float(z))
        if at_id is not None: 
            self.at_id = int(at_id)
        else:
            self.at_id = None
        if name is not None: 
            self.name = str(at_id)
        else:
            self.name = None        
    def set_id(self, new_id:int) -> None:
        self.at_id = int(new_id)
    def set_name(self, new_name:str) -> None:
        self.name = str(new_name)

class Molecule:
    def __init__(self, atoms:list[Atom],  mol_id:int=None) -> None:
        self.coordinate_type = None
        self.atoms = {atom.at_id : atom for atom in atoms}
        if mol_id is not None: 
            self.mol_id = int(mol_id)
        else:
            self.mol_id = None
    def set_id(self, new_id:int) -> None:
        self.mol_id = int(new_id)
    def get_atom_list(self) -> list[Atom]:
        return [atom for atom in list(self.atoms.values())]
    def map_name_id(self) -> dict:
        return {atom.at_id : atom.name for atom in self.get_atom_list()}
    def set_coordinate_type(self, coordinate_type: str) -> None:
        '''raw or shifted'''
        self.coordinate_type = coordinate_type
    def shifted_to_raw(self, shift_vec:tuple[float, float, float]) -> 'Molecule':
        '''returns a copy of the molecule in raw coords'''
        if self.coordinate_type != "shifted":
            raise ValueError
        atom_list_raw = []
        for atom in self.get_atom_list():
            atom_raw_xyz = tuple(at_sh_coord + sh_vec_coord for at_sh_coord, sh_vec_coord in zip (atom.xyz, shift_vec))
            atom_raw = Atom(*atom_raw_xyz, at_id=atom.at_id)
            atom_list_raw.append(atom_raw)
        mol_raw = Molecule (atoms=atom_list_raw, mol_id=self.mol_id)
        mol_raw.set_coordinate_type('raw')
        return mol_raw

class Optimization:
    def __init__(self) -> None:
        # molecules as written in input
        self.initial_mol_raw = None
        self.intermediate_mols_raw = []
        self.final_mol_raw = None
        # molecules shifted by a vector by Gaussian
        self.initial_mol_shifted = None
        self.intermediate_mols_shifted = []
        self.final_mol_shifted = None
        self.shift_vec = None
    def calc_shift_vec(self) -> None:
        mol_raw = self.initial_mol_raw
        mol_shifted = self.initial_mol_shifted
        if not mol_raw and mol_shifted:
            raise ValueError
        atom1_raw = mol_raw.atoms[1]
        atom1_shifted = mol_shifted.atoms[1]
        shift_vec = tuple(coord_raw - coord_shifted for coord_raw, coord_shifted in zip(atom1_raw.xyz, atom1_shifted.xyz))
        self.shift_vec = shift_vec
    def calc_mols_raw(self) -> None:
        if not self.shift_vec:
            raise ValueError
        self.intermediate_mols_raw = [mol_sh.shifted_to_raw(self.shift_vec) for mol_sh in self.intermediate_mols_shifted]
        self.final_mol_raw = self.final_mol_shifted.shifted_to_raw(self.shift_vec)
def parse_intermediate_atomic_line(line:str) -> Atom:
    line = line.strip(r'\n')
    at_id = line.split()[0]
    atom_xyz = tuple(line.split()[3:6])
    return Atom(*atom_xyz, at_id=at_id)

def parse_initial_atomic_line(line:str) -> Atom:
    line = line.strip(r'\n')
    #check if QMMM or isolOpt:
    if re.findall('--', line):
        at_name = line.split('--')[0]
        atom_xyz = line.split()[2:5]
    else:
        at_name = line.split()[0]
        atom_xyz = line.split()[1:4]
    return Atom(*atom_xyz, name=at_name)

def get_initial_molecule(logfile:os.PathLike) -> Molecule:
    initial_atomic_line = []
    with open (input1, mode='r') as f:
        for line in f:
            #set iterator is at the start of initial atomic geometry
            if re.findall(r'Symbolic Z-matrix:', line):
                cur_line = next(f)
                while re.findall(r'Charge', cur_line):
                    cur_line = next(f)
            else: continue
            #read initial atomic geometry until end of high-level molecule/end of opt part
            while re.findall(r'H\s\n', cur_line) or re.findall(r'[0-9]\s\n', cur_line):
                initial_atomic_line.append(cur_line)
                cur_line = next(f)
            break
    #parse initial molecule atom list and add their ids based on their order (starts at id = 1)
    parsed_initial_atoms = [parse_initial_atomic_line(line) for line in initial_atomic_line]
    for i, atom in enumerate(parsed_initial_atoms):
        atom.set_id(i+1)
    initial_molecule = Molecule (atoms=parsed_initial_atoms, mol_id=1)
    initial_molecule.set_coordinate_type('raw')
    return (initial_molecule)

def get_next_intermediate_molecule(file, id_name_mapping: dict) -> Molecule | int:
        for line in file:
            #set iterator is at the start of next intermediate atomic geometry, skipping some lines
            if re.findall(r'Standard orientation:', line):
                cur_line = next(file)
                cur_line = next(file)
                cur_line = next(file)
                cur_line = next(file)
                cur_line = next(file)
            else: continue
            #read intermediate atomic geometry based on the number of atoms in the initial molecule
            intermediate_atomic_line = []
            for i in range(len(id_name_mapping)):
                intermediate_atomic_line.append(cur_line)
                cur_line = next(file)
            parsed_intermediate_atoms = [parse_intermediate_atomic_line(at_line) for at_line in intermediate_atomic_line]
            #name the atoms correctly based on the provided mapping
            for atom in parsed_intermediate_atoms:
                atom.set_name(id_name_mapping[atom.at_id])
            intermediate_molecule = Molecule (atoms=parsed_intermediate_atoms)
            intermediate_molecule.set_coordinate_type('shifted')
            return (intermediate_molecule)

def get_all_molecules(logfile:os.PathLike) -> tuple[Molecule, list[Molecule]]:
    '''returns (inital molecule in raw coords, all molecules in Gaussian coordinates)'''
    initial_mol = get_initial_molecule(logfile)
    id_name_mapping = (initial_mol.map_name_id())
    molecule_list = []
    with open (input1, mode='r') as f:
        while cur_molecule := get_next_intermediate_molecule(f, id_name_mapping):
            # add a molecule ID
            cur_molecule.set_id(len(molecule_list) + 1)
            molecule_list.append(cur_molecule)
    return (initial_mol, molecule_list)

if __name__ == '__main__':
    inital_mol_raw, all_mols_shifted = get_all_molecules(input1)
    calc = Optimization()
    calc.initial_mol_raw = inital_mol_raw
    calc.initial_mol_shifted = all_mols_shifted[0]
    calc.intermediate_mols_shifted = all_mols_shifted[1:-1]
    calc.final_mol_shifted = all_mols_shifted[-1]
    calc.calc_shift_vec()
    calc.calc_mols_raw()
    pass
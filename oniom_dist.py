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
            return (intermediate_molecule)
        #return -1 if EOF
        # else: raise EOFError

if __name__ == '__main__':
    initial_mol = get_initial_molecule(input1)
    id_name_mapping = (initial_mol.map_name_id())
    molecule_list = [initial_mol]
    with open (input1, mode='r') as f:
        while cur_molecule := get_next_intermediate_molecule(f, id_name_mapping):
            # add a molecule ID
            cur_molecule.set_id(len(molecule_list) + 1)
            molecule_list.append(cur_molecule)

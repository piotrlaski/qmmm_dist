import re
import os

input1 = r'C:\Users\piotr\Documents\working_dirs_lapek\oniom_dist\example_files\cam_6311_lan_T1_qmmm.log'

class Atom:
    def __init__(self, x:float, y:float, z:float, name:str=None,  at_id:int=None) -> None:
        self.at_id = at_id
        self.name = name
        self.xyz = (x,y,z)
    def set_id(self, new_id:int) -> None:
        self.at_id = new_id
    def set_name(self, new_name:str) -> None:
        self.name = new_name

class Molecule:
    def __init__(self, atoms:list[Atom],  mol_id:int=None) -> None:
        self.mol_id = mol_id
        self.atoms = {atom.at_id : atom for atom in atoms}
    def set_id(self, new_id:int) -> None:
        self.mol_id = new_id
    def get_atom_list(self) -> list[Atom]:
        return [atom for atom in list(self.atoms.values())]
    def map_name_id(self) -> dict:
        return {atom.at_id : atom.name for atom in self.get_atom_list()}

def parse_intermediate_atomic_line(line:str) -> Atom:
    line = line.strip(r'\n')
    at_id = line.split()[0]
    atom_xyz = tuple(line.split()[3:6])
    return Atom(*atom_xyz, name=at_id)

def parse_initial_atomic_line(line:str) -> Atom:
    line = line.strip(r'\n')
    at_name = line.split('--')[0]
    atom_xyz = line.split()[2:5]
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
            #read initial atomic geometry until end of high-level molecule
            while re.findall(r'H\s\n', cur_line):
                initial_atomic_line.append(cur_line)
                cur_line = next(f)
            break
    #parse initial molecule atom list and add their ids based on their order (starts at id = 1)
    parsed_initial_atoms = [parse_initial_atomic_line(line) for line in initial_atomic_line]
    for i, atom in enumerate(parsed_initial_atoms):
        atom.set_id(i+1)
    initial_molecule = Molecule (atoms=parsed_initial_atoms, mol_id=1)
    return (initial_molecule)


if __name__ == '__main__':
    initial_mol = get_initial_molecule(input1)
    print(initial_mol.atoms[1].xyz)
    id_name_mapping = (initial_mol.map_name_id())
    print (id_name_mapping[1])
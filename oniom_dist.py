import re
import os
from math import sqrt, acos, degrees
import numpy as np
import matplotlib.pyplot as plt

class Atom:
    def __init__(self, x:float, y:float, z:float, name:str=None,  at_id:int=None) -> None:
        self.xyz = (float(x),float(y),float(z))
        if at_id is not None: 
            self.at_id = int(at_id)
        else:
            self.at_id = None
        if name is not None: 
            self.name = str(name)
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
            atom_raw = Atom(*atom_raw_xyz, name=atom.name, at_id=atom.at_id)
            atom_list_raw.append(atom_raw)
        mol_raw = Molecule (atoms=atom_list_raw, mol_id=self.mol_id)
        mol_raw.set_coordinate_type('raw')
        return mol_raw
    def calc_distance(self, at_id1:int, at_id2:int) -> float:
        dist = distance(self.atoms[at_id1].xyz, self.atoms[at_id2].xyz)
        return (dist)
    def calc_angle(self, at_id1:int, at_id2:int, at_id3:int) -> float:
        ang = degrees(angle(self.atoms[at_id1].xyz, self.atoms[at_id2].xyz, self.atoms[at_id3].xyz))
        return (ang)
    
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
        self.map_name_id = None
    def calc_shift_vec(self) -> None:
        mol_raw = self.initial_mol_raw
        mol_shifted = self.initial_mol_shifted
        if not mol_raw and mol_shifted:
            raise ValueError
        #calculate the general Gaussian to raw vector based on the first atom shift
        atom1_raw = mol_raw.atoms[1]
        atom1_shifted = mol_shifted.atoms[1]
        shift_vec = tuple(coord_raw - coord_shifted for coord_raw, coord_shifted in zip(atom1_raw.xyz, atom1_shifted.xyz))
        self.shift_vec = shift_vec
    def calc_mols_raw(self) -> None:
        if not self.shift_vec:
            raise ValueError
        self.intermediate_mols_raw = [mol_sh.shifted_to_raw(self.shift_vec) for mol_sh in self.intermediate_mols_shifted]
        self.final_mol_raw = self.final_mol_shifted.shifted_to_raw(self.shift_vec)
    def set_map_name_id(self, map_name_id) -> None:
        self.map_name_id = map_name_id
    def get_distance_list(self, at_id1:int, at_id2:int) -> list[float]:
        all_mols = self.get_mol_list()
        dists = [mol.calc_distance(at_id1, at_id2) for mol in all_mols]
        return (dists)
    def get_angle_list(self, at_id1:int, at_id2:int, at_id3:int) -> list[float]:
        all_mols = self.get_mol_list()
        angs = [mol.calc_angle(at_id1, at_id2, at_id3) for mol in all_mols]
        return (angs)
    def get_mol_list(self) -> list[Molecule]:
        all_mols = self.intermediate_mols_shifted
        all_mols.insert(0, self.initial_mol_shifted)
        all_mols.append(self.final_mol_shifted)
        return (all_mols)
    def save_xyz_at_step(self, step_nb:int, output_path:os.PathLike) -> None:
        nb_atoms = len(self.initial_mol_raw.get_atom_list())
        mol = self.get_mol_list()[step_nb]
        with open(output_path, '+w') as f:
            f.write(f'{str(nb_atoms)}\n')
            for atom in mol.get_atom_list():
                f.write(f'{atom.name} {atom.xyz[0]:.6f} {atom.xyz[1]:.6f} {atom.xyz[2]:.6f}\n')
    
def parse_intermediate_atomic_line(line:str) -> Atom:
    line = line.strip(r'\n')
    at_id = line.split()[0]
    atom_xyz = tuple(line.split()[3:6])
    return Atom(*atom_xyz, at_id=at_id)

def parse_initial_atomic_line(line:str) -> Atom:
    line = line.strip(r'\n')
    #check if QMMM or isolOpt:
    if re.findall('--', line):
        at_name = line.split('--')[0].split()[0]
        atom_xyz = line.split()[2:5]
    else:
        at_name = line.split()[0]
        atom_xyz = line.split()[1:4]
    return Atom(*atom_xyz, name=at_name)

def get_initial_molecule(logfile:os.PathLike) -> Molecule:
    initial_atomic_line = []
    with open (logfile, mode='r') as f:
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
        atom.set_id(i + 1)
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
    with open (logfile, mode='r') as f:
        while cur_molecule := get_next_intermediate_molecule(f, id_name_mapping):
            # add a molecule ID
            cur_molecule.set_id(len(molecule_list) + 1)
            molecule_list.append(cur_molecule)
    return (initial_mol, molecule_list)

def get_calc(logfile:os.PathLike) -> Optimization:
    '''fully parses all info in a logfile into a complete Optimiziation instance'''
    inital_mol_raw, all_mols_shifted = get_all_molecules(logfile)
    calc = Optimization()
    calc.initial_mol_raw = inital_mol_raw
    calc.initial_mol_shifted = all_mols_shifted[0]
    calc.intermediate_mols_shifted = all_mols_shifted[1:-1]
    calc.final_mol_shifted = all_mols_shifted[-1]
    calc.calc_shift_vec()
    calc.calc_mols_raw()
    calc.set_map_name_id(inital_mol_raw.map_name_id())
    return (calc)

def distance(r1:tuple[float,float,float], r2:tuple[float,float,float]) -> float:
    return (sqrt((r1[0] - r2[0])**2 + (r1[1] - r2[1])**2 + (r1[2] - r2[2])**2))

def vector_from_points(r1:tuple[float,float,float], r2:tuple[float,float,float]) -> float:
    return (r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2])

def dot_product(v1:tuple[float,float,float], v2:tuple[float,float,float]) -> float:
    return (sum((a*b) for a, b in zip(v1, v2)))

def vec_length(v:tuple[float,float,float]) -> float:
    return sqrt(dot_product(v, v))

def angle(r1:tuple[float,float,float], r2:tuple[float,float,float], r3:tuple[float,float,float]) -> float:
    v1 = vector_from_points(r2, r1)
    v2 = vector_from_points(r2, r3)
    return acos(dot_product(v1, v2) / (vec_length(v1) * vec_length(v2)))

if __name__ == '__main__':
    input_s0 = r'C:\Users\piotr\Documents\VS_Code\qmmm_dist\files\cam_6311_lan_S0_qmmm.log'
    input_s1 = r'C:\Users\piotr\Documents\VS_Code\qmmm_dist\files\cam_6311_lan_S1_qmmm.log'
    input_t1 = r'C:\Users\piotr\Documents\VS_Code\qmmm_dist\files\cam_6311_lan_T1_qmmm.log'
    calc_s0  = get_calc(input_s0)
    calc_s1  = get_calc(input_s1)
    calc_t1  = get_calc(input_t1)
    calc_s0.save_xyz_at_step(output_path='s0_final.xyz', step_nb=-1)
    calc_s1.save_xyz_at_step(output_path='s1_final.xyz', step_nb=-3)
    calc_t1.save_xyz_at_step(output_path='t1_final.xyz', step_nb=-1)



    # x1 = range (1, len(d1) + 1)
    # plt.plot(x1, d1, marker='.', markersize='2', linewidth=0.5, color='black')
    # plt.ylim(3.3, 3.6)
    # plt.axhline(y=d1[-1], color='r', linestyle=':', label=f'{d1[-1]:.3f}')
    # plt.title(r'S0 (CAM-B3lYP/6-31G**+Lanl2dz)')
    # plt.grid(which='major', linestyle='-.', color='grey', alpha=0.3)
    # plt.legend()
    # plt.xlabel(r'Iteration number')
    # plt.ylabel(r'Rh...Rh distance / A')
    # plt.savefig('S0.png', dpi = 600)
    # plt.close()

    # x2 = range (1, len(d2) + 1)
    # plt.plot(x2, d2, marker='.', markersize='2', linewidth=0.5, color='black')
    # plt.ylim(3.0, 3.5)
    # plt.axhline(y=d2[-1], color='r', linestyle=':', label=f'{d2[-1]:.3f}')
    # plt.title(r'S1 (CAM-B3lYP/6-31G**+Lanl2dz)')
    # plt.grid(which='major', linestyle='-.', color='grey', alpha=0.3)
    # plt.legend()
    # plt.xlabel(r'Iteration number')
    # plt.ylabel(r'Rh...Rh distance / A')
    # plt.savefig('S1.png', dpi = 600)
    # plt.close()

    # x3 = range (1, len(d3) + 1)
    # plt.plot(x3, d3, marker='.', markersize='2', linewidth=0.5, color='black')
    # plt.ylim(3.1, 3.4)
    # plt.axhline(y=d3[-1], color='r', linestyle=':', label=f'{d3[-1]:.3f}')
    # plt.title(r'T1 (CAM-B3lYP/6-31G**+Lanl2dz)')
    # plt.grid(which='major', linestyle='-.', color='grey', alpha=0.3)
    # plt.legend()
    # plt.xlabel(r'Iteration number')
    # plt.ylabel(r'Rh...Rh distance / A')
    # plt.savefig('T1.png', dpi = 600)
    # plt.close()
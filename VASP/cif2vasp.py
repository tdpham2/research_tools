from ase.io import read, write
import sys
from ase.build import sort
import subprocess
import os

pbe_path = '/home/tdp797/package/pseudopotentials/PBE/potpaw_PBE.54'

pbe_path = 'path/to/pbe/pseudopotentials' # E.g '/home/pseudopotentials/PBE/potpaw_PBE.54'

# Pseudopotentials that I used as suggested by VASP
pbe = {'H':'H', 'He':'He', 'Li':'Li_sv', 'Be':'Be', 'B':'B', 'C':'C', 'N':'N', 'O':'O', 'F':'F', 'Ne':'Ne',
        'Na':'Na_pv', 'Mg':'Mg', 'Al':'Al', 'Si':'Si', 'P':'P', 'S':'S', 'Cl':'Cl', 'Ar':'Ar',
        'K':'K_sv', 'Ca':'Ca_sv', 'Sc': 'Sc_sv', 'Ti':'Ti_sv', 'V':'V_sv', 'Cr': 'Cr_pv', 'Mn':'Mn_pv', 'Fe':'Fe', 'Co':'Co', 'Ni':'Ni', 'Cu':'Cu', 'Zn':'Zn', 'Ga':'Ga_d', 'Ge':'Ge_d', 'As':'As', 'Se':'Se', 'Br':'Br', 'Kr':'Kr',
        'Rb':'Rb_sv', 'Sr':'Sr_sv', 'Y':'Y_sv', 'Zr':'Zr_sv', 'Nb':'Nb_sv', 'Mo':'Mo_sv', 'Tc':'Tc_pv', 'Ru':'Ru_pv', 'Rh':'Rh_pv', 'Pd':'Pd', 'Ag':'Ag', 'Cd':'Cd', 'In':'In_d', 'Sn':'Sn_d', 'Sb':'Sb', 'Te':'Te', 'I': 'I', 'Xe':'Xe',
        'Cs':'Cs_sv', 'Ba':'Ba_sv', 
        'La':'La', 'Ce':'Ce', 'Pr':'Pr_3', 'Nd':'Nd_3', 'Pm':'Pm_3', 'Sm':'Sm_3', 'Eu':'Eu_2', 'Gd':'Gd_3', 'Tb':'Tb_3', 'Dy': 'Dy_3', 'Ho': 'Ho_3', 'Er': 'Er_3', 'Tm':'Tm_3', 'Yb':'Yb_2'
        }

#TODO: Finish the dictionary

pbe_valency = {'H':1, 'He':2, 'Li':3, 'Be':2, 'B':3, 'C':4, 'N':5, 'O':6, 'F':7, 'Ne':8,
        'Na':7, 'Mg':2, 'Al':3, 'Si':4, 'P':5, 'S':6, 'Cl':7, 'Ar':8,
        'K':9, 'Ca':10, 'Sc': 11, 'Ti':12, 'V':13, 'Cr':12, 'Mn':13, 'Fe':8, 'Co':9, 'Ni':10, 'Cu':11, 'Zn':12, 'Ga':13, 'Ge':14, 'As':5, 'Se':6, 'Br':7, 'Kr':8,
        'Rb':9, 'Sr':10, 'Y':11, 'Zr':12, 'Nb':13, 'Mo':14, 'Tc':13, 'Ru':14, 'Rh':15, 'Pd':10, 'Ag':11, 'Cd':12, 'In':13, 'Sn':14, 'Sb':5, 'Te':6, 'I': 7, 'Xe':8,
        'Cs':9, 'Ba':10, 
        'La':11, 'Ce':12, 'Pr':11, 'Nd':11, 'Pr':11, 'Nd':11, 'Pm':11, 'Sm':11, 'Eu': 8, 'Gd':9, 'Tb': 9, 'Dy': 9, 'Ho': 9, 'Er': 9, 'Tm':9, 'Yb':8
        }

atoms = read(sys.argv[1])
atoms_sorted = sort(atoms)

write('POSCAR', atoms_sorted)
elements = subprocess.check_output('cat POSCAR | head -1', shell=True).decode('ascii').strip().split()

if os.path.isfile('POTCAR'):
    pass
else:
    for element in elements:
        if element in pbe:
            path = pbe[element]
            subprocess.run("cat {}/{}/POTCAR >> POTCAR".format(pbe_path, path), shell=True)
        else:
            print("Missing POTCAR for {}".format(element))
            subprocess.run("rm POTCAR", shell=True)
            break
            
def calculate_nelect(path_to_POSCAR, framework_charge=0):

    # Calculate number of electrons
    nelect = 0
    with open(path_to_POSCAR, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                atoms = line.strip().split()
            if index == 5:
                natoms = line.strip().split()
                break
    natoms = [int(i) for i in natoms]
    for i, j in zip(atoms, natoms):
        nelect +=  pbe_valency[i] * j

    nelect = nelect - framework_charge
    print(nelect)
    return nelect
calculate_nelect('POSCAR', framework_charge=-12)

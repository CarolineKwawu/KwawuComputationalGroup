# this script works for when you want to replace the Fe atom that the CO2 is adsorbed on with any of the metals
from ase.io import read, write
from ase import Atoms
import os

# Directory where the original CIF file with the slab adsorbe with CO2 is stored
input_file =r"C:\Users\Amanda\Desktop\My_scripts _for _Fe100\Slab+CO2.cif"

# Directory where the output CIF files will be saved
output_dir =r"C:\Users\Amanda\Desktop\My_scripts _for _Fe100\slabs+CO2+various_metals" 

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# List of metal symbols used for replacement
metals = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 
    'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 
    'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Th', 'Pa', 'U', 'Np', 
    'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
    'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]

# Read the original CIF file
slab = read(input_file)

# Atom index to be replaced
atom_index = 22

# Function to replace atom and save new CIF
def adsorbed_slab(slab, atom_index, metal_symbol, output_dir):
    modified_slab = slab.copy()
    modified_slab[atom_index].symbol = metal_symbol
    output_filename = os.path.join(output_dir, f'Fe_slab_with_{metal_symbol}_+CO2.cif')
    write(output_filename, modified_slab)
    print(f"File created: {output_filename}")

# Loop through each metal and replace the atom
for metal in metals:
    adsorbed_slab(slab, atom_index, metal, output_dir)

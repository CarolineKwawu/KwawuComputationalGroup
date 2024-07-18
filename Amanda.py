#This script creates an Fe 100 slab and dope the slab with each transition metal at position 22 and then adsorbed CO2 to the metal at that position.



from ase.build import bcc100
from ase.io import write, read
from ase.visualize import view
from ase import Atoms
import numpy as np
import os

metals = [
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Rf', 'Db', 'Sg',
    'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Pb', 'Bi', 'Sn', 'Sb'
]

# Function to modify a slab structure by adding a CO2 molecule
def modify_slab_with_CO2(slab):
    # Create a CO2 molecule and adjust it to have a bent shape
    CO2 = Atoms('CO2', positions=[
        [0.0, 0.0, 0.0],      # Carbon atom
        [1.16, 1.0, 0.0],    # Oxygen atom
        [-1.16, 1.0, 0.0]    # Oxygen atom
    ])

    # Set the atom number where the CO2 will be adsorbed
    atom_number = 22

    # Get the position of the specified atom in the slab
    adsorption_site = slab[atom_number].position

    # Set the height above the surface where CO2 will be placed
    adsorption_height = 2.5  # Adjust this value as needed

    # Translate the CO2 molecule to the adsorption site
    CO2.translate(adsorption_site + np.array([0, 0, adsorption_height]))

    # Ensure the CO2 molecule is above the surface by checking the z-coordinate
    if CO2.positions[:, 2].min() < slab.positions[:, 2].max():
        z_shift = slab.positions[:, 2].max() - CO2.positions[:, 2].min() + 1.0
        CO2.translate([0, 0, z_shift])

    # Rotate CO2 molecule to ensure O atoms are bent upwards
    angle = 90.0
    CO2.rotate(angle, v='x', center='COM')

    # Add the CO2 molecule to the slab
    slab += CO2
    return slab

# Create the initial iron slab
fe_slab = bcc100(symbol='Fe', size=(3, 3, 3), a=2.86, vacuum=30, orthogonal=True)

# Define the output directories
initial_output_dir = 'C:/Users/Amanda/Desktop/initial'
adsorbed_output_dir = 'C:/Users/Amanda/Desktop/final'
os.makedirs(initial_output_dir, exist_ok=True)
os.makedirs(adsorbed_output_dir, exist_ok=True)

# Loop through each metal and create doped slabs
for metal in metals:
    modified_slab = fe_slab.copy()
    modified_slab[22].symbol = metal
    initial_filename = f'Fe_slab100_dopedwith{metal}.cif'
    initial_output_filepath = os.path.join(initial_output_dir, initial_filename)
    write(initial_output_filepath, modified_slab)
    
    # Modify the slab by adding a CO2 molecule
    slab_with_CO2 = modify_slab_with_CO2(modified_slab)
    
    # Save the modified slab with CO2 to the adsorbed output directory
    adsorbed_filename = f'Fe_slab100_dopedwith{metal}+CO2.cif'
    adsorbed_output_filepath = os.path.join(adsorbed_output_dir, adsorbed_filename)
    write(adsorbed_output_filepath, slab_with_CO2)
    print(f'Adsorption done for {metal}')
    view(slab_with_CO2)#this line of code is optional

print("Fe slabs have been successfully doped with each transition metal and modified with CO2.")

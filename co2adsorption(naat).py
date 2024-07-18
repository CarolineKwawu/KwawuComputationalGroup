from ase import Atoms
from ase.io import read, write
import numpy as np
import os

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
    adsorption_height = 2.8  # Adjust this value as needed

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

# Define the input and output directories
input_dir = '/Users/mac/Desktop/pristineslabs'
output_dir = '/Users/mac/Desktop/generatedClusters'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Loop through all files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".cif"):
        input_filepath = os.path.join(input_dir, filename)
        
        # Load the slab structure
        slab = read(input_filepath)
        
        # Modify the slab by adding a CO2 molecule
        modified_slab = modify_slab_with_CO2(slab)
        
        # Save the modified slab to the output directory
        output_filepath = os.path.join(output_dir, filename)
        write(output_filepath, modified_slab)
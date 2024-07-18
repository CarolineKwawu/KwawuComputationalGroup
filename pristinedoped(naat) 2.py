import os
from ase.io import write
from ase.build import bcc111

# List of transition metals (excluding La and Ac)
transition_metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                     'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                     'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Rf', 'Db', 'Sg',
                     'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Pb', 'Bi','Sn', 'Sb']

# Define the path to save the output files
output_dir = '/Users/mac/Desktop/outputcifs'
os.makedirs(output_dir, exist_ok=True)

# Loop through each transition metal and create a modified CIF file
for metal in transition_metals:
    # Create the bcc (111) slab of iron
    original_atoms = bcc111(symbol= "Fe", size=(3, 3, 3), a=2.86, vacuum=30, orthogonal=False)

    # Make a copy of the original atoms
    atoms = original_atoms.copy()

    # Replace a specific atom (e.g., the 5th atom) with the transition metal
    atoms[22].symbol = metal
    
    # Write the modified structure to a new CIF file
    output_filename = os.path.join(output_dir, f'Fe_dopedwith_{metal}.cif')
    write(output_filename, atoms)

print("Fe atoms have been successfully been replaced with each transition metal and saved to separate CIF files.")
from ase import Atoms, Atom
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.visualize import view


# Load the structure from a CIF file
input_cif = '/Users/mac/Desktop/generatedClusters/Fe_dopedwith_Cd.cif'
atoms = read(input_cif)

# Identify the surface atoms at the negative z-axis
surface_threshold = 1.5 # Adjust this threshold based on your surface definition
min_z = min(atoms.positions[:, 2])
surface_atoms = [atom.index for atom in atoms if atom.position[2] > min_z + surface_threshold]

# Constrain all atoms except the surface atoms
mask = [atom.index not in surface_atoms for atom in atoms]
constraint = FixAtoms(mask=mask)
atoms.set_constraint(constraint)

# Save the modified structure to a new CIF file
output_cif = 'Fe_output_constrained.cif'
write(output_cif, atoms)

# To view the output structure with ASE GUI and see the constraints
# Uncomment the following lines:
from ase.visualize import view
view(atoms)
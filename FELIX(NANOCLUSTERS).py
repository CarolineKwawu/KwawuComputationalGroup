import ase
import ase.atoms
from ase.cluster import wulff_construction

surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)] # surface of the clusters
esurf = [2.25, 2.29, 2.48] #energy  surface of the Febcc (use your own)
lc = 2.56

# Function to create and store clusters
def create_clusters_and_store(start_size, end_size):
    for size in range(start_size, end_size + 1):
        atoms = wulff_construction('Fe', surfaces, esurf, size, 'bcc', rounding='closest', latticeconstant=lc) ##(it works for fcc and other types of clusteres)

        # Check and iterate until we get the desired number of atoms
        while len(atoms) < size:
            size += 1
            atoms = wulff_construction('Fe', surfaces, esurf, size, 'bcc', rounding='closest', latticeconstant=lc) #(u can chage the closest to above or below)

        atoms = ase.Atoms(atoms)
        atoms.center(vacuum=6.0)
        atoms.set_pbc(True)
        filename = f"Febcc{len(atoms)}.cif" # Use the size as part of the filename
        atoms.write(filename)
        # print(f"Cluster with {len(atoms)} atoms stored in {filename}")
print("Creating of clusters  is in process ......")

# Define the range for cluster sizes
start_size = 10 # Starting cluster size
end_size = 840 # Ending cluster size

# Create and store clusters
create_clusters_and_store(start_size, end_size) 
 

print("All clusters created and stored successfully.")

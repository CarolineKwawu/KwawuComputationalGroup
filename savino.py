from ase import Atoms
from ase.visualize import view
from ase.build import bcc110, add_adsorbate
from ase.io import write, Trajectory
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.calculators.lj import LennardJones
import matplotlib.pyplot as plt


####################### BUILDING SLAB AND CONSTRAINTS
#created the bent CO2 molecule straight away
adsorbate = Atoms('CO2', positions=[[0,0,0],[0,0.76,0.58],[0,-0.76,0.58]])

#slab, added a Pbc
adsorbent = bcc110 ('Fe', size=(3,3,3), a=2.86, vacuum=10.0, orthogonal=False, periodic=True)
#slab.center(axis=2, vacuum=4.0)
#appying constraints to second and third layer (each layer has a tag(0,1,2))
constraints = FixAtoms(mask=[atom.tag > 1 for atom in adsorbent])
adsorbent.set_constraint(constraints)

#view adsorbent structure after setting constraints
view(adsorbent)


########################## DOPING and ADSORPTION

#get a target atom for doping (and adsorption site)
target_atom = adsorbent[7]

#dope target atom
target_atom.symbol = 'Pt'

#view adsorbent after doping
view(adsorbent)

#adsorption site is the cordinates of the target atom
add_adsorbate(adsorbent, adsorbate, 1.5, position=(target_atom.position[0], target_atom.position[1]),)

#write to a cif file
write('doped.cif', adsorbent)

#VIEW after adsorption
view(adsorbent)

##########################################
#DOPING WITH OTHER TRANSITION METALS AND CREATING CIFS
transition_metals = [
   "Sc", "Ti", "V", "Cr", 
    "Mn", "Co", "Ni", 
    "Cu", "Zn", "Y", "Zr", 
    "Nb", "Mo", "Tc", "Ru", 
    "Rh", "Pd", "Ag", "Cd", 
    "Hf", "Ta", "W", "Re", 
    "Os", "Ir", "Pt", "Au", 
    "Hg" 
]

for i in transition_metals:
    new_slab = adsorbent
    for x in new_slab:
        if x.symbol != 'Fe' and x.symbol != 'C' and x.symbol !='O':
            x.symbol = i
            view(new_slab)
    write(f'./structures/doped/Fe_doped_with_Trans/Fe_doped_with_{i}.cif', new_slab)



###### STRUCTURE OPTIMIZATION
#adsorbent.set_calculator(LennardJones())

#optimize_system = BFGS(adsorbent, trajectory = 'opt_traj.traj')
#optimize_system.run(fmax=0.05)

#view(adsorbent)


############################### PLOTTING GRAPH FOR THE OPTIMIZATION
#traj = Trajectory('opt_traj.traj')
#energies = [atoms.get_potential_energy() for atoms in traj]

# Plot energy vs. step number
#plt.plot(energies, marker='o')
#plt.xlabel('Optimization Step')
#plt.ylabel('Potential Energy (eV)')
#plt.title('Energy vs. Optimization Step')
#plt.grid(True)
#plt.show()



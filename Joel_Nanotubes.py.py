# Loading Important libraries and Modules
import os
import ase
from ase.build import nanotube
from ase.io import write, read
from ase.build import molecule
from ase import Atoms


# Set required parameters for the nanotubes
initial_length = 1
final_length = 50
Vacuum = 100

# Generating carbon dioxide molecules for adsorption
co2 = molecule("CO2")
second_co2 = molecule("CO2")
third_co2 = molecule("CO2")
fourth_co2 = molecule("CO2")
fifth_co2 = molecule("CO2")
sixth_co2 = molecule("CO2")
seventh_co2 = molecule("CO2")

# Starting with pristine nanotubes and their adsorption first
# Atoms for Doping
Dopants = ['Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Zr','Mo','Pd','Ag','Cd','Pt']

# Folder for pristine iron nanotubes
pristine_nanotube_directory = "/home/joel/Videos/Nanotubes/Pristine_Fe_nanotubes"
os.makedirs(pristine_nanotube_directory, exist_ok=True)

# Folder for CO2 adsorbed iron nanotubes on edge, surface and within
CO2_adsorbed_pristine_nanotube_directory = "/home/joel/Pictures/Adsorbed_Pristine_Fe_nanotubes"
os.makedirs(CO2_adsorbed_pristine_nanotube_directory, exist_ok=True)

# Folder for CO2 edge adsorbed iron nanotubes
CO2_edge_adsorbed_pristine_nanotube_directory = "/home/joel/Pictures/Adsorbed_Pristine_Fe_nanotubes/Edge_Adsorbed_Fe_nanotubes"
os.makedirs(CO2_edge_adsorbed_pristine_nanotube_directory, exist_ok=True)

# Folder for CO2 surface adsorbed iron nanotubes
CO2_surface_adsorbed_pristine_nanotube_directory = "/home/joel/Pictures/Adsorbed_Pristine_Fe_nanotubes/Surface_Adsorbed_Fe_nanotubes"
os.makedirs(CO2_surface_adsorbed_pristine_nanotube_directory, exist_ok=True)

# Folder for CO2 interior adsorbed iron nanotubes
CO2_interior_adsorbed_pristine_nanotube_directory = "/home/joel/Pictures/Adsorbed_Pristine_Fe_nanotubes/Interior_Adsorbed_Fe_nanotubes"
os.makedirs(CO2_interior_adsorbed_pristine_nanotube_directory, exist_ok=True)




def Pristine_nanotube_generation():
    # Creating pristine Fe nanotubes
    for i in range(initial_length, final_length):
        tube = nanotube(10,0, length = i, bond = 2.40, symbol= 'Fe', vacuum= Vacuum)
        filename = f"Pristine_Fe_nanotube_{i}.cif"
        filepath = os.path.join(pristine_nanotube_directory,filename)
        write(filepath, tube)


    # Adsorbing CO2 on pristine Fe nanotubes

        # Adsorb CO2 at the rim of the nanotube and put them in the new directory
        coordinates = tube[0].position
        co2[0].position = coordinates
        co2[0].position[0] += 1.2

        co2[1].position = co2[0].position 
        co2[1].position[1] += 1
        co2[1].position[0] += 0.8

        co2[2].position = co2[0].position
        co2[2].position[1] += -1
        co2[2].position[0] += 0.8

        tube.extend(co2)

        filename_1_edge_adsorbed = f"Fe_nanotube_{i}_edge_adsorbed.cif"
        filepath_1_edge_adsorbed = os.path.join(CO2_edge_adsorbed_pristine_nanotube_directory,filename_1_edge_adsorbed)
        write(filepath_1_edge_adsorbed, tube)

        # Adsorb CO2 on the surface of the different nanotube and put them in the new directory
        tube2 = nanotube(10,0,length = i, bond =2.40, symbol = 'Fe', vacuum= Vacuum)
        optimal_positions = []
        x_initial = coordinates[0]
        y_initial = coordinates[1]
        z_initial = coordinates[2]
        for j in range(1,len(tube2), 1):
            test_atom = tube[j].position
            x_final = test_atom[0]
            y_final = test_atom[1]
            z_final = test_atom[2]
            if x_initial == x_final and y_initial == y_final:
                best_position = [x_final] + [y_final] + [z_final]
                optimal_positions.append(best_position)
                
        second_adsorption_index = len(optimal_positions)//2
        second_adsorption_coordinates = optimal_positions[second_adsorption_index]
        
        # Previous code involves getting multiple nanotube atoms with same x and y coordinates but varying z coordinates to get an atom in line with nanotube atom number 0

        second_co2[0].position = second_adsorption_coordinates
        second_co2[0].position[0] += 1.5

        second_co2[1].position = second_co2[0].position
        second_co2[1].position[2] += 1
        second_co2[1].position[0] += 0.8

        second_co2[2].position = second_co2[0].position
        second_co2[2].position[2] += -1
        second_co2[2].position[0] += 0.8

        tube2.extend(second_co2)
        
        filename_2_surface_adsorbed = f"Fe_Nanotube_{i}_surface_adsorbed.cif"
        filepath_2_surface_adsorbed = os.path.join(CO2_surface_adsorbed_pristine_nanotube_directory, filename_2_surface_adsorbed)
        write(filepath_2_surface_adsorbed, tube2)


        # Adsorb within nanotube and save them in same directory as other adsorbed nanotubes
        tube3 = nanotube(10,0,length = i, bond =2.40, symbol = 'Fe', vacuum= Vacuum)
        optimal_positions_within = []
        x_initial_within = coordinates[0]
        y_initial_within = coordinates[1]
        z_initial_within = coordinates[2]
        for k in range(1,len(tube3), 1):
            test_atom_within = tube3[k].position
            x_final_within = test_atom_within[0]
            y_final_within = test_atom_within[1]
            z_final_within = test_atom_within[2]
            if x_initial_within == x_final_within and y_initial_within == y_final_within:
                best_position_within = [x_final_within] + [y_final_within] + [z_final_within]
                optimal_positions_within.append(best_position_within)
                
        third_adsorption_index = len(optimal_positions_within)//2
        third_adsorption_coordinates = optimal_positions_within[third_adsorption_index]
        
        # Same procedure as before

        third_co2[0].position = third_adsorption_coordinates
        third_co2[0].position[0] += -1.5

        third_co2[1].position = third_co2[0].position
        third_co2[1].position[1] += 0.8
        third_co2[1].position[0] += -0.8

        third_co2[2].position = third_co2[0].position
        third_co2[2].position[1] += -0.8
        third_co2[2].position[0] += -0.8

        tube3.extend(third_co2)
        
        filename_3_adsorbed_within = f"Fe_Nanotube_{i}_adsorbed_within.cif"
        filepath_3_adsorbed_within = os.path.join(CO2_interior_adsorbed_pristine_nanotube_directory, filename_3_adsorbed_within)
        write(filepath_3_adsorbed_within, tube3)



def Doped_nanotube_data_manipulation():
    # Create directories for all each doped nanotubes
    for a in Dopants:
        # Folder for doped nanotubes
        doped_nanotube_directory = f"/home/joel/Videos/Nanotubes/ALL_DOPED_FE_NANOTUBES/{a}_doped_Fe_nanotube"
        os.makedirs(doped_nanotube_directory, exist_ok=True)

        # Folder for edge doped nanotubes
        edge_doped_nanotube_directory = f"/home/joel/Videos/Nanotubes/ALL_DOPED_FE_NANOTUBES/{a}_doped_Fe_nanotube/{a}_edge_doped_Fe_nanotube"
        os.makedirs(edge_doped_nanotube_directory, exist_ok=True)

        # Folder for surface doped nanotubes
        surface_doped_nanotube_directory = f"/home/joel/Videos/Nanotubes/ALL_DOPED_FE_NANOTUBES/{a}_doped_Fe_nanotube/{a}_surface_doped_Fe_nanotube"
        os.makedirs(surface_doped_nanotube_directory, exist_ok=True)

        # Folder for Adsorbed doped nanotubes
        CO2_adsorbed_doped_nanotube_directory = f"/home/joel/Pictures/Adsorbed_{a}_doped_Fe_nanotubes"
        os.makedirs(CO2_adsorbed_doped_nanotube_directory, exist_ok=True)

        # Folder for CO2 edge adsorbed doped iron nanotubes
        CO2_edge_adsorbed_doped_nanotube_directory = f"/home/joel/Pictures/Adsorbed_{a}_doped_Fe_nanotubes/Edge_Adsorbed_{a}_doped_Fe_nanotubes"
        os.makedirs(CO2_edge_adsorbed_doped_nanotube_directory, exist_ok=True)

        # Folder for CO2 surface adsorbed doped iron nanotubes
        CO2_surface_adsorbed_doped_nanotube_directory = f"/home/joel/Pictures/Adsorbed_{a}_doped_Fe_nanotubes/Surface_Adsorbed_{a}_doped_Fe_nanotubes"
        os.makedirs(CO2_surface_adsorbed_doped_nanotube_directory, exist_ok=True)

        # Folder for CO2 interior adsorbed doped iron nanotubes
        CO2_interior_adsorbed_doped_nanotube_directory = f"/home/joel/Pictures/Adsorbed_{a}_doped_Fe_nanotubes/Interior_Adsorbed_{a}_doped_Fe_nanotubes"
        os.makedirs(CO2_interior_adsorbed_doped_nanotube_directory, exist_ok=True)

        # Folder for CO2 interior adsorbed doped iron nanotubes for both edge 
        CO2_edge_interior_adsorbed_doped_nanotube_directory = f"/home/joel/Pictures/Adsorbed_{a}_doped_Fe_nanotubes/Interior_Adsorbed_{a}_doped_Fe_nanotubes/EDGE"
        os.makedirs(CO2_edge_interior_adsorbed_doped_nanotube_directory, exist_ok=True)

        # Folder for CO2 interior adsorbed doped iron nanotubes for both edge
        CO2_surface_interior_adsorbed_doped_nanotube_directory = f"/home/joel/Pictures/Adsorbed_{a}_doped_Fe_nanotubes/Interior_Adsorbed_{a}_doped_Fe_nanotubes/SURFACE"
        os.makedirs(CO2_surface_interior_adsorbed_doped_nanotube_directory, exist_ok=True)


    # Create doped nanotubes on the edge and dope
        for j in range(initial_length, final_length):
            doped_tube_edge = nanotube(10,0,length = j, bond = 2.40, symbol = 'Fe', vacuum=Vacuum)
            doped_tube_edge[0].symbol = a
            filename_doped_edge = f"{a}_edge_doped_Fe_nanotube_{j}.cif"
            filepath_doped_edge = os.path.join(edge_doped_nanotube_directory,filename_doped_edge)
            write(filepath_doped_edge,doped_tube_edge)

        # Adsorb CO2 on edge doped nanotubes
            coordinates_edge = doped_tube_edge[0].position
            fourth_co2[0].position = coordinates_edge
            fourth_co2[0].position[0] += 1.6

            fourth_co2[1].position = fourth_co2[0].position 
            fourth_co2[1].position[1] += 1
            fourth_co2[1].position[0] += 0.8

            fourth_co2[2].position = fourth_co2[0].position
            fourth_co2[2].position[1] += -1
            fourth_co2[2].position[0] += 0.8

            duplicate_doped_tube_edge = doped_tube_edge.copy()
            duplicate_doped_tube_edge.extend(fourth_co2)

            filename_doped_edge_adsorbed = f"{a}_doped_Fe_nanotube_{j}_edge_adsorbed.cif"
            filepath_doped_edge_adsorbed = os.path.join(CO2_edge_adsorbed_doped_nanotube_directory,filename_doped_edge_adsorbed)
            write(filepath_doped_edge_adsorbed, duplicate_doped_tube_edge)

        # Adsorb CO2 on edge INTERIOR doped naotubes
            sixth_co2[0].position = coordinates_edge
            sixth_co2[0].position[0] += -1.6

            sixth_co2[1].position = sixth_co2[0].position 
            sixth_co2[1].position[1] += 1
            sixth_co2[1].position[0] += -0.8

            sixth_co2[2].position = sixth_co2[0].position
            sixth_co2[2].position[1] += -1
            sixth_co2[2].position[0] += -0.8

            duplicate_doped_tube_edge_for_interior_adsorption = doped_tube_edge.copy()
            duplicate_doped_tube_edge_for_interior_adsorption.extend(sixth_co2)

            filename_doped_edge_interior_adsorbed = f"{a}_doped_Fe_nanotube_{j}_edge_interior_adsorbed.cif"
            filepath_doped_edge_interior_adsorbed = os.path.join(CO2_edge_interior_adsorbed_doped_nanotube_directory,filename_doped_edge_interior_adsorbed)
            write(filepath_doped_edge_interior_adsorbed, duplicate_doped_tube_edge_for_interior_adsorption)


    # Create doped nanotubes on the surface and dope
        for k in range(initial_length, final_length):
            doped_tube_surface = nanotube(10,0,length = k, bond = 2.40, symbol = 'Fe', vacuum=Vacuum)
            optimal_doping_positions = []
            edge_coordinate = doped_tube_surface[0].position
            x_surf_coord = edge_coordinate[0]
            y_surf_coord = edge_coordinate[1]
            z_surf_coord = edge_coordinate[2]
            for l in range(1, len(doped_tube_surface), 1):
                test_atom = doped_tube_surface[l].position
                test_x = test_atom[0]
                test_y = test_atom[1]
                test_z = test_atom[2]
                if x_surf_coord == test_x and y_surf_coord == test_y:
                    optimal_doping_positions.append(l)
            doping_index = len(optimal_doping_positions)//2
            surface_doping_index = optimal_doping_positions[doping_index]
            doped_tube_surface[surface_doping_index].symbol = a

            filename_doped_surface = f"{a}_doped_Fe_nanotube_{k}.cif"
            filepath_doped_surface = os.path.join(surface_doped_nanotube_directory,filename_doped_surface)
            write(filepath_doped_surface,doped_tube_surface)

        # Adsorb CO2 on surface doped nanotubes
            duplicate_doped_tube_surface = doped_tube_surface.copy()
            surface_doping_index_coordinates = duplicate_doped_tube_surface[surface_doping_index].position
            fifth_co2[0].position = surface_doping_index_coordinates
            fifth_co2[0].position[0] += 1.5

            fifth_co2[1].position = fifth_co2[0].position
            fifth_co2[1].position[2] += 1
            fifth_co2[1].position[0] += 0.8

            fifth_co2[2].position = fifth_co2[0].position
            fifth_co2[2].position[2] += -1
            fifth_co2[2].position[0] += 0.8

            duplicate_doped_tube_surface.extend(fifth_co2)
            
            filename_doped_surface_adsorbed = f"{a}_doped_Fe_Nanotube_{k}_surface_adsorbed.cif"
            filepath_doped_surface_adsorbed = os.path.join(CO2_surface_adsorbed_doped_nanotube_directory, filename_doped_surface_adsorbed)
            write(filepath_doped_surface_adsorbed, duplicate_doped_tube_surface)

        # Adsorb CO2 on surface INTERIOR doped nanotubes
            duplicate_doped_tube_surface_interior = doped_tube_surface.copy()
            seventh_co2[0].position = surface_doping_index_coordinates
            seventh_co2[0].position[0] += -1.5

            seventh_co2[1].position = seventh_co2[0].position
            seventh_co2[1].position[2] += 1
            seventh_co2[1].position[0] += -0.8

            seventh_co2[2].position = seventh_co2[0].position
            seventh_co2[2].position[2] += -1
            seventh_co2[2].position[0] += -0.8

            duplicate_doped_tube_surface_interior.extend(seventh_co2)
            
            filename_doped_surface_interior_adsorbed = f"{a}_doped_Fe_Nanotube_{k}_surface_interior_adsorbed.cif"
            filepath_doped_surface_interior_adsorbed = os.path.join(CO2_surface_interior_adsorbed_doped_nanotube_directory, filename_doped_surface_interior_adsorbed)
            write(filepath_doped_surface_interior_adsorbed, duplicate_doped_tube_surface_interior)






Pristine_nanotube_generation()
Doped_nanotube_data_manipulation()






        











            
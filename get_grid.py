import numpy as np
import MDAnalysis as mda
import os

path = "/media/didymos/Projects/ribosome_tunnels/"
os.chdir(path)

bac_id = ['4ybb','abau','bbur','bsub','cacn','drad','ecoli','efae','fjoh','linn','llac','lmon','mpne','msme','mtub','paer','pura','saur','tthe']
euk_id = ['calb','cthe','dmel','drer','ecun','egra','glam','hsap','klac','ldon','mmus','ncra','ntab','ocun','pfal','ploc','scer','slop','slyc','spom','sscr','taes','tbru','tcru','tgon','tthe','tvag']
arc_id = ['hmar','pfur','saci','tkod']



# Parameters
grid_resolution = 1.0  # Grid resolution in Angstroms

global_min_dim = dict()
global_max_dim = dict()
for item in ["10", "20", "30", "40", "60"]:
    all_min_dims = []
    all_max_dims = []
    # BACTERIA
    for name in bac_id:
        print(name)
        min_dims = []
        max_dims = []
        os.chdir(path+"methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_FME_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        all_min_dims.append(np.min(min_dims, axis=0))
        all_max_dims.append(np.max(max_dims, axis=0))
    # EUKARYOTA
    for name in euk_id:
        print(name)
        min_dims = []
        max_dims = []
        os.chdir(path+"eukaryota/methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_MET_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        all_min_dims.append(np.min(min_dims, axis=0))
        all_max_dims.append(np.max(max_dims, axis=0))
    # ARCHAEA
    for name in arc_id:
        print(name)
        min_dims = []
        max_dims = []
        os.chdir(path+"archaea/methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_MET_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for # Parameters
grid_resolution = 1.0  # Grid resolution in Angstroms

global_min_dim = dict()
global_max_dim = dict()
for item in ["10", "20", "30", "40", "60"]:
    all_min_dims = []
    all_max_dims = []
    # BACTERIA
    for name in bac_id:
        print(name)
        min_dims = []
        max_dims = []
        os.chdir(path+"methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_FME_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        all_min_dims.append(np.min(min_dims, axis=0))
        all_max_dims.append(np.max(max_dims, axis=0))
    # EUKARYOTA
    for name in euk_id:
        print(name)
        min_dims = []
        max_dims = []
        os.chdir(path+"eukaryota/methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_MET_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        all_min_dims.append(np.min(min_dims, axis=0))
        all_max_dims.append(np.max(max_dims, axis=0))
    # ARCHAEA
    for name in arc_id:
        print(name)
        min_dims = []
        max_dims = []
        os.chdir(path+"archaea/methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_MET_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        all_min_dims.append(np.min(min_dims, axis=0))
        all_max_dims.append(np.max(max_dims, axis=0))
    global_min_dim[item] = np.min(all_min_dims, axis=0)
    global_max_dim[item] = np.max(all_max_dims, axis=0)
the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Determine the overall bounding box for the trajectory
        for ts in u.trajectory:
            min_dims.append(u.atoms.positions.min(axis=0))
            max_dims.append(u.atoms.positions.max(axis=0))
        all_min_dims.append(np.min(min_dims, axis=0))
        all_max_dims.append(np.max(max_dims, axis=0))
    global_min_dim[item] = np.min(all_min_dims, axis=0)
    global_max_dim[item] = np.max(all_max_dims, axis=0)


file = open(path+"grid_data.dat","w")
for key in global_min_dim.keys():
    min_str = ' '.join(map(str, global_min_dim[key].tolist()))
    max_str = ' '.join(map(str, global_max_dim[key].tolist()))
    # Write the key and the joined string to the file, separated by a space
    file.write(key + ' ' + min_str + '\n')
    file.write(key + ' ' + max_str + '\n')

file.close()


# Define the grid based on overall dimensions
grid_shape = np.ceil((global_max_dim - global_min_dim) / grid_resolution).astype(int)
occupancy_grid = np.zeros(grid_shape)

import MDAnalysis as mda
import numpy as np
import mrcfile
import os

#######################################################################################
# This script takes a trajectory and convert it into 3D map that represents occupancy #
#######################################################################################

def map_traj(u,grid_shape,grid_dim_min,grid_resolution):
    # Function maps atoms to the grid for each frame
    # Calculates occupancy as a fraction of the total number of frames
    occupancy_grid = np.zeros(grid_shape)
    for ts in u.trajectory:
        # Adjust positions based on the overall grid
        positions = (u.atoms.positions - grid_dim_min) / grid_resolution
        # Convert positions to grid indices
        indices = np.round(positions).astype(int)
        # Ensure indices are within grid bounds
        indices = np.clip(indices, 0, grid_shape - 1)
        # Update occupancy grid - using a set to avoid double counting in a single frame
        unique_indices = {tuple(ind) for ind in indices}
        for ind in unique_indices:
            occupancy_grid[ind] += 1
    occupancy_grid /= len(u.trajectory)
    return occupancy_grid

def write_map(name,occupancy_grid,grid_resolution,grid_dim_min):
    with mrcfile.new(name, overwrite=True) as mrc:
        mrc.set_data(occupancy_grid.T.astype(np.float32))
        mrc.voxel_size = grid_resolution
        mrc.header.origin.x = grid_dim_min[item][0]
        mrc.header.origin.y = grid_dim_min[item][1]
        mrc.header.origin.z = grid_dim_min[item][2]

################
path = "/media/didymos/Projects/ribosome_tunnels/"
os.chdir(path)

grid_file = np.loadtxt(path+"grid_data.dat")
grid_dim_min = dict()
grid_dim_max = dict()
ii = 0
for i in grid_file:
    if (~ii%2): grid_dim_min[str(int(i[0]))] = i[1:]
    else: grid_dim_max[str(int(i[0]))] = i[1:]
    ii+=1



bac_id = ['4ybb','abau','bbur','bsub','cacn','drad','ecoli','efae','fjoh','linn','llac','lmon','mpne','msme','mtub','paer','pura','saur','tthe']
euk_id = ['calb','cthe','dmel','drer','ecun','egra','glam','hsap','klac','ldon','mmus','ncra','ntab','ocun','pfal','ploc','scer','slop','slyc','spom','sscr','taes','tbru','tcru','tgon','tthe','tvag']
arc_id = ['hmar','pfur','saci','tkod']

# Parameters
grid_resolution = 1.0  # Grid resolution in Angstroms


for item in ["10", "20", "30", "40", "60"]:
    # Define the grid based on overall dimensions
    grid_shape = np.ceil((grid_dim_max[item] - grid_dim_min[item]) / grid_resolution).astype(int)
    # BACTERIA
    for name in bac_id:
        print(name,item)
        os.chdir(path+"methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_FME_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "ref_1/occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "ref_2/occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)
    # EUKARYOTA
    for name in euk_id:
        print(name,item)
        os.chdir(path+"eukaryota/methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_MET_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "ref_1/occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "ref_2/occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)
    # ARCHAEA
    for name in arc_id:
        print(name,item)
        os.chdir(path+"archaea/methionine/"+name+"/"+item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        topology_file = 'NC_MET_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "ref_1/occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate occupancy grid
        occupancy_grid = map_traj(u,grid_shape,grid_dim_min[item],grid_resolution)
        # Save the occupancy grid to an MRC file
        name = "ref_2/occupancy_map.mrc"
        write_map(name,occupancy_grid,grid_resolution,grid_dim_min)

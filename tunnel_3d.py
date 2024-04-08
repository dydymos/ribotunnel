import MDAnalysis as mda
import numpy as np
import mrcfile

# Parameters
trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
topology_file = 'NC_FME_40.pdb'  # Change to your topology file path
grid_resolution = 1.0  # Grid resolution in Angstroms

# Load the MD trajectory
u = mda.Universe(topology_file, trajectory_file)

# Determine the overall bounding box for the trajectory
all_min_dims = []
all_max_dims = []

for ts in u.trajectory:
    all_min_dims.append(u.atoms.positions.min(axis=0))
    all_max_dims.append(u.atoms.positions.max(axis=0))

global_min_dim = np.min(all_min_dims, axis=0)
global_max_dim = np.max(all_max_dims, axis=0)

# Define the grid based on overall dimensions
grid_shape = np.ceil((global_max_dim - global_min_dim) / grid_resolution).astype(int)
occupancy_grid = np.zeros(grid_shape)

# Map atoms to the grid for each frame
for ts in u.trajectory:
    # Adjust positions based on the overall grid
    positions = (u.atoms.positions - global_min_dim) / grid_resolution
    # Convert positions to grid indices
    indices = np.floor(positions).astype(int)
    # Ensure indices are within grid bounds
    indices = np.clip(indices, 0, grid_shape - 1)
    # Update occupancy grid - using a set to avoid double counting in a single frame
    unique_indices = {tuple(ind) for ind in indices}
    for ind in unique_indices:
        occupancy_grid[ind] += 1

# Calculate occupancy as a fraction of the total number of frames
occupancy_grid /= len(u.trajectory)

# Save the occupancy grid to an MRC file
with mrcfile.new("new_map.mrc", overwrite=True) as mrc:
    mrc.set_data(occupancy_grid.T.astype(np.float32))
    mrc.voxel_size = grid_resolution
    mrc.header.origin.x = global_min_dim[0]
    mrc.header.origin.y = global_min_dim[1]
    mrc.header.origin.z = global_min_dim[2]

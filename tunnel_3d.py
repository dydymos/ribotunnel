import MDAnalysis as mda
import numpy as np

# Parameters
trajectory_file = 'your_trajectory_file.xtc'  # Change to your trajectory file path
topology_file = 'your_topology_file.gro'  # Change to your topology file path
grid_resolution = 1.0  # Grid resolution in Angstroms

# Load the MD trajectory
u = mda.Universe(topology_file, trajectory_file)

# Define the grid based on system dimensions
min_dim = np.floor(u.atoms.positions.min(axis=0))
max_dim = np.ceil(u.atoms.positions.max(axis=0))
grid_shape = np.ceil((max_dim - min_dim) / grid_resolution).astype(int)
occupancy_grid = np.zeros(grid_shape)

# Map atoms to the grid for each frame
for ts in u.trajectory:
    # Adjust positions based on the grid
    positions = (u.atoms.positions - min_dim) / grid_resolution
    # Convert positions to grid indices
    indices = np.floor(positions).astype(int)
    # Ensure indices are within grid bounds
    indices = np.clip(indices, 0, grid_shape - 1)
    # Update occupancy grid
    occupancy_grid[tuple(indices.T)] += 1

# Calculate occupancy as a fraction of the total number of frames
occupancy_grid /= len(u.trajectory)

# occupancy_grid now contains the occupancy values for each grid cell

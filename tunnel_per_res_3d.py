import MDAnalysis as mda
import numpy as np
import mrcfile
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import splprep, splev
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


def calculate_asphericity_from_grid(occupancy_grid):
    # Identify occupied points
    occupied_indices = np.argwhere(occupancy_grid > 0)
    # Calculate centroid of occupied points
    centroid = np.mean(occupied_indices, axis=0)
    # Compute distances from each occupied point to the centroid
    distances = np.sqrt(np.sum((occupied_indices - centroid)**2, axis=1))
    # Calculate mean and standard deviation of the distances
    mean_distance = np.mean(distances)
    std_distance = np.std(distances)
    # Calculate asphericity
    asphericity = std_distance / mean_distance if mean_distance != 0 else 0
    return asphericity


def map_res_traj(u,output_name,grid_resolution):
    # Initialize dictionaries to store the overall min/max positions for each residue
    resid_min = {r: np.inf * np.ones(3) for r in u.residues.resids}
    resid_max = {r: -np.inf * np.ones(3) for r in u.residues.resids}
    path_coords = []
    # Pre-select atoms for each residue to avoid redundant selections
    residue_atoms = {r: u.select_atoms(f"resid {r} and name N CA C O") for r in u.residues.resids}
    # Iterate through the trajectory to update the overall min/max positions
    for ts in u.trajectory:
        for r, atoms in residue_atoms.items():
            positions = atoms.positions
            resid_min[r] = np.minimum(resid_min[r], positions.min(axis=0))
            resid_max[r] = np.maximum(resid_max[r], positions.max(axis=0))
    # Define the grid based on overall dimensions
    grid_shape = {r: np.ceil((resid_max[r] - resid_min[r]) / grid_resolution).astype(int) for r in u.residues.resids}
    occupancy_grid = {r: np.zeros(grid_shape[r]) for r in u.residues.resids}
    # Generate occupancy grid for each residue
    for ts in u.trajectory:
        for r, atoms in residue_atoms.items():
            # Adjust positions based on the overall grid
            positions = (atoms.positions - resid_min[r]) / grid_resolution
            # Convert positions to grid indices
            indices = np.round(positions).astype(int)
            # Ensure indices are within grid bounds
            indices = np.clip(indices, 0, grid_shape[r] - 1)
            # Update occupancy grid - using a set to avoid double counting in a single frame
            unique_indices = {tuple(ind) for ind in indices}
            for ind in unique_indices:
                occupancy_grid[r][ind] += 1
    # Find the maximum occupancy for each residue and put an atom there
    # Plus calculate asphericity and mean distance
    file = open(output_name,'w')
    asphericity = []
    path_coords = []
    for r in occupancy_grid:
        occupancy_grid[r] /= len(u.trajectory)
        # Find the max occupancy value and its grid indices for each residue
        max_index = np.unravel_index(np.argmax(occupancy_grid[r]), occupancy_grid[r].shape)
        # Convert grid indices to real-space coordinates
        max_coords = np.array(max_index) * grid_resolution + resid_min[r]
        path_coords.append(max_coords)
        file.write(pdb_format.format(r, atom_type, max_coords[0], max_coords[1], max_coords[2], atom_type))
        # Asphericity
        asphericity.append(calculate_asphericity_from_grid(occupancy_grid[r]))
    file.close()
    # Caclulating RMSF
    ca_atoms = u.select_atoms("name CA")
    R=rms.RMSF(ca_atoms).run()
    return R.results.rmsf,asphericity,path_coords


def write_map(name,occupancy_grid,grid_resolution,grid_dim_min):
    with mrcfile.new(name, overwrite=True) as mrc:
        mrc.set_data(occupancy_grid.T.astype(np.float32))
        mrc.voxel_size = grid_resolution
        mrc.header.origin.x = grid_dim_min[item][0]
        mrc.header.origin.y = grid_dim_min[item][1]
        mrc.header.origin.z = grid_dim_min[item][2]


def get_tunnel_path(path_coords,item,s):
    x, y, z = np.array(path_coords).T
    # Fit the spline to the data points
    tck, u = splprep([x, y, z], s=s)
    # Generate new interpolated points from the spline representation
    tunnel_path = splev(np.linspace(0, 1, int(int(item)*3)), tck)
    return tunnel_path


def plot_tunnel_path(path_coords,tunnel_path):
    fig = plt.figure()
    x = np.array(path_coords).T[0]
    y = np.array(path_coords).T[1]
    z = np.array(path_coords).T[2]
    x_spread = x.max() - x.min()
    y_spread = y.max() - y.min()
    z_spread = z.max() - z.min()
    max_spread = max(x_spread, y_spread, z_spread)
    # Calculate the center of each dimension
    x_center = 0.5 * (x.max() + x.min())
    y_center = 0.5 * (y.max() + y.min())
    z_center = 0.5 * (z.max() + z.min())
    # Calculate the new limits for each axis, centered around the respective means
    x_limits = [x_center - max_spread / 2, x_center + max_spread / 2]
    y_limits = [y_center - max_spread / 2, y_center + max_spread / 2]
    z_limits = [z_center - max_spread / 2, z_center + max_spread / 2]
    # Plotting
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, facecolors='none', edgecolors='k', label='Average position of the nascent chain residue')
    ax.plot(tunnel_path[0], tunnel_path[1], tunnel_path[2], lw=2, label='Fitted Curve at s=20')
    # Set the calculated limits for each axis
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_zlim(z_limits)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()
    plt.show()


def occupancy_edge(occupancy_map,epsilon):
    # list of occupied cells
    occ_indices = np.array(np.where(occupancy_map>epsilon)).T
    vectors = []
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            for z in [-1, 0, 1]:
                vectors.append((x, y, z))
    vectors.remove((0, 0, 0))
    edges = []
    for indx in occ_indices:
        # Calculate absolute indices for all neighbors
        neighbor_indices = indx + vectors
        # Checking how many neighbours have occupancy > epsilon
        occupied_neighbors = np.sum(occupancy_map[tuple(neighbor_indices.T)] > epsilon)
        if occupied_neighbors < 26:
            edges.append(indx)
    return np.array(edges)

################
################
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

bac_name = ['E.coli','A.baumannii','B.burgdorferi','B.subtilis','C.acnes','D.radiodurans','E.coli','E.faecalis','F.johnsoniae','L.innocua','L.lactis','L.monocytogenes','M.pneumoniae','M.smegmatis','M.tuberculosis','P.aeruginosa','P.urativorans','S.aureus','T.thermophilus']

bac_dict = dict()
i=0
for name in bac_id:
    bac_dict[name] = dict()
    bac_dict[name]['name'] = bac_name[i]
    bac_dict[name]['rmsf'] = dict()
    bac_dict[name]['asph'] = dict()
    bac_dict[name]['tunnel'] = dict()
    bac_dict[name]['tunnel path'] = dict()
    bac_dict[name]['tunnel interpol'] = dict()
    bac_dict[name]['tunnel ind'] = dict()
    i+=1


# PDB format template
pdb_format = "ATOM  {:5d}  {:<4s}MOL     1    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           {}\n"
atom_type = "C"

# Parameters
grid_resolution = 1.0  # Grid resolution in Angstroms

for item in ["10", "20", "30", "40", "60"]:
    grid_shape = np.ceil((grid_dim_max[item] - grid_dim_min[item]) / grid_resolution).astype(int)
    # BACTERIA
    for name in bac_id:
        print(name,item)
        trajectory_file = 'fitted.xtc'  # Change to your trajectory file path
        bac_dict[name]['rmsf'][item] = dict()
        bac_dict[name]['asph'][item] = dict()
        bac_dict[name]['tunnel'][item] = dict()
        os.chdir(path+"methionine/"+name+"/"+item)
        topology_file = 'NC_FME_'+item+'.pdb'  # Change to your topology file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate max occupancy path and residue based properties
        output_name = path+"methionine/"+name+"/"+item+"/"+"occ_path.pdb"
        rmsf,asph,tunnel = map_res_traj(u,output_name,grid_resolution)
        bac_dict[name]['rmsf'][item]['0'] = rmsf
        bac_dict[name]['asph'][item]['0'] = asph
        bac_dict[name]['tunnel'][item]['0'] = tunnel
        # REF_1
        trajectory_file = 'ref_1/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate max occupancy path and residue based properties
        output_name = path+"methionine/"+name+"/"+item+"/ref_1/"+"occ_path.pdb"
        rmsf,asph,tunnel = map_res_traj(u,output_name,grid_resolution)
        bac_dict[name]['rmsf'][item]['1'] = rmsf
        bac_dict[name]['asph'][item]['1'] = asph
        bac_dict[name]['tunnel'][item]['1'] = tunnel
        # REF_2
        trajectory_file = 'ref_2/fitted.xtc'  # Change to your trajectory file path
        # Load the MD trajectory
        u = mda.Universe(topology_file, trajectory_file)
        # Calculate max occupancy path and residue based properties
        output_name = path+"methionine/"+name+"/"+item+"/ref_2/"+"occ_path.pdb"
        rmsf,asph,tunnel = map_res_traj(u,output_name,grid_resolution)
        bac_dict[name]['rmsf'][item]['2'] = rmsf
        bac_dict[name]['asph'][item]['2'] = asph
        bac_dict[name]['tunnel'][item]['2'] = tunnel




# Calculate spline interpolation of the tunnel path
s = 100
for name in bac_id:
    for item in ["10", "20", "30", "40", "60"]:
        bac_dict[name]['tunnel path'][item] = dict()
        bac_dict[name]['tunnel path'][item]['0'] = get_tunnel_path(bac_dict[name]['tunnel'][item]['0'],item,s)
        bac_dict[name]['tunnel path'][item]['1'] = get_tunnel_path(bac_dict[name]['tunnel'][item]['1'],item,s)
        bac_dict[name]['tunnel path'][item]['2'] = get_tunnel_path(bac_dict[name]['tunnel'][item]['2'],item,s)



for name in bac_id:
    for item in ["10", "20", "30", "40", "60"]:
        bac_dict[name]['tunnel interpol'][item] = dict()
        bac_dict[name]['tunnel ind'][item] = dict()
        bac_dict[name]['tunnel ind'][item]['0'] = [np.argmin(np.sum((np.array(bac_dict[name]['tunnel path'][item]['0']).T-point)**2,1)) for point in bac_dict[name]['tunnel'][item]['0']]
        bac_dict[name]['tunnel ind'][item]['1'] = [np.argmin(np.sum((np.array(bac_dict[name]['tunnel path'][item]['1']).T-point)**2,1)) for point in bac_dict[name]['tunnel'][item]['1']]
        bac_dict[name]['tunnel ind'][item]['2'] = [np.argmin(np.sum((np.array(bac_dict[name]['tunnel path'][item]['2']).T-point)**2,1)) for point in bac_dict[name]['tunnel'][item]['2']]
        bac_dict[name]['tunnel interpol'][item]['0'] = np.array(bac_dict[name]['tunnel path'][item]['0']).T[np.array(bac_dict[name]['tunnel ind'][item]['0'])]
        bac_dict[name]['tunnel interpol'][item]['1'] = np.array(bac_dict[name]['tunnel path'][item]['1']).T[np.array(bac_dict[name]['tunnel ind'][item]['1'])]
        bac_dict[name]['tunnel interpol'][item]['2'] = np.array(bac_dict[name]['tunnel path'][item]['2']).T[np.array(bac_dict[name]['tunnel ind'][item]['2'])]

n = bac_dict[name]['tunnel ind'][item]['0'][1]
vektor = np.array(bac_dict[name]['tunnel path'][item]['0']).T[n+1] - np.array(bac_dict[name]['tunnel path'][item]['0']).T[n-1]
vec_from_point = edges_coord - np.array(bac_dict[name]['tunnel path'][item]['0']).T[n]
# Calculate the dot product
dot_products = np.dot(vec_from_point, vektor)
# Define a small threshold
threshold = 2  # Adjust this value based on your specific requirements

# Find indices where the absolute value of the dot product is less than the threshold
perpendicular_indices = np.where(np.abs(dot_products) < threshold)[0]

# Select the edge points that are close to being perpendicular to 'vektor'
perpendicular_edge_points = edges_coord[perpendicular_indices]

atom_type = "P"
output_name = path+"methionine/"+name+"/"+item+"/"+"occ_edge_perp.pdb"
file = open(output_name,"w")
for r in range(0,len(perpendicular_edge_points)):
    file.write(pdb_format.format(r+1, atom_type, perpendicular_edge_points[r][2], perpendicular_edge_points[r][1], perpendicular_edge_points[r][0], atom_type))


file.close()

# RMSD between tunnel paths from max points
for name in bac_id:
    bac_dict[name]['rmsd_avg'] = dict()
    for item in ["10", "20", "30", "40", "60"]:
        rmsd = []
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel'][item]['0'])-np.array(bac_dict[name]['tunnel'][item]['1']))**2,1)))
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel'][item]['0'])-np.array(bac_dict[name]['tunnel'][item]['2']))**2,1)))
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel'][item]['1'])-np.array(bac_dict[name]['tunnel'][item]['2']))**2,1)))
        bac_dict[name]['rmsd_avg'][item] = np.mean(np.array(rmsd),0)


# RMSD between tunnel paths interpolations
for name in bac_id:
    bac_dict[name]['rmsd_avg_int'] = dict()
    for item in ["10", "20", "30", "40", "60"]:
        rmsd = []
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel path'][item]['0'])-np.array(bac_dict[name]['tunnel path'][item]['1']))**2,0)))
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel path'][item]['0'])-np.array(bac_dict[name]['tunnel path'][item]['2']))**2,0)))
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel path'][item]['1'])-np.array(bac_dict[name]['tunnel path'][item]['2']))**2,0)))
        bac_dict[name]['rmsd_avg_int'][item] = np.mean(np.array(rmsd),0)

# RMSD between tunnel paths via points interpolation
for name in bac_id:
    bac_dict[name]['rmsd_avg_int_point'] = dict()
    for item in ["10", "20", "30", "40", "60"]:
        rmsd = []
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel interpol'][item]['0'])-np.array(bac_dict[name]['tunnel interpol'][item]['1']))**2,1)))
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel interpol'][item]['0'])-np.array(bac_dict[name]['tunnel interpol'][item]['2']))**2,1)))
        rmsd.append(np.sqrt(np.sum((np.array(bac_dict[name]['tunnel interpol'][item]['1'])-np.array(bac_dict[name]['tunnel interpol'][item]['2']))**2,1)))
        bac_dict[name]['rmsd_avg_int_point'][item] = np.mean(np.array(rmsd),0)



########### Occupancy MAP ##########
for name in bac_id:
    for item in ["10", "20", "30", "40", "60"]:



occupancy_file = mrcfile.open(path+"methionine/"+name+"/"+item+"/occupancy_map.mrc")
occupancy_map = occupancy_file.data
occupancy_origin = np.array(occupancy_file.header['origin'].tolist())
occupancy_voxel = occupancy_file.voxel_size.tolist()[0]


# Calculating edge of the occupancy map
edges = occupancy_edge(occupancy_map,0.001)
edges_coord = np.array(edges)*occupancy_voxel+occupancy_origin[::-1]








atom_type = "N"
output_name = path+"methionine/"+name+"/"+item+"/"+"occ_edge.pdb"
file = open(output_name,"w")
edges_coord = np.array(edges)*occupancy_voxel+occupancy_origin[::-1]
for r in range(0,len(edges)):
    file.write(pdb_format.format(r+1, atom_type, edges_coord[r][2], edges_coord[r][1], edges_coord[r][0], atom_type))


file.close()



########### PLOTTING #########

# RMSF plot for each residue
plt.plot(bac_dict[name]['rmsf'][item]['0'][:-1])
plt.plot(bac_dict[name]['rmsf'][item]['1'][:-1])
plt.plot(bac_dict[name]['rmsf'][item]['2'][:-1])
plt.show()

# Asphericity plot for each residue
plt.plot(bac_dict[name]['asph'][item]['0'][:-1])
plt.plot(bac_dict[name]['asph'][item]['1'][:-1])
plt.plot(bac_dict[name]['asph'][item]['2'][:-1])
plt.show()


plt.plot(bac_dict[name]['rmsd_avg_int'][item])
plt.plot(bac_dict[name]['rmsd_avg'][item])
plt.plot(bac_dict[name]['rmsd_avg_int_point'][item])


#### 3D PLOTs
fig = plt.figure()
path_coords = bac_dict[name]['tunnel'][item]['2']
x = np.array(path_coords).T[0]
y = np.array(path_coords).T[1]
z = np.array(path_coords).T[2]
x_spread = x.max() - x.min()
y_spread = y.max() - y.min()
z_spread = z.max() - z.min()
max_spread = max(x_spread, y_spread, z_spread)
# Calculate the center of each dimension
x_center = 0.5 * (x.max() + x.min())
y_center = 0.5 * (y.max() + y.min())
z_center = 0.5 * (z.max() + z.min())
# Calculate the new limits for each axis, centered around the respective means
x_limits = [x_center - max_spread / 2, x_center + max_spread / 2]
y_limits = [y_center - max_spread / 2, y_center + max_spread / 2]
z_limits = [z_center - max_spread / 2, z_center + max_spread / 2]
# Plotting
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.array(bac_dict[name]['tunnel'][item]['0'])[:,0], np.array(bac_dict[name]['tunnel'][item]['0'])[:,1], np.array(bac_dict[name]['tunnel'][item]['0'])[:,2])
ax.scatter(np.array(bac_dict[name]['tunnel'][item]['1'])[:,0], np.array(bac_dict[name]['tunnel'][item]['1'])[:,1], np.array(bac_dict[name]['tunnel'][item]['1'])[:,2])
ax.scatter(np.array(bac_dict[name]['tunnel'][item]['2'])[:,0], np.array(bac_dict[name]['tunnel'][item]['2'])[:,1], np.array(bac_dict[name]['tunnel'][item]['2'])[:,2])
ax.plot(bac_dict[name]['tunnel path'][item]['0'][0], bac_dict[name]['tunnel path'][item]['0'][1], bac_dict[name]['tunnel path'][item]['0'][2], lw=2, label='Fitted Curve at s='+str(s))
ax.plot(bac_dict[name]['tunnel path'][item]['1'][0], bac_dict[name]['tunnel path'][item]['1'][1], bac_dict[name]['tunnel path'][item]['1'][2], lw=2, label='Fitted Curve at s='+str(s))
ax.plot(bac_dict[name]['tunnel path'][item]['2'][0], bac_dict[name]['tunnel path'][item]['2'][1], bac_dict[name]['tunnel path'][item]['2'][2], lw=2, label='Fitted Curve at s='+str(s))
# Set the calculated limits for each axis
ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
ax.set_zlim(z_limits)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.legend()
plt.show()






















###########################################
# Calculations of the whole occupancy map #
###########################################
for item in ["10", "20", "30", "40", "60"]:
    grid_shape = np.ceil((grid_dim_max[item] - grid_dim_min[item]) / grid_resolution).astype(int)
    # BACTERIA
    for name in bac_id:
        print(name)
        os.chdir(path+"methionine/"+name+"/"+item)
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
        print(name)
        os.chdir(path+"eukaryota/methionine/"+name+"/"+item)
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
        print(name)
        os.chdir(path+"archaea/methionine/"+name+"/"+item)
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

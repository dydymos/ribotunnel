import os
import glob
import numpy as np
import itertools
import mrcfile
import tqdm
from tqdm import *
import MDAnalysis as mda
from sklearn.cluster import DBSCAN
from scipy.interpolate import splprep, splev
from MDAnalysis.lib.distances import capped_distance

# PDB format templates
# pdb_format = "ATOM  {:5d}  {:<4s}MOL     1    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           {}\n"
pdb_format_beta = "ATOM  {:5d}  {:<4s}MOL     1    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}            {}\n"
pdb_format = "ATOM  {:5d}  {:<4s}MOL  {:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           {}\n"

# Standard amino acid chemical classifications
polar = ['SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN']
hydrophobic = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP']
positive = ['LYS', 'ARG', 'HIS']
negative = ['ASP', 'GLU']
special = ['GLY', 'PRO']

##########################################
# Generates dictionary to store all data #
##########################################

def get_dictionary(organisms_dict):
    """
    Generates a nested dictionary to store all properties for given organisms.
    Expects a dictionary format: {'org_id': 'org_name'}
    """
    # Grouping arguments across multiple lines makes it easier to read/edit
    arguments = [
        "occupancy map avg", "res occ path", "res occ path avg", 
        "spline path", "spline path avg", "edges coord", 
        "edges coord avg", "interior coord avg", "path length", 
        "dist mean", "dist asph", "dist mean avg", "dist asph avg", 
        "proteins", "rna", "proteins_matrix", "rna_matrix"
    ]
    org_dict = {}
    for org_id, org_name in organisms_dict.items():
        org_dict[org_id] = {'name': org_name}
        for arg in arguments:
            org_dict[org_id][arg] = {}       
    return org_dict


###########################################################
# Reading grids from previously calculated occupancy maps #
###########################################################
def read_grid(grid_file_path):
    """
    Reads grid dimensions from a data file. 
    Assumes alternating rows correspond to minimum and maximum grid dimensions.
    """
    grid_data = np.loadtxt(grid_file_path)
    grid_dim_min = {}
    grid_dim_max = {}
    # enumerate() automatically provides an index (0, 1, 2...) for each row
    for index, row in enumerate(grid_data):
        key = str(int(row[0])) # The identifier
        coords = row[1:]       # The grid coordinates
        # Even indices (0, 2, 4...) are min; Odd indices (1, 3, 5...) are max
        if index % 2 == 0: 
            grid_dim_min[key] = coords
        else: 
            grid_dim_max[key] = coords
    return grid_dim_min, grid_dim_max


###############################################
# Reading occupancy pathways for each residue #
###############################################

def read_occupancy_path(base_path, org_dict, org_ids, short=False):
    """
    Reads PDB files containing occupancy paths for different nascent chain lengths.
    Averages coordinates across the main trajectory and two replicates.
    """
    if not short:
        chain_lengths = ["10", "20", "30", "40", "60"]
    else: 
        chain_lengths = ["10", "20", "30", "40"]
    for name in org_ids:
        for length in chain_lengths:
            org_dict[name]['res occ path'][length] = {}
            # os.path.join safely constructs paths regardless of operating system
            path_0 = os.path.join(base_path, name, length, "occ_max_res.pdb")
            path_1 = os.path.join(base_path, name, length, "ref_1", "occ_max_res.pdb")
            path_2 = os.path.join(base_path, name, length, "ref_2", "occ_max_res.pdb")
            # Load Universes
            u_0 = mda.Universe(path_0)
            u_1 = mda.Universe(path_1)
            u_2 = mda.Universe(path_2)
            # Extract positions
            pos_0 = u_0.atoms.positions
            pos_1 = u_1.atoms.positions
            pos_2 = u_2.atoms.positions
            # Store in dictionary
            org_dict[name]['res occ path'][length]['0'] = pos_0
            org_dict[name]['res occ path'][length]['1'] = pos_1
            org_dict[name]['res occ path'][length]['2'] = pos_2
            # Avg position of the path (axis=0 averages across the 3 replicates)
            org_dict[name]['res occ path avg'][length] = np.mean([pos_0, pos_1, pos_2], axis=0)



#################################################
# Calculates the spline function for every path #
#################################################

def get_spline_path(path_coords, points, s):
    """
    Fits a 3D spline to coordinate data and returns interpolated points.
    Removes consecutive duplicate coordinates before fitting.
    """
    # removing duplicates from the list if any
    new_list = [path_coords[0]]  
    for i in range(1, len(path_coords)):
        if not np.array_equal(path_coords[i], path_coords[i-1]):  
            new_list.append(path_coords[i])    
    path_coords = np.array(new_list) 
    x, y, z = path_coords.T
    # Fit the spline to the data points
    tck, u = splprep([x, y, z], s=s)
    # Generate new interpolated points from the spline representation
    tunnel_path = splev(np.linspace(0, 1, int(points)), tck)
    return tunnel_path

def calculate_spline(org_dict, org_ids, s, NN, short=False):
    """
    Calculates splines for individual trajectories and averaged paths.
    Handles coordinate truncation when the nascent chain exits the tunnel.
    """
    if not short:
        chain_lengths = ["10", "20", "30", "40", "60"]
    else: 
        chain_lengths = ["10", "20", "30", "40"]
    for name in org_ids:
        for length in chain_lengths:
            org_dict[name]['spline path'][length] = {}
            org_dict[name]['spline path avg'][length] = {}
            # Create local references for cleaner code
            res_path = org_dict[name]['res occ path'][length]
            avg_path = org_dict[name]['res occ path avg'][length]
            points = int(length) * 3
            # Calculate splines for the 3 individual trajectories
            org_dict[name]['spline path'][length]['0'] = get_spline_path(res_path['0'], points, s)
            org_dict[name]['spline path'][length]['1'] = get_spline_path(res_path['1'], points, s)
            org_dict[name]['spline path'][length]['2'] = get_spline_path(res_path['2'], points, s)
            if not short and length == "60":
                # NN is a cutoff as with length = 60 we are outside the ribosome tunnel
                mean_outside = np.mean(avg_path[:NN], axis=0)
                # np.vstack cleanly stacks the 1D mean_outside array on top of the 2D remaining path
                extended = np.vstack((mean_outside, avg_path[NN:]))
                org_dict[name]['spline path avg'][length] = get_spline_path(extended, (int(length) - NN + 1) * 3, s)
            elif short and length == "40":
                org_dict[name]['spline path avg'][length] = get_spline_path(avg_path[NN:], (int(length) - NN) * 3, s)
            else: 
                org_dict[name]['spline path avg'][length] = get_spline_path(avg_path, points, s)

#################
# Writing paths #
#################

def write_paths(base_path, org_dict, org_ids, NN, short=False):
    """
    Writes the averaged occupancy paths and calculated splines 
    into standard PDB format for visualization.
    """
    if not short:
        chain_lengths = ["10", "20", "30", "40", "60"]
    else: 
        chain_lengths = ["10", "20", "30", "40"]
    atom_type = "H"
    for name in org_ids:
        for length in chain_lengths:
            # Safely join paths
            output_name_occ = os.path.join(base_path, name, length, "occ_max_res_avg.pdb")
            output_name_spline = os.path.join(base_path, name, length, "spline_path_avg.pdb")
            # Determine which pathway data to use based on chain length cutoff
            if not short and length == "60":
                pathway = org_dict[name]['res occ path avg'][length][NN:]
            else:
                pathway = org_dict[name]['res occ path avg'][length]
            spline = org_dict[name]['spline path avg'][length]
            # Use 'with' to automatically handle file opening and closing
            with open(output_name_occ, "w") as file_occ:
                for r in range(len(pathway[:,0])):
                    # Using format: Atom ID, Atom Name, Res ID, X, Y, Z, Element Symbol
                    file_occ.write(pdb_format.format(r+1, atom_type, 1, pathway[r][0], pathway[r][1], pathway[r][2], atom_type))
            with open(output_name_spline, "w") as file_spline:
                for r in range(len(spline[0])):
                    file_spline.write(pdb_format.format(r+1, atom_type, 1, spline[0][r], spline[1][r], spline[2][r], atom_type))


#################################################################
# Calculating where are edges and interior in the occupancy map #
#################################################################

def get_occupancy_edge_and_interior(occupancy_map, epsilon):
    """
    Identifies edge and interior voxels of the occupancy map.
    A voxel is considered 'interior' if it has 23 or more occupied neighbors.
    """
    # np.argwhere perfectly replaces np.array(np.where(...)).T
    occ_indices = np.argwhere(occupancy_map >= epsilon)
    # Generate the 26 neighbor vectors using itertools (cleaner than 3 nested loops)
    vectors = [vec for vec in itertools.product([-1, 0, 1], repeat=3) if vec != (0, 0, 0)]
    vectors = np.array(vectors)
    edges = []
    interior = []
    # Get grid boundaries to prevent wrap-around indexing
    z_max, y_max, x_max = occupancy_map.shape
    for indx in occ_indices:
        neighbor_indices = indx + vectors
        # SAFETY CHECK: Filter out neighbors that fall outside the grid bounds
        valid_z = (neighbor_indices[:, 0] >= 0) & (neighbor_indices[:, 0] < z_max)
        valid_y = (neighbor_indices[:, 1] >= 0) & (neighbor_indices[:, 1] < y_max)
        valid_x = (neighbor_indices[:, 2] >= 0) & (neighbor_indices[:, 2] < x_max)
        valid_neighbors = neighbor_indices[valid_z & valid_y & valid_x]
        # Check how many valid neighbors have occupancy >= epsilon
        occupied_neighbors = np.sum(occupancy_map[tuple(valid_neighbors.T)] >= epsilon)
        if occupied_neighbors < 23:
            edges.append(indx)
        else: 
            interior.append(indx)
    return np.array(edges), np.array(interior)

#################################################################################
# Getting edges and interior of the occupancy maps at certain map threshold     #
# Generating average occupancy map and corresponding edges.                     #
#################################################################################

def get_edge_from_occupancy(base_path, grid_dim_min, grid_resolution, org_dict, organisms_id, threshold, short=False):
    if not short:
        chain_lengths = ["10", "20", "30", "40", "60"]
    else: 
        chain_lengths = ["10", "20", "30", "40"]    
    atom_type = "N"
    # tqdm will give a nice progress bar in your notebook!
    for name in tqdm(organisms_id, desc="Processing Occupancy Maps"):
        for length in chain_lengths:
            org_dict[name]['edges coord'][length] = {}
            sum_occupancy_map = None
            # Use a dictionary to loop through the 3 replicates effortlessly
            replicates = {"0": "", "1": "ref_1", "2": "ref_2"}
            for rep_key, rep_folder in replicates.items():
                mrc_path = os.path.join(base_path, name, length, rep_folder, "occupancy_map.mrc")
                # 'with' ensures the MRC file is properly closed after reading
                with mrcfile.open(mrc_path) as occupancy_file:
                    occupancy_map = occupancy_file.data
                    occupancy_origin = np.array(occupancy_file.header['origin'].tolist())
                    occupancy_voxel = occupancy_file.voxel_size.tolist()[0]
                # Accumulate for the average map
                if sum_occupancy_map is None:
                    sum_occupancy_map = np.copy(occupancy_map)
                else:
                    sum_occupancy_map += occupancy_map
                # Calculate edges
                edges, _ = get_occupancy_edge_and_interior(occupancy_map, threshold)
                # Convert to physical coords. Note: [::-1] flips ZYX to XYZ for PDB formatting
                edges_coord = np.array(edges) * occupancy_voxel + occupancy_origin[::-1]
                # We can cleanly extract columns rather than typing edges_coord[:,2] repeatedly
                org_dict[name]['edges coord'][length][rep_key] = edges_coord[:, [2, 1, 0]]
                # Write edges PDB
                output_name = os.path.join(base_path, name, length, rep_folder, "occ_edge.pdb")
                with open(output_name, "w") as file:
                    for r, coord in enumerate(edges_coord):
                        file.write(pdb_format.format(r+1, atom_type, 1, coord[2], coord[1], coord[0], atom_type))
            # --- Process the Averaged Map ---
            avg_occupancy_map = sum_occupancy_map / 3.0
            org_dict[name]['occupancy map avg'] = avg_occupancy_map
            # Calculate edge and interior of the averaged map
            edges, interior = get_occupancy_edge_and_interior(avg_occupancy_map, threshold)
            edges_coord = np.array(edges) * occupancy_voxel + occupancy_origin[::-1]
            interior_coord = np.array(interior) * occupancy_voxel + occupancy_origin[::-1]
            org_dict[name]['edges coord avg'][length] = edges_coord[:, [2, 1, 0]]
            org_dict[name]['interior coord avg'][length] = interior_coord[:, [2, 1, 0]]
            # Write averaged edges PDB
            output_name = os.path.join(base_path, name, length, "occ_edge_avg.pdb")
            with open(output_name, "w") as file:
                for r, coord in enumerate(edges_coord):
                    file.write(pdb_format.format(r+1, atom_type, 1, coord[2], coord[1], coord[0], atom_type))
            # Write averaged interior PDB
            output_name = os.path.join(base_path, name, length, "occ_inter_avg.pdb")
            with open(output_name, "w") as file:
                for r, coord in enumerate(interior_coord):
                    file.write(pdb_format.format(r+1, "H", 1, coord[2], coord[1], coord[0], "H"))
            # Write averaged MRC map
            output_name = os.path.join(base_path, name, length, "occupancy_map_avg.mrc")
            write_map(output_name, avg_occupancy_map, grid_resolution, grid_dim_min[length])

################
# Writing maps #
################
def write_map(name, occupancy_grid, grid_resolution, grid_dim_min):
    """Saves a NumPy array as an MRC file with proper origin mapping."""
    with mrcfile.new(name, overwrite=True) as mrc:
        mrc.set_data(occupancy_grid.astype(np.float32))
        mrc.voxel_size = grid_resolution
        mrc.header.origin.x = grid_dim_min[0]
        mrc.header.origin.y = grid_dim_min[1]
        mrc.header.origin.z = grid_dim_min[2]

import os
import glob
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm


############################################################################
# Retrieving part of the ribosome structure that is close to the all edges #
# Along with main properties of this tunnel                                #
############################################################################

# Generate selection strings dynamically for readability
aa_names = hydrophobic + polar + positive + negative + special
aa_names_q = [res + '?' for res in aa_names] # Handle custom naming like ALA?

prot_sel_str = "resname " + " ".join(aa_names + aa_names_q)
# RNA selection: Everything that is NOT a protein residue
rna_sel_str = f"not ({prot_sel_str})"

def get_tunnel(base_path, org_dict, org_ids, dist_threshold=5.0, short=False):
    if not short:
        chain_lengths = ["10", "20", "30", "40", "60"]
    else: 
        chain_lengths = ["10", "20", "30", "40"]
    pbar = tqdm(org_ids, desc="Analyzing Tunnels")
    for name in pbar:
        pbar.set_description(f"Analyzing {name}")
        rib_path = os.path.join(base_path, name)
        rib_file = glob.glob(os.path.join(rib_path, '*_tunnel_fix.pdb'))[0]
        ribosome = mda.Universe(rib_file)
        all_atoms = ribosome.atoms
        all_positions = all_atoms.positions
        org_dict[name]["proteins"] = {}
        org_dict[name]["rna"] = {}
        for length in chain_lengths:
            edge_coords = org_dict[name]['edges coord avg'][length]
            # MEMORY-SAFE FAST CALCULATION:
            # capped_distance uses a KD-tree. It only calculates distances up to 'max_cutoff'.
            # Returns 'pairs' (indices of points that are close) and 'distances' (actual distances)
            pairs, _ = capped_distance(edge_coords, all_positions, max_cutoff=dist_threshold)
            # pairs[:, 1] contains the indices of the atoms from 'all_positions'
            # np.unique ensures we don't count an atom twice if it's near two different edge points
            close_atom_indices = np.unique(pairs[:, 1])
            # Extract unique atoms and their parent residues
            unique_atoms = all_atoms[close_atom_indices]
            unique_residues = unique_atoms.residues 
            # Get all atoms belonging to these residues to keep residues fully intact
            full_residue_atoms = unique_residues.atoms
            # Classify protein vs RNA
            protein_atoms = full_residue_atoms.select_atoms(prot_sel_str)
            rna_atoms = full_residue_atoms.select_atoms(rna_sel_str)
            # Extract the 3-letter codes for the protein residues
            protein_resnames = np.array([res.resname[:3] for res in protein_atoms.residues])
            # Count physico-chemical properties
            polar_n = np.sum(np.isin(protein_resnames, polar))
            hydrophobic_n = np.sum(np.isin(protein_resnames, hydrophobic))
            positive_n = np.sum(np.isin(protein_resnames, positive))
            negative_n = np.sum(np.isin(protein_resnames, negative))
            special_n = np.sum(np.isin(protein_resnames, special))            
            org_dict[name]["proteins"][length] = [polar_n, hydrophobic_n, positive_n, negative_n, special_n]
            org_dict[name]["rna"][length] = len(rna_atoms.residues)
            # Write the isolated tunnel residues to a new PDB file
            output_pdb = os.path.join(rib_path, length, 'tunnel.pdb')
            with mda.Writer(output_pdb, multiframe=False) as PDB:
                PDB.write(full_residue_atoms)

####################################################
# Function to write total properties to text files #
####################################################

def write_tunnel_prop(org_dict, org_ids, theme, output_path, short=False):
    """Calculates means and standard deviations of properties and writes to .dat files."""
    chain_lengths = ["10", "20", "30", "40"] if short else ["10", "20", "30", "40", "60"]
    means = []
    stds = []
    # Collect the data
    for length in chain_lengths:
        data_matrix = np.array([org_dict[name]["proteins"][length] for name in org_ids])
        means.append(np.mean(data_matrix, axis=0))
        stds.append(np.std(data_matrix, axis=0))
    # Construct safe file paths
    output_name_mean = os.path.join(output_path, f"{theme}_mean_composition.dat")
    output_name_std = os.path.join(output_path, f"{theme}_std_composition.dat")
    # np.savetxt is the fastest, safest way to write 2D arrays to text files.
    # fmt='%.4f' keeps the decimal lengths clean and uniform.
    np.savetxt(output_name_mean, means, fmt='%.4f')
    np.savetxt(output_name_std, stds, fmt='%.4f')


##################################################
# Calculating tunnel rings                       #
# Which are rings build from the edge points     #
# That are perpendicular to the main tunnel path #
##################################################

def get_rings(base_path, org_dict, org_ids, dot_threshold, ref_tunnel_path, short=False):
    """
    Calculates cross-sectional 'rings' of the tunnel by finding edge points
    that are strictly perpendicular to the local tangent of the reference spline.
    Saves these rings as PDB files for visualization.
    """
    atom_type = "C"
    chain_lengths = ["10", "20", "30", "40"] if short else ["10", "20", "30", "40", "60"]  
    # The reference path is likely shape (3, N). We transpose it to (N, 3) for coordinate math.
    ref_path_coords = ref_tunnel_path.T 
    # 1. Get tangent vectors (Difference between consecutive points on the spline)
    tangents = np.diff(ref_path_coords, axis=0)
    # 2. Normalize the tangent vectors
    tangents_n = tangents / np.linalg.norm(tangents, axis=1)[:, np.newaxis]
    # 3. Calculate distance along the spline path
    n_points = len(tangents)
    tangents_len = np.linalg.norm(tangents, axis=1)
    path_l = np.zeros(n_points + 1)
    # Reversing the cumsum tracks distance from the *end* of the tunnel!
    path_l[:-1] = np.cumsum(tangents_len[::-1])[::-1]
    pbar = tqdm(org_ids, desc="Calculating Rings")
    for name in pbar:
        pbar.set_description(f"Calculating Rings: {name}")
        for length in chain_lengths:
            edges_avg = org_dict[name]['edges coord avg'][length]
            for i in range(n_points):
                spline_point = ref_path_coords[i]
                # Vector from the spline to ALL edge points
                vec_to_edge = edges_avg - spline_point
                # Calculate the dot product
                dot_products = np.dot(vec_to_edge, tangents_n[i])
                # Find indices where the absolute value of the dot product < threshold
                perpendicular_indices = np.where(np.abs(dot_products) < dot_threshold)[0]
                if perpendicular_indices.size > 0:
                    perpendicular_edge_points = edges_avg[perpendicular_indices]
                    # Format distance nicely (e.g., "45.20")
                    dist_str = f"{path_l[i]:.2f}"
                    output_name = os.path.join(base_path, name, length, f"ring_avg_{dist_str}.pdb")
                    # Write ring to PDB
                    with open(output_name, "w") as file:
                        for r, point in enumerate(perpendicular_edge_points):
                            file.write(pdb_format.format(
                                r+1, atom_type, r+1, 
                                point[0], point[1], point[2], 
                                atom_type
                            ))


##################################################
# Calculating tunnel discs                       #
# Which are disc build from occupancy map voxels #
# That are perpendicular to the main tunnel path #
##################################################

def cluster_points(perpendicular_interior_points):
    """
    Groups contiguous interior points into clusters. 
    Uses DBSCAN for high-speed spatial clustering based on a 2.0 Angstrom threshold.
    """
    # eps=2.0 matches your original sqrt(3) voxel diagonal tolerance.
    # min_samples=1 ensures every point gets a cluster (no noise points).
    clustering = DBSCAN(eps=2.0, min_samples=1).fit(perpendicular_interior_points)
    # DBSCAN labels start at 0, your original code started at 1. We add 1 to match.
    return clustering.labels_ + 1

def get_interior_edge(perpendicular_interior_points, edge_points):
    """Finds edge points that are immediately adjacent to the main interior cluster."""
    # Fast vectorized distance calculation between two point clouds
    from scipy.spatial.distance import cdist
    distances = cdist(edge_points, perpendicular_interior_points) 
    # Find edge points where the minimum distance to ANY interior point is < 2.0
    close_edges = edge_points[np.min(distances, axis=1) < 2.0]
    return close_edges

def compute_cross_section_asphericity(disc):
    """Calculates asphericity based on the 2D gyration tensor."""
    centroid = disc.mean(axis=0)
    centered = disc - centroid
    # PCA to find best-fit plane
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # Project to 2D plane defined by the two largest eigenvectors
    basis_2D = eigenvectors[:, 1:]
    projected_2D = centered @ basis_2D 
    # Compute 2D gyration tensor
    gyration_tensor = np.cov(projected_2D.T)
    lambda1, lambda2 = np.linalg.eigvalsh(gyration_tensor)
    # Avoid division by zero if the disc is perfectly linear/singular
    if (lambda1 + lambda2) == 0:
        return 0.0
    A = ((lambda1 - lambda2)**2) / (lambda1 + lambda2)**2
    return A

def get_discs(base_path, org_dict, org_ids, dot_threshold, ref_tunnel_path, short=False):
    chain_lengths = ["10", "20", "30", "40"] if short else ["10", "20", "30", "40", "60"]
    atom_type = "P"
    ref_path_coords = ref_tunnel_path.T 
    tangents = np.diff(ref_path_coords, axis=0)
    tangents_n = tangents / np.linalg.norm(tangents, axis=1)[:, np.newaxis]
    n_points = len(tangents)
    tangents_len = np.linalg.norm(tangents, axis=1)
    path_l = np.zeros(n_points + 1)
    path_l[:-1] = np.cumsum(tangents_len[::-1])[::-1]
    pbar = tqdm(org_ids, desc="Processing Discs")
    for name in pbar:
        pbar.set_description(f"Discs: {name}")
        # Initialize dictionary arrays
        org_dict[name].update({
            'path length': {}, 'dist mean avg': {}, 'dist asph avg': {}, 
            'cross section': {}, 'proteins_matrix': {}, 'rna_matrix': {}
        })
        rib_path = os.path.join(base_path, name)
        rib_file = glob.glob(os.path.join(rib_path, '*_tunnel_fix.pdb'))[0]
        ribosome = mda.Universe(rib_file)
        all_atoms = ribosome.atoms
        all_positions = all_atoms.positions
        for length in chain_lengths:
            # Pre-allocate lists and arrays
            org_dict[name]['path length'][length] = path_l
            org_dict[name]['dist mean avg'][length] = []
            org_dict[name]['dist asph avg'][length] = []
            org_dict[name]['cross section'][length] = []
            org_dict[name]["proteins_matrix"][length] = np.zeros([n_points, 5])
            org_dict[name]["rna_matrix"][length] = np.zeros([n_points])
            for i in range(n_points):
                spline_point = ref_path_coords[i]
                vec_to_interior = org_dict[name]['interior coord avg'][length] - spline_point
                dot_products = np.dot(vec_to_interior, tangents_n[i])
                perpendicular_indices = np.where(np.abs(dot_products) < dot_threshold)[0]
                if perpendicular_indices.size > 0:
                    perp_interior = org_dict[name]['interior coord avg'][length][perpendicular_indices]
                    # 1. Cluster Identification
                    cluster_id = cluster_points(perp_interior)
                    # Distance from spline to all interior points
                    dists_to_spline = np.linalg.norm(perp_interior - spline_point, axis=1)
                    if np.min(dists_to_spline) <= 8.0:
                        # Find the cluster that contains the point closest to the spline
                        main_cluster_id = cluster_id[np.argmin(dists_to_spline)]
                        main_interior = perp_interior[cluster_id == main_cluster_id]
                        # Load pre-calculated edge points
                        ring_file = os.path.join(base_path, name, length, f"ring_avg_{path_l[i]:.2f}.pdb")
                        edge_points = mda.Universe(ring_file).atoms.positions
                        # 2. Get active edges
                        interior_edge = get_interior_edge(main_interior, edge_points)
                        # 3. Assemble whole disc and calculate properties
                        whole_disc = np.vstack((main_interior, interior_edge))
                        org_dict[name]['cross section'][length].append(len(main_interior) + len(interior_edge) / 3.0)
                        org_dict[name]['dist mean avg'][length].append(np.mean(np.linalg.norm(interior_edge - spline_point, axis=1)))
                        org_dict[name]['dist asph avg'][length].append(compute_cross_section_asphericity(whole_disc))
                        # Write disc PDB safely
                        output_name = os.path.join(base_path, name, length, f"disc_avg_{path_l[i]:.2f}.pdb")
                        with open(output_name, "w") as file:
                            for r, point in enumerate(whole_disc):
                                file.write(pdb_format.format(r+1, "C", 1, point[0], point[1], point[2], "C"))
                        # 4. FAST Physico-chemical property extraction
                        pairs, _ = capped_distance(interior_edge, all_positions, max_cutoff=5.0)
                        close_atom_indices = np.unique(pairs[:, 1])
                        if len(close_atom_indices) > 0:
                            unique_atoms = all_atoms[close_atom_indices]
                            full_residue_atoms = unique_atoms.residues.atoms
                            # Safely save the tunnel slice
                            slice_out_dir = os.path.join(rib_path, length, "tunnel")
                            os.makedirs(slice_out_dir, exist_ok=True)
                            # Only write if it's not totally empty to avoid MDA errors
                            if len(full_residue_atoms) > 0:
                                full_residue_atoms.write(os.path.join(slice_out_dir, f"{path_l[i]:.2f}.pdb"))
                            # Count using the safe selectors we built earlier
                            protein_atoms = full_residue_atoms.select_atoms(prot_sel_str)
                            rna_atoms = full_residue_atoms.select_atoms(rna_sel_str)
                            protein_resnames = np.array([res.resname[:3] for res in protein_atoms.residues])
                            org_dict[name]["proteins_matrix"][length][i,] = [
                                np.sum(np.isin(protein_resnames, polar)),
                                np.sum(np.isin(protein_resnames, hydrophobic)),
                                np.sum(np.isin(protein_resnames, positive)),
                                np.sum(np.isin(protein_resnames, negative)),
                                np.sum(np.isin(protein_resnames, special))
                            ]
                            org_dict[name]["rna_matrix"][length][i] = len(rna_atoms.residues)
                    else:
                        _append_nans(org_dict[name], length)
                else:
                    _append_nans(org_dict[name], length)

def _append_nans(org_data, length):
    """Helper to keep logic clean when no disc is found."""
    org_data['dist mean avg'][length].append(np.nan)
    org_data['dist asph avg'][length].append(np.nan)
    org_data['cross section'][length].append(np.nan)


def write_tunnel_path(org_dict, org_ids, theme, short=False):
    """Calculates and saves the mean and std of cross-sectional areas."""
    chain_lengths = ["10", "20", "30", "40"] if short else ["10", "20", "30", "40", "60"]
    mean_results = []
    std_results = []
    for item in chain_lengths:
        # Safely extract data into a 2D float array
        path_matrix = np.array([org_dict[org]['cross section'][item] for org in org_ids], dtype=float)
        # Calculate stats ignoring NaNs
        mean_results.append(np.nanmean(path_matrix, axis=0))
        std_results.append(np.nanstd(path_matrix, axis=0))
    # Save the entire 2D array to a space-separated text file instantly
    np.savetxt(f"{theme}_avg_path.dat", mean_results, fmt='%.6f')
    np.savetxt(f"{theme}_std_path.dat", std_results, fmt='%.6f')


def write_tunnel_asph(org_dict, org_ids, theme, short=False):
    """Calculates and saves the mean and std of tunnel asphericity."""
    chain_lengths = ["10", "20", "30", "40"] if short else ["10", "20", "30", "40", "60"]
    mean_results = []
    std_results = []
    for item in chain_lengths:
        # Safely extract data into a 2D float array
        path_matrix = np.array([org_dict[org]['dist asph avg'][item] for org in org_ids], dtype=float)
        # Calculate stats ignoring NaNs
        mean_results.append(np.nanmean(path_matrix, axis=0))
        std_results.append(np.nanstd(path_matrix, axis=0))
    # Save the entire 2D array to a space-separated text file instantly
    np.savetxt(f"{theme}_avg_asph.dat", mean_results, fmt='%.6f')
    np.savetxt(f"{theme}_std_asph.dat", std_results, fmt='%.6f')
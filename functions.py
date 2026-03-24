import os
import glob
import mrcfile
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
from scipy.interpolate import splprep, splev, interp1d
from tqdm import tqdm

# PDB format template
#pdb_format = "ATOM  {:5d}  {:<4s}MOL     1    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           {}\n"
pdb_format_beta = "ATOM  {:5d}  {:<4s}MOL     1    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}           {}\n"
pdb_format = "ATOM  {:5d}  {:<4s}MOL  {:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           {}\n"

# Amino acid properties
hydrophobic = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP']
polar = ['GLN', 'ASN', 'SER', 'THR', 'TYR']
positive = ['LYS', 'ARG', 'HIS']
negative = ['ASP', 'GLU']
special = ['GLY','PRO','CYS']

##########################################
# Generates dictionary to story all data #
##########################################
def get_dictionary(org_id,org_name):
    org_dict = dict()
    arguments = ["occupancy map avg","res occ path","res occ path avg","spline path","spline path avg","edges coord","edges coord avg","interior coord avg","path length","dist mean","dist asph","dist mean avg","dist asph avg","proteins","rna","proteins_matrix","rna_matrix"]
    i=0
    for name in org_id:
        org_dict[name] = dict()
        org_dict[name]['name'] = org_name[i]
        for j in arguments:
            org_dict[name][j] = dict()
        i+=1
    return org_dict

###########################################################
# Reading grids from previously calculated occupancy maps #
###########################################################
def read_grid(grid_file_path):
    grid_file = np.loadtxt(grid_file_path)
    grid_dim_min = dict()
    grid_dim_max = dict()
    ii = 0
    for i in grid_file:
        if (~ii%2): grid_dim_min[str(int(i[0]))] = i[1:]
        else: grid_dim_max[str(int(i[0]))] = i[1:]
        ii+=1
    return grid_dim_min,grid_dim_max


###############################################
# Reading occupancy pathways for each residue #
###############################################
def read_occupancy_path(path,org_dict,org_id,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]
    for name in org_id:
        for item in list:
            org_dict[name]['res occ path'][item] = dict()
            input_name = path+name+"/"+item+"/occ_max_res.pdb"
            input_name_1 = path+name+"/"+item+"/ref_1/occ_max_res.pdb"
            input_name_2 = path+name+"/"+item+"/ref_2/occ_max_res.pdb"
            u = mda.Universe(input_name)
            u_1 = mda.Universe(input_name_1)
            u_2 = mda.Universe(input_name_2)
            org_dict[name]['res occ path'][item]['0'] = u.atoms.positions
            org_dict[name]['res occ path'][item]['1'] = u_1.atoms.positions
            org_dict[name]['res occ path'][item]['2'] = u_2.atoms.positions
            # Avg position of the path
            org_dict[name]['res occ path avg'][item] = np.mean([np.array(org_dict[name]['res occ path'][item]['0']),np.array(org_dict[name]['res occ path'][item]['1']),np.array(org_dict[name]['res occ path'][item]['2'])],0)

#################################################
# Calculates the spline function for every path #
#################################################
def get_spline_path(path_coords,points,s):
    # removing duplicates from the list if any
    new_list = [path_coords[0]]  
    for i in range(1, np.shape(path_coords)[0]):
        if not np.array_equal(path_coords[i], path_coords[i-1]):  
            new_list.append(path_coords[i])
    path_coords = np.array(new_list) 
    x, y, z = np.array(path_coords).T
    # Fit the spline to the data points
    tck, u = splprep([x, y, z], s=s)
    # Generate new interpolated points from the spline representation
    tunnel_path = splev(np.linspace(0, 1, int(points)), tck)
    return tunnel_path

def calculate_spline(org_dict,org_id,s,NN,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]
    for name in org_id:
        for item in list:
            org_dict[name]['spline path'][item] = dict()
            org_dict[name]['spline path avg'][item] = dict()
            org_dict[name]['spline path'][item]['0'] = get_spline_path(org_dict[name]['res occ path'][item]['0'],int(item)*3,s)
            org_dict[name]['spline path'][item]['1'] = get_spline_path(org_dict[name]['res occ path'][item]['1'],int(item)*3,s)
            org_dict[name]['spline path'][item]['2'] = get_spline_path(org_dict[name]['res occ path'][item]['2'],int(item)*3,s)
            if (short==False and item=="60"):
                # NN is a cutoff as with length = 60 we are outside the ribosome tunnel
                # org_dict[name]['spline path avg'][item] = get_spline_path(org_dict[name]['occ path avg'][item][NN:],(int(item)-NN)*3,s)
                mean_outside = np.mean(org_dict[name]['res occ path avg'][item][:NN],axis=0)
                extended = np.zeros([60-NN+1,3])
                extended[1:,0] = org_dict[name]['res occ path avg'][item][NN:,0]
                extended[1:,1] = org_dict[name]['res occ path avg'][item][NN:,1]
                extended[1:,2] = org_dict[name]['res occ path avg'][item][NN:,2]
                extended[0][0] = mean_outside[0]
                extended[0][1] = mean_outside[1]
                extended[0][2] = mean_outside[2]
                org_dict[name]['spline path avg'][item] = get_spline_path(extended,(int(item)-NN+1)*3,s)
            elif (short==True and item=="40"):
                org_dict[name]['spline path avg'][item] = get_spline_path(org_dict[name]['res occ path avg'][item][NN:],(int(item)-NN)*3,s)
            else: 
                org_dict[name]['spline path avg'][item] = get_spline_path(org_dict[name]['res occ path avg'][item],int(item)*3,s)


#################
# Writing paths #
#################
def write_paths(path,org_dict,org_id,NN,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]
    atom_type = "H"
    for name in org_id:
        for item in list:
            output_name_occ = path+name+"/"+item+"/"+"occ_max_res_avg.pdb"
            output_name_spline = path+name+"/"+item+"/"+"spline_path_avg.pdb"
            file_occ = open(output_name_occ,"w")
            file_spline = open(output_name_spline,"w")
            if (short==False and item=="60"):
                pathway = org_dict[name]['res occ path avg'][item][NN:]
                spline = org_dict[name]['spline path avg'][item]
            else:
                pathway = org_dict[name]['res occ path avg'][item]
                spline = org_dict[name]['spline path avg'][item]
            for r in range(0,len(pathway[:,0])):
                file_occ.write(pdb_format.format(r+1, atom_type, 1, pathway[r][0], pathway[r][1], pathway[r][2], atom_type))
            for r in range(0,len(spline[0])):
                file_spline.write(pdb_format.format(r+1, atom_type, 1, spline[0][r], spline[1][r], spline[2][r], atom_type))
            file_spline.close()
            file_occ.close()

#################################################################
# Calculating where are edges and interior in the occupancy map #
# Saving the aerage occupancy map                               #
#################################################################
def get_occupancy_edge_and_interior(occupancy_map,epsilon):
    # list of occupied cells
    occ_indices = np.array(np.where(occupancy_map>=epsilon)).T
    vectors = []
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            for z in [-1, 0, 1]:
                vectors.append((x, y, z))
    vectors.remove((0, 0, 0))
    edges = []
    interior = []
    for indx in occ_indices:
        # Calculate absolute indices for all neighbors
        neighbor_indices = indx + vectors
        # Checking how many neighbours have occupancy > epsilon
        occupied_neighbors = np.sum(occupancy_map[tuple(neighbor_indices.T)] >= epsilon)
        if occupied_neighbors < 23:
            edges.append(indx)
        else: interior.append(indx)
    return np.array(edges),np.array(interior)

def get_edge_from_occupancy(path,grid_dim_min,grid_resolution,org_dict,organisms_id,threshold,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]    
    atom_type = "N"
    for name in tqdm(organisms_id):
        for item in list:
            org_dict[name]['edges coord'][item] = dict()
            occupancy_file = mrcfile.open(path+name+"/"+item+"/occupancy_map.mrc")
            occupancy_map = occupancy_file.data
            occupancy_origin = np.array(occupancy_file.header['origin'].tolist())
            occupancy_voxel = occupancy_file.voxel_size.tolist()[0]
            org_dict[name]['occupancy map avg'] = np.copy(occupancy_map)
            # Calculating edge of the occupancy map
            edges,_ = get_occupancy_edge_and_interior(occupancy_map,threshold)
            edges_coord = np.array(edges)*occupancy_voxel+occupancy_origin[::-1]
            org_dict[name]['edges coord'][item]['0'] = np.array([edges_coord[:,2],edges_coord[:,1],edges_coord[:,0]]).T
            # writing edges
            output_name = path+name+"/"+item+"/"+"occ_edge.pdb"
            file = open(output_name,"w")
            for r in range(0,len(edges)):
                file.write(pdb_format.format(r+1, atom_type, 1, edges_coord[r][2], edges_coord[r][1], edges_coord[r][0], atom_type))
            file.close()
            # ref_1
            occupancy_file = mrcfile.open(path+name+"/"+item+"/ref_1/occupancy_map.mrc")
            occupancy_map = occupancy_file.data
            org_dict[name]['occupancy map avg'] += occupancy_map
            # Calculating edge of the occupancy map
            edges,_ = get_occupancy_edge_and_interior(occupancy_map,threshold)
            edges_coord = np.array(edges)*occupancy_voxel+occupancy_origin[::-1]
            org_dict[name]['edges coord'][item]['1'] = np.array([edges_coord[:,2],edges_coord[:,1],edges_coord[:,0]]).T
            # writing edges
            output_name = path+name+"/"+item+"/"+"ref_1/occ_edge.pdb"
            file = open(output_name,"w")
            for r in range(0,len(edges)):
                file.write(pdb_format.format(r+1, atom_type, 1, edges_coord[r][2], edges_coord[r][1], edges_coord[r][0], atom_type))
            file.close()
            # ref_2
            occupancy_file = mrcfile.open(path+name+"/"+item+"/ref_2/occupancy_map.mrc")
            occupancy_map = occupancy_file.data
            org_dict[name]['occupancy map avg'] += occupancy_map
            org_dict[name]['occupancy map avg'] /= 3
            # Calculating edge of the occupancy map
            edges,_ = get_occupancy_edge_and_interior(occupancy_map,threshold)
            edges_coord = np.array(edges)*occupancy_voxel+occupancy_origin[::-1]
            org_dict[name]['edges coord'][item]['2'] = np.array([edges_coord[:,2],edges_coord[:,1],edges_coord[:,0]]).T
            # writing edges
            output_name = path+name+"/"+item+"/"+"ref_2/occ_edge.pdb"
            file = open(output_name,"w")
            for r in range(0,len(edges)):
                file.write(pdb_format.format(r+1, atom_type, 1, edges_coord[r][2], edges_coord[r][1], edges_coord[r][0], atom_type))
            file.close()
            # Calculating edge and interior of the averaged occupancy map
            edges,interior = get_occupancy_edge_and_interior(org_dict[name]['occupancy map avg'],threshold)
            edges_coord = np.array(edges)*occupancy_voxel+occupancy_origin[::-1]
            interior_coord = np.array(interior)*occupancy_voxel+occupancy_origin[::-1]
            org_dict[name]['edges coord avg'][item] = np.array([edges_coord[:,2],edges_coord[:,1],edges_coord[:,0]]).T
            org_dict[name]['interior coord avg'][item] = np.array([interior_coord[:,2],interior_coord[:,1],interior_coord[:,0]]).T
            # writing averaged edges
            output_name = path+name+"/"+item+"/"+"/occ_edge_avg.pdb"
            file = open(output_name,"w")
            for r in range(0,len(edges)):
                file.write(pdb_format.format(r+1, atom_type, 1, edges_coord[r][2], edges_coord[r][1], edges_coord[r][0], atom_type))
            file.close()
            # writing averaged interior
            output_name = path+name+"/"+item+"/"+"/occ_inter_avg.pdb"
            file = open(output_name,"w")
            for r in range(0,len(interior)):
                file.write(pdb_format.format(r+1, "H", 1,interior_coord[r][2], interior_coord[r][1], interior_coord[r][0], atom_type))
            file.close()
            # writing averaged map
            output_name = path+name+"/"+item+"/"+"/occupancy_map_avg.mrc"
            write_map(output_name,org_dict[name]['occupancy map avg'],grid_resolution,grid_dim_min[item])

################
# Writing maps #
################
def write_map(name,occupancy_grid,grid_resolution,grid_dim_min):
    with mrcfile.new(name, overwrite=True) as mrc:
        mrc.set_data(occupancy_grid.astype(np.float32))
        mrc.voxel_size = grid_resolution
        mrc.header.origin.x = grid_dim_min[0]
        mrc.header.origin.y = grid_dim_min[1]
        mrc.header.origin.z = grid_dim_min[2]

############################################################################
# Retrieving part of the ribosome structure that is close to the all edges #
# Along with main properties of this tunnel                                #
############################################################################
def get_tunnel(path,org_dict,org_id,dist_threshold,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]
    for name in tqdm(org_id):
        print(name)
        rib_path = path+name+"/"
        rib_file = glob.glob(os.path.join(rib_path,'*_tunnel_fix.pdb'))[0]
        ribosome = mda.Universe(rib_file)
        # Select all atoms in the ribosome
        all_atoms = ribosome.select_atoms("all")
        # Threshold distance
        threshold = dist_threshold
        org_dict[name]["proteins"] = dict()
        org_dict[name]["rna"] = dict()
        for item in list:
            # Find atoms within 5A of the given coordinates
            close_atoms = []
            for coord in org_dict[name]['edges coord avg'][item]:
                distances = distance_array(coord, all_atoms.positions)
                within_threshold = np.where(distances <= threshold)[1]
                for atom_index in within_threshold:
                    close_atoms.append(all_atoms[atom_index])
            uniq_atoms = np.unique(close_atoms)
            # Get unique residues containing these atoms
            unique_residues_id = np.unique([atom.resindex for atom in uniq_atoms])
            # Select the unique residues
            unique_residues = ribosome.residues[unique_residues_id]
            # Get all atoms in the unique residues
            unique_atoms = unique_residues.atoms
            # Select protein and RNA atoms from the unique atoms
            protein_atoms = unique_atoms.select_atoms('resname ALA GLY PRO CYS TRP TYR ILE LEU LYS VAL GLN GLU ASP ASN MET THR PHE HIS SER ARG ALA? GLY? PRO? CYS? TRP? TYR? ILE? LEU? LYS? VAL? GLN? GLU? MET? ASN? ASP? THR? PHE? HIS? SER? ARG?')
            rna_atoms = unique_atoms.select_atoms('not resname ALA GLY PRO CYS TRP TYR ILE LEU LYS VAL GLN GLU ASP ASN MET THR PHE HIS SER ARG ALA? GLY? PRO? CYS? TRP? TYR? ILE? LEU? LYS? VAL? GLN? GLU? MET? ASN? ASP? THR? PHE? HIS? SER? ARG?')
            protein_resnames = np.array([x[:3] for x in protein_atoms.residues.resnames])
            #if protein_atoms.n_atoms == 0:
            #    protein_atoms = unique_atoms.select_atoms('resname ALA? GLY? PRO? CYS? TRP? TYR? ILE? LEU? LYS? VAL? GLN? GLU? ASP? ASN? ASP? THR? PHE? HIS? SER? ARG?')
            #    rna_atoms = unique_atoms.select_atoms('not resname ALA? GLY? PRO? CYS? TRP? TYR? ILE? LEU? LYS? VAL? GLN? GLU? ASP? ASN? ASP? THR? PHE? HIS? SER? ARG?')
            #    protein_resnames = np.array([x[:3] for x in protein_atoms.residues.resnames])
            polar_n  = np.sum(np.isin(protein_resnames,polar))
            hydrophobic_n = np.sum(np.isin(protein_resnames,hydrophobic))
            positive_n = np.sum(np.isin(protein_resnames,positive))
            negative_n = np.sum(np.isin(protein_resnames,negative))
            special_n = np.sum(np.isin(protein_resnames,special))
            org_dict[name]["proteins"][item] = [polar_n,hydrophobic_n,positive_n,negative_n,special_n]
            org_dict[name]["rna"][item] = len(rna_atoms.residues)
            # Write the unique residues to a new PDB file
            output_pdb = os.path.join(rib_path+item+"/", 'tunnel.pdb')
            with mda.Writer(output_pdb, multiframe=False) as PDB:
               PDB.write(unique_residues.atoms)

def write_tunnel_prop(org_dict,org_id,theme,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]
    output_name0 = theme+"_mean_composition.dat"
    output_name1 = theme+"_std_composition.dat"
    proteins_prop = dict()
    for item in list:
        proteins_prop[item] = np.zeros([len(org_id),5])
    i=0
    for name in org_id:
        for item in list:
            proteins_prop[item][i] = org_dict[name]["proteins"][item]
        i+=1    
    with open(output_name0, 'w') as f:
        f.write(' '.join(map(str, np.mean(proteins_prop["10"],0)))+"\n")
    for item in list[1:]:
        with open(output_name0, 'a') as f:
            f.write(' '.join(map(str, np.mean(proteins_prop[item],0)))+"\n")
    with open(output_name1, 'w') as f:
        f.write(' '.join(map(str, np.std(proteins_prop["10"],0)))+"\n")
    for item in list[1:]:
        with open(output_name1, 'a') as f:
            f.write(' '.join(map(str, np.std(proteins_prop[item],0)))+"\n")

def write_tunnel_path(org_dict,org_id,theme,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
        tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
    else: 
        list = ["10", "20", "30", "40"]
        tunnel_length = org_dict[org_id[0]]['path length']['40'][:-1]
    path_matrix = np.zeros([len(org_id),tunnel_length.size,])
    with open(theme+"_avg_path.dat", 'w') as f:
            for i,org in enumerate(org_id):
                path_matrix[i] = org_dict[org]['cross section']['10']
            f.write(' '.join(map(str, np.nanmean(path_matrix,0)))+"\n")
    for item in list[1:]:
            for i,org in enumerate(org_id):
                path_matrix[i] = org_dict[org]['cross section'][item]
            with open(theme+"_avg_path.dat", 'a') as f:
                f.write(' '.join(map(str, np.nanmean(path_matrix,0)))+"\n")
    with open(theme+"_std_path.dat", 'w') as f:
            for i,org in enumerate(org_id):
                path_matrix[i] = org_dict[org]['cross section']['10']
            f.write(' '.join(map(str, np.nanstd(path_matrix,0)))+"\n")
    for item in list[1:]:
            for i,org in enumerate(org_id):
                path_matrix[i] = org_dict[org]['cross section'][item]
            with open(theme+"_std_path.dat", 'a') as f:
                f.write(' '.join(map(str, np.nanstd(path_matrix,0)))+"\n")

def write_tunnel_asph(org_dict,org_id,theme,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]
    tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
    path_matrix = np.zeros([len(org_id),tunnel_length.size,])
    with open(theme+"_avg_asph.dat", 'w') as f:
            for i,org in enumerate(org_id):
                path_matrix[i] = org_dict[org]['dist asph avg']['10']
            f.write(' '.join(map(str, np.nanmean(path_matrix,0)))+"\n")
    for item in list[1:]:
            for i,org in enumerate(org_id):
                path_matrix[i] = org_dict[org]['dist asph avg'][item]
            with open(theme+"_avg_asph.dat", 'a') as f:
                f.write(' '.join(map(str, np.nanmean(path_matrix,0)))+"\n")
    with open(theme+"_std_asph.dat", 'w') as f:
            for i,org in enumerate(org_id):
                path_matrix[i] = org_dict[org]['dist asph avg']['10']
            f.write(' '.join(map(str, np.nanstd(path_matrix,0)))+"\n")
    for item in list[1:]:
            for i,org in enumerate(org_id):
                path_matrix[i] = org_dict[org]['dist asph avg'][item]
            with open(theme+"_std_asph.dat", 'a') as f:
                f.write(' '.join(map(str, np.nanstd(path_matrix,0)))+"\n")

##################################################
# Calculating tunnel rings                       #
# Which are rings build from the edge points     #
# That are perpendicular to the main tunnel path #
##################################################
def get_rings(path,org_dict,org_id,dot_threshold,ref_tunnel_path,short=False):
    atom_type = "C"
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]
    # Calculations for average map
    for name in tqdm(org_id):
        for item in list:
            # Getting tangent vectors
            tangents = np.diff(ref_tunnel_path.T,axis=0)
            # normalize tangents
            tangents_n = tangents / np.linalg.norm(tangents,axis=1)[:,np.newaxis]
            # distance along the spline_path
            n = np.shape(tangents)[0]
            tangents_len = np.sqrt(np.sum(tangents**2,1))
            path_l = np.zeros(n+1)
            path_l[:-1] = np.cumsum(tangents_len[::-1])[::-1]
            for i in range(0,n):
                # Vector from spline to all edge points
                spline_point = ref_tunnel_path.T[i]
                vec_to_edge = org_dict[name]['edges coord avg'][item] - spline_point
                # Calculate the dot product
                dot_products = np.dot(vec_to_edge, tangents_n[i])
                # Find indices where the absolute value of the dot product is less than the threshold
                perpendicular_indices = np.where(np.abs(dot_products) < dot_threshold)[0]
                if perpendicular_indices.size > 0 :
                    # Select the edge points that are close to being perpendicular to 'vektor'
                    perpendicular_edge_points = org_dict[name]['edges coord avg'][item][perpendicular_indices]
                    # Writing a ring
                    output_name = path+name+"/"+item+"/"+"ring_avg_"+str(np.round(path_l[i],2))+".pdb"
                    file = open(output_name,"w")
                    for r in range(0,len(perpendicular_edge_points)):
                        file.write(pdb_format.format(r+1, atom_type, r+1, perpendicular_edge_points[r][0], perpendicular_edge_points[r][1], perpendicular_edge_points[r][2], atom_type))
                    file.close()


##################################################
# Calculating tunnel discs                       #
# Which are disc build from occupancy map voxels #
# That are perpendicular to the main tunnel path #
##################################################
def cluster_points(perpendicular_interior_points):
    # Initialize cluster IDs and checked list
    n_points = len(perpendicular_interior_points)
    cluster_id = np.zeros(n_points)
    checked_list = []
    cluster_name = 1
    while any(cluster_id == 0):
        # Start with the first point as the seed for the first cluster
        sample_idx = np.where(cluster_id == 0)[0][0]
        cluster_id[sample_idx] = cluster_name
        checked_list.append(sample_idx)
        while checked_list:
            # Take the first unchecked point from the list
            sample_idx = checked_list.pop(0)
            sample = perpendicular_interior_points[sample_idx]
            # Find points within the threshold distance and assign them to the same cluster
            distances = np.linalg.norm(perpendicular_interior_points - sample, axis=1)
            neighbors_idx = np.where(distances < 2)[0] # as the voxel size is 1A it could be 1, sqrt(2) or sqrt(3) so 2 for simplicity
            for idx in neighbors_idx:
                if cluster_id[idx] == 0:  # If not yet assigned to any cluster
                    cluster_id[idx] = cluster_name
                    checked_list.append(idx)
        cluster_name +=1
    return cluster_id

def get_interior_edge(perpendicular_interior_points,edge_points):
    int_edge = []
    for i in edge_points:
        distance = np.linalg.norm(perpendicular_interior_points - i, axis=1)
        if np.min(distance)<2:
            int_edge.append(i)
    return np.array(int_edge)

def compute_cross_section_asphericity(disc):
    # Center the data
    centroid = disc.mean(axis=0)
    centered = disc - centroid
    # Perform PCA to find best-fit plane
    # Get eigenvectors of the covariance matrix
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # Plane normal is the eigenvector with smallest eigenvalue
    normal = eigenvectors[:, 0]
    # The other two eigenvectors define the in-plane directions
    basis_2D = eigenvectors[:, 1:]
    # Project points to 2D plane
    projected_2D = centered @ basis_2D 
    # Compute 2D gyration tensor
    gyration_tensor = np.cov(projected_2D.T)
    lambda1, lambda2 = np.linalg.eigvalsh(gyration_tensor)
    # Compute asphericity
    A = ((lambda1 - lambda2)**2) / (lambda1 + lambda2)**2
    return A


def get_discs(path,org_dict,org_id,dot_threshold,ref_tunnel_path,short=False):
    if short==False:
        list = ["10", "20", "30", "40", "60"]
    else: list = ["10", "20", "30", "40"]
    atom_type = "P"
    # Calculations for average map
    for name in tqdm(org_id):
        org_dict[name]['path length'] = dict()
        org_dict[name]['dist mean avg'] = dict()
        org_dict[name]['dist asph avg'] = dict()
        org_dict[name]['cross section'] = dict()
        org_dict[name]["proteins_matrix"] = dict()
        org_dict[name]["rna_matrix"] = dict()
        rib_path = path+name+"/"
        rib_file = glob.glob(os.path.join(rib_path,'*_tunnel_fix.pdb'))[0]
        ribosome = mda.Universe(rib_file)
        # Select all atoms in the ribosome
        all_atoms = ribosome.select_atoms("all")
        # Getting tangent vectors for the whole tunnel defined with item = 60
        tangents = np.diff(ref_tunnel_path.T,axis=0)
        # normalize tangents
        tangents_n = tangents / np.linalg.norm(tangents,axis=1)[:,np.newaxis]
        # distance along the spline_path
        n = np.shape(tangents)[0]
        tangents_len = np.sqrt(np.sum(tangents**2,1))
        path_l = np.zeros(n+1)
        path_l[:-1] = np.cumsum(tangents_len[::-1])[::-1]
        for item in list:
            org_dict[name]['path length'][item] = path_l
            org_dict[name]['dist mean avg'][item] = []
            org_dict[name]['dist asph avg'][item] = []
            org_dict[name]['cross section'][item] = []
            org_dict[name]["proteins_matrix"][item] = np.zeros([n,5])
            org_dict[name]["rna_matrix"][item] = np.zeros([n])
            for i in range(0,n):
                # Vector from splint to all perpendicular_interior_points points
                spline_point = ref_tunnel_path.T[i]
                vec_to_interior = org_dict[name]['interior coord avg'][item] - spline_point
                # Calculate the dot product
                dot_products = np.dot(vec_to_interior, tangents_n[i])
                # Find indices where the absolute value of the dot product is less than the threshold
                perpendicular_indices = np.where(np.abs(dot_products) < dot_threshold)[0]
                if perpendicular_indices.size > 0 :
                    # Select the interior points that are close to being perpendicular to 'vektor'
                    perpendicular_interior_points = org_dict[name]['interior coord avg'][item][perpendicular_indices]
                    # checking for clusters in perpendicular_interior_points
                    cluster_id = cluster_points(perpendicular_interior_points)
                    spline_mind_dist = np.min(np.linalg.norm(perpendicular_interior_points - spline_point, axis=1))
                    close_atoms = []
                    # testing if we even have a disc that is crossed by spline
                    if spline_mind_dist <= 8.0:
                        main_cluster_id = cluster_id[np.argmin(np.linalg.norm(perpendicular_interior_points - spline_point, axis=1))]
                        main_interior = np.where(cluster_id == main_cluster_id)[0]
                        edge_points = mda.Universe(path+name+"/"+item+"/"+"ring_avg_"+str(np.round(path_l[i],2))+".pdb").atoms.positions
                        # getting edge around the main interior cluster
                        interior_edge = get_interior_edge(perpendicular_interior_points[main_interior],edge_points)
                        # Construction of the whole disc
                        i0 = np.shape(perpendicular_interior_points[main_interior])[0]
                        i1 = np.shape(interior_edge)[0]
                        whole_disc = np.zeros([i0+i1,3])
                        whole_disc[:i0] = perpendicular_interior_points[main_interior]
                        whole_disc[i0:] = interior_edge
                        # Calculating properties of the disc
                        org_dict[name]['cross section'][item].append(main_interior.size+interior_edge.size/3.)
                        dist = np.linalg.norm(interior_edge - spline_point, axis=1)
                        org_dict[name]['dist mean avg'][item].append(np.mean(dist))
                        org_dict[name]['dist asph avg'][item].append(compute_cross_section_asphericity(whole_disc))
                        # Writing a disc
                        output_name = path+name+"/"+item+"/"+"disc_avg_"+str(np.round(path_l[i],2))+".pdb"
                        file = open(output_name,"w")
                        atom_type = "C"
                        for r in range(0,i0):
                            file.write(pdb_format.format(r+1, atom_type, 1, whole_disc[r][0], whole_disc[r][1], whole_disc[r][2], atom_type))
                        for r in range(i0,len(whole_disc)):
                            file.write(pdb_format.format(r+1, atom_type, 1, whole_disc[r][0], whole_disc[r][1], whole_disc[r][2], atom_type))
                        file.close()                 
                        # Find ribosome atoms within 5A of the given coordinates
                        # Threshold distance
                        threshold = 5.0
                        for coord in interior_edge:
                            distances = distance_array(coord, all_atoms.positions)
                            within_threshold = np.where(distances <= threshold)[1]
                            for atom_index in within_threshold:
                                close_atoms.append(all_atoms[atom_index])
                        uniq_atoms = np.unique(close_atoms)
                        if len(uniq_atoms) > 0:
                            # Get unique residues containing these atoms
                            unique_residues_id = np.unique([atom.resindex for atom in uniq_atoms])
                            # Select the unique residues
                            unique_residues = ribosome.residues[unique_residues_id]
                            # Writing uniue residue
                            output = rib_path + item + "/tunnel/"
                            if not os.path.exists(output):
                                os.makedirs(output)
                            unique_residues.atoms.write(output+str(np.round(path_l[i],2))+".pdb")
                            # Get all atoms in the unique residues 
                            unique_atoms = unique_residues.atoms
                            # Select protein and RNA atoms from the unique atoms
                            protein_atoms = unique_atoms.select_atoms('resname ALA GLY PRO CYS TRP TYR ILE LEU LYS VAL GLN GLU ASP ASN ASP THR PHE HIS SER ARG')
                            rna_atoms = unique_atoms.select_atoms('not resname ALA GLY PRO CYS TRP TYR ILE LEU LYS VAL GLN GLU ASP ASN ASP THR PHE HIS SER ARG')
                            protein_resnames = protein_atoms.residues.resnames
                            if protein_atoms.n_atoms == 0:
                                protein_atoms = unique_atoms.select_atoms('resname ALA? GLY? PRO? CYS? TRP? TYR? ILE? LEU? LYS? VAL? GLN? GLU? ASP? ASN? ASP? THR? PHE? HIS? SER? ARG?')
                                rna_atoms = unique_atoms.select_atoms('not resname ALA? GLY? PRO? CYS? TRP? TYR? ILE? LEU? LYS? VAL? GLN? GLU? ASP? ASN? ASP? THR? PHE? HIS? SER? ARG?')
                                protein_resnames = np.array([x[:3] for x in protein_atoms.residues.resnames])
                            polar_n  = np.sum(np.isin(protein_resnames,polar))
                            hydrophobic_n = np.sum(np.isin(protein_resnames,hydrophobic))
                            positive_n = np.sum(np.isin(protein_resnames,positive))
                            negative_n = np.sum(np.isin(protein_resnames,negative))
                            special_n = np.sum(np.isin(protein_resnames,special))
                            org_dict[name]["proteins_matrix"][item][i,] = [polar_n,hydrophobic_n,positive_n,negative_n,special_n]
                            org_dict[name]["rna_matrix"][item][i] = len(rna_atoms.residues)
                    else:
                        org_dict[name]['dist mean avg'][item].append(np.nan)
                        org_dict[name]['dist asph avg'][item].append(np.nan)
                        org_dict[name]['cross section'][item].append(np.nan)
                else:
                    org_dict[name]['dist mean avg'][item].append(np.nan)
                    org_dict[name]['dist asph avg'][item].append(np.nan)
                    org_dict[name]['cross section'][item].append(np.nan)


def get_discs_without_clusters(path,org_dict,org_id,dot_threshold,ref_org):
    atom_type = "P"
    # Calculations for average map
    for name in tqdm(org_id):
        org_dict[name]['path length'] = dict()
        org_dict[name]['dist mean avg'] = dict()
        org_dict[name]['dist asph avg'] = dict()
        org_dict[name]['cross section'] = dict()
        org_dict[name]["proteins_matrix"] = dict()
        org_dict[name]["rna_matrix"] = dict()
        rib_path = path+name+"/"
        rib_file = glob.glob(os.path.join(rib_path,'*_tunnel_fix.pdb'))[0]
        ribosome = mda.Universe(rib_file)
        # Select all atoms in the ribosome
        all_atoms = ribosome.select_atoms("all")
        # Getting tangent vectors for the whole tunnel defined with item = 60
        tangents = np.diff(np.array(euk_dict[ref_org]['spline path avg']["60"]).T,axis=0)
        # normalize tangents
        tangents_n = tangents / np.linalg.norm(tangents,axis=1)[:,np.newaxis]
        # distance along the spline_path
        n = np.shape(tangents)[0]
        tangents_len = np.sqrt(np.sum(tangents**2,1))
        path_l = np.zeros(n+1)
        path_l[:-1] = np.cumsum(tangents_len[::-1])[::-1]
        for item in ["10", "20", "30", "40", "60"]:
            org_dict[name]['path length'][item] = path_l
            org_dict[name]['dist mean avg'][item] = []
            org_dict[name]['dist asph avg'][item] = []
            org_dict[name]['cross section'][item] = []
            org_dict[name]["proteins_matrix"][item] = np.zeros([n,5])
            org_dict[name]["rna_matrix"][item] = np.zeros([n])
            for i in range(0,n):
                # Vector from splint to all perpendicular_interior_points points
                spline_point = np.array(euk_dict[ref_org]['spline path avg']["60"]).T[i]
                vec_to_interior = org_dict[name]['interior coord avg'][item] - spline_point
                # Calculate the dot product
                dot_products = np.dot(vec_to_interior, tangents_n[i])
                # Find indices where the absolute value of the dot product is less than the threshold
                perpendicular_indices = np.where(np.abs(dot_products) < dot_threshold)[0]
                if perpendicular_indices.size > 0 :
                    # Select the interior points that are close to being perpendicular to 'vektor'
                    perpendicular_interior_points = org_dict[name]['interior coord avg'][item][perpendicular_indices]
                    close_atoms = []
                    # testing if we even have a disc that is crossed by spline
                    edge_points = mda.Universe(path+name+"/"+item+"/"+"ring_avg_"+str(np.round(path_l[i],2))+".pdb").atoms.positions
                    # getting edge around the main interior cluster
                    interior_edge = get_interior_edge(perpendicular_interior_points,edge_points)
                    # Construction of the whole disc
                    i0 = np.shape(perpendicular_interior_points)[0]
                    i1 = np.shape(interior_edge)[0]
                    whole_disc = np.zeros([i0+i1,3])
                    whole_disc[:i0] = perpendicular_interior_points
                    whole_disc[i0:] = interior_edge
                    # Calculating properties of the disc
                    org_dict[name]['cross section'][item].append(perpendicular_interior_points.size/3+interior_edge.size/3.)
                    dist = np.linalg.norm(interior_edge - spline_point, axis=1)
                    org_dict[name]['dist mean avg'][item].append(np.mean(dist))
                    org_dict[name]['dist asph avg'][item].append(np.min(dist)/np.max(dist))
                    # Writing a disc
                    output_name = path+name+"/"+item+"/"+"disc_avg_"+str(np.round(path_l[i],2))+".pdb"
                    file = open(output_name,"w")
                    atom_type = "C"
                    for r in range(0,i0):
                        file.write(pdb_format.format(r+1, atom_type, whole_disc[r][0], whole_disc[r][1], whole_disc[r][2], atom_type))
                    for r in range(i0,len(whole_disc)):
                        file.write(pdb_format.format(r+1, atom_type, whole_disc[r][0], whole_disc[r][1], whole_disc[r][2], atom_type))
                    file.close()                 
                    # Find ribosome atoms within 5A of the given coordinates
                    # Threshold distance
                    threshold = 5.0
                    for coord in interior_edge:
                        distances = distance_array(coord, all_atoms.positions)
                        within_threshold = np.where(distances <= threshold)[1]
                        for atom_index in within_threshold:
                            close_atoms.append(all_atoms[atom_index])
                    uniq_atoms = np.unique(close_atoms)
                    if len(uniq_atoms) > 0:
                        # Get unique residues containing these atoms
                        unique_residues_id = np.unique([atom.resindex for atom in uniq_atoms])
                        # Select the unique residues
                        unique_residues = ribosome.residues[unique_residues_id]
                        # Writing uniue residue
                        output = rib_path + item + "/tunnel/"
                        if not os.path.exists(output):
                            os.makedirs(output)
                        unique_residues.atoms.write(output+str(np.round(path_l[i],2))+".pdb")
                        # Get all atoms in the unique residues 
                        unique_atoms = unique_residues.atoms
                        # Select protein and RNA atoms from the unique atoms
                        protein_atoms = unique_atoms.select_atoms('resname ALA GLY PRO CYS TRP TYR ILE LEU LYS VAL GLN GLU ASP ASN ASP THR PHE HIS SER ARG')
                        rna_atoms = unique_atoms.select_atoms('not resname ALA GLY PRO CYS TRP TYR ILE LEU LYS VAL GLN GLU ASP ASN ASP THR PHE HIS SER ARG')
                        protein_resnames = protein_atoms.residues.resnames
                        if protein_atoms.n_atoms == 0:
                            protein_atoms = unique_atoms.select_atoms('resname ALA? GLY? PRO? CYS? TRP? TYR? ILE? LEU? LYS? VAL? GLN? GLU? ASP? ASN? ASP? THR? PHE? HIS? SER? ARG?')
                            rna_atoms = unique_atoms.select_atoms('not resname ALA? GLY? PRO? CYS? TRP? TYR? ILE? LEU? LYS? VAL? GLN? GLU? ASP? ASN? ASP? THR? PHE? HIS? SER? ARG?')
                            protein_resnames = np.array([x[:3] for x in protein_atoms.residues.resnames])
                        polar_n  = np.sum(np.isin(protein_resnames,polar))
                        hydrophobic_n = np.sum(np.isin(protein_resnames,hydrophobic))
                        positive_n = np.sum(np.isin(protein_resnames,positive))
                        negative_n = np.sum(np.isin(protein_resnames,negative))
                        special_n = np.sum(np.isin(protein_resnames,special))
                        org_dict[name]["proteins_matrix"][item][i,] = [polar_n,hydrophobic_n,positive_n,negative_n,special_n]
                        org_dict[name]["rna_matrix"][item][i] = len(rna_atoms.residues)
                    else:
                        org_dict[name]['dist mean avg'][item].append(np.nan)
                        org_dict[name]['dist asph avg'][item].append(np.nan)
                        org_dict[name]['cross section'][item].append(np.nan)
                else:
                    org_dict[name]['dist mean avg'][item].append(np.nan)
                    org_dict[name]['dist asph avg'][item].append(np.nan)
                    org_dict[name]['cross section'][item].append(np.nan)


#############################################################
# Function to interpolate cross-sections on the common path #
#############################################################
def path_interp(org_dict,org_id):
    interp_cs = dict()
    for item in ["10", "20", "30", "40","60"]:
        average_p = []
        length_p = []
        for name in org_id:
            average_p.append(org_dict[name]['cross section'][item])
            length_p.append(org_dict[name]['path length'][item][:-1])
            len_max = int(np.ceil(np.max(length_p)))
        new_length = np.arange(0,len_max+1)
        interp_cs[item] = np.zeros([len(org_id),len_max+1])
        # Interpolate cross sections on new_length
        for i in range(0,len(org_id)):
            mask = ~np.isnan(average_p[i])
            interp_func1 = interp1d(length_p[i][mask], np.array(average_p[i])[mask], kind='linear', fill_value="extrapolate")
            y_interp = interp_func1(new_length)
            y_min = int(np.floor(length_p[i][mask][-1]))
            y_max = int(np.ceil(length_p[i][mask][0]))
            y_interp[:y_min] = np.nan
            y_interp[y_max:] = np.nan
            interp_cs[item][i] = y_interp
    return interp_cs,new_length

def path_interp_prop(org_dict,org_id,length):
    item = str(length)
    polar = []
    hydro = [] 
    negat = []
    posit = []
    speci = []
    rna = []
    length_p = []
    for name in org_id:
        polar.append(org_dict[name]["proteins_matrix"][item][:,0])
        hydro.append(org_dict[name]["proteins_matrix"][item][:,1])
        posit.append(org_dict[name]["proteins_matrix"][item][:,2])
        negat.append(org_dict[name]["proteins_matrix"][item][:,3])
        speci.append(org_dict[name]["proteins_matrix"][item][:,4])
        rna.append(org_dict[name]["rna_matrix"][item])
        length_p.append(org_dict[name]['path length'][item][:-1])
        len_max = int(np.ceil(np.max(length_p)))
    new_length = np.arange(0,len_max+1)
    interp_polar = np.zeros([len(org_id),len_max+1])
    interp_hydro = np.zeros([len(org_id),len_max+1])
    interp_negat = np.zeros([len(org_id),len_max+1])
    interp_posit = np.zeros([len(org_id),len_max+1])
    interp_speci = np.zeros([len(org_id),len_max+1])
    interp_rna = np.zeros([len(org_id),len_max+1])
    # Interpolate various properties on the on new_length
    for i in range(0,len(org_id)):
        interp_func1 = interp1d(length_p[i], np.array(polar[i]), kind='linear', fill_value="extrapolate")
        interp_func2 = interp1d(length_p[i], np.array(hydro[i]), kind='linear', fill_value="extrapolate")
        interp_func3 = interp1d(length_p[i], np.array(negat[i]), kind='linear', fill_value="extrapolate")
        interp_func4 = interp1d(length_p[i], np.array(posit[i]), kind='linear', fill_value="extrapolate")
        interp_func5 = interp1d(length_p[i], np.array(speci[i]), kind='linear', fill_value="extrapolate")
        interp_func6 = interp1d(length_p[i], np.array(rna[i]), kind='linear', fill_value="extrapolate")
        interp_polar[i] = interp_func1(new_length)
        interp_hydro[i] = interp_func2(new_length)
        interp_negat[i] = interp_func3(new_length)
        interp_posit[i] = interp_func4(new_length)
        interp_speci[i] = interp_func5(new_length)
        interp_rna[i] = interp_func6(new_length)
        interp_polar[i][interp_polar[i]<0] = 0
        interp_hydro[i][interp_hydro[i]<0] = 0
        interp_negat[i][interp_negat[i]<0] = 0
        interp_posit[i][interp_posit[i]<0] = 0
        interp_speci[i][interp_speci[i]<0] = 0
        interp_rna[i][interp_rna[i]<0] = 0
    return interp_polar,interp_hydro,interp_posit,interp_negat,interp_speci,interp_rna,new_length

import numpy as np
import matplotlib.pyplot as plt
import mrcfile
import MDAnalysis as mda
from tqdm import tqdm
import os

##############################
##### Ribosome interface #####
##############################

# First we need to find all atoms of the ribosome that where in the contact with NC so they build the tunnelimport numpy as n

# Load trajectory and ribosome structure
traj = mda.Universe('40/NC_FME_40.pdb', '40/fitted.xtc')
rib = mda.Universe('6v39_tunnel_fix.pdb')

# Parameters
N = traj.trajectory.n_frames
cut_off = 3.0

# Get ribosome atom positions
rib_positions = rib.atoms.positions

# Set to store ribosome indexes in contact with the protein
rib_contact_indices = set()

# Preselect atoms for each residue once, not every time in the loop
residue_atom_groups = [residue.atoms for residue in traj.residues]

# Loop over the frames with a progress bar
for i in tqdm(range(N), desc='Processing frames'):
    traj.trajectory[i]  # update the current frame
    for atom_group in residue_atom_groups:
        # Calculate the distance from each atom in the residue to all atoms in the ribosome
        distances = mda.lib.distances.distance_array(atom_group.positions, rib_positions)
        # Check if any atoms are within the cutoff distance and update the set
        if np.any(distances.min(axis=1) < cut_off):
            rib_contact_indices.update(np.where(distances < cut_off)[1])

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import matplotlib.pyplot as plt

# Define your lists of IDs
bac_id = ['abau','bbur','bsub','drad','ecoli','efae','fjoh','linn','llac','lmon','mpne','msme','mtub','paer', 'pgin','pura','saur', 'sfra','tthe', 'vnat']
euk_id = ['atha','calb','cele', 'cpor', 'cser', 'cthe','dmel','drer','ecun','egra','ehis','ggal','glam','hgla','hsap','klac','ldon','lmaj','mmus','ncra','ntab','ocun','pfal','ploc','rnor','scer','slop','slyc','spom','sscr','taes','tbru','tcru','tgon','tthe','tvag','vnec','xlae']
arc_id = ['hmar', 'mace','pfur','saci','tkod','pcal']
#chloro_id = ['chloro2']

# 1. Master dictionary to eliminate repetitive if/elif loops!
domain_setup = {
    "bac": {"ids": bac_id, "path": "/home/twlodarski/Projects/ribosome_tunnels/new_code/data_new/bacteria/"},
    "arc": {"ids": arc_id, "path": "/home/twlodarski/Projects/ribosome_tunnels/new_code/data_new/archaea/"},
    "euk": {"ids": euk_id, "path": "/home/twlodarski/Projects/ribosome_tunnels/new_code/data_new/eukaryota/"}
    #"chloro": {"ids": chloro_id, "path": "/home/twlodarski/Projects/ribosome_tunnels/chloro/"}
}

# Prepare a dictionary to hold exit values
exit_lengths = {key: [] for key in domain_setup}

# 2. The Unified Loop
for domain, info in domain_setup.items():
    base_path = info["path"]
    for name in info["ids"]:
        # Safe path joining
        disc_file = os.path.join(base_path, name, "60", "disc_avg_91.72.pdb")
        spline_file = os.path.join(base_path, name, "60", "spline_path_avg.pdb")
        disc = mda.Universe(disc_file)
        spline = mda.Universe(spline_file)
        disc_coords = disc.atoms.positions
        path_coords = spline.atoms.positions
        # Compute distance array and find closest points
        D = distance_array(disc_coords, path_coords)
        disc_atom_index, path_atom_index = np.unravel_index(np.argmin(D), D.shape)
        # 3. CRITICAL BUG FIX: Ensure cum_path has the exact same length as path_coords
        n_points = len(path_coords)
        # Fast numpy vector length calculation
        path_l = np.linalg.norm(np.diff(path_coords, axis=0), axis=1)
        # Initialize an array of exact size to prevent IndexErrors
        cum_path = np.zeros(n_points)
        cum_path[:-1] = np.cumsum(path_l[::-1])[::-1]
        # Store the value
        exit_lengths[domain].append(cum_path[path_atom_index])

# ---------------------------------------------------------
# Plotting
# ---------------------------------------------------------
data = [exit_lengths["bac"], exit_lengths["arc"], exit_lengths["euk"]]
labels = ["Bacteria", "Archaea", "Eukaryota"]
colors = ["#D36B8B", "#D4AF37", "#66B2B2"]

plt.figure(figsize=(8, 6))

# Create the boxplot with slightly transparent boxes
bp = plt.boxplot(data, labels=labels, patch_artist=True, showmeans=False, showfliers=False,
                 medianprops=dict(color='black', linewidth=2))

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6) # Transparency lets the scatter points stand out

# Add individual data points on top (Jitter Plot)
for i, d in enumerate(data):
    # np.random.normal adds a tiny bit of horizontal scatter so points don't overlap perfectly
    x_jitter = np.random.normal(i + 1, 0.05, size=len(d))
    plt.scatter(x_jitter, d, color=colors[i], edgecolors='black', linewidth=0.5, zorder=3, alpha=0.6,s=50)

# Formatting
plt.ylabel("Length of the ribosome exit tunnel [Å]", fontsize=16)
#plt.title("Ribosome Exit Tunnel Lengths Across Domains", fontsize=16, pad=15)

plt.xticks(fontsize=16)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig("tunnel_lengths_boxplot.svg")
plt.show()
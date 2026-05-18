import numpy as np
import matplotlib.pyplot as plt

# Load data
b = np.loadtxt("bacteria_mean_composition.dat")
e = np.loadtxt("eukaryota_core_mean_composition.dat")
a = np.loadtxt("archaea_mean_composition.dat")
b_std = np.loadtxt("bacteria_std_composition.dat")
e_std = np.loadtxt("eukaryota_core_std_composition.dat")
a_std = np.loadtxt("archaea_std_composition.dat")

labels = ['Bacteria', 'Eukaryota', 'Archaea'] # Capitalized for aesthetics
features = ["Polar", "Hydrophobic", "Positive", "Negative", "Special"]
colors = ["pink", "orange", "red", "blue", "silver"]

num_features = len(features)
w = 0.15
x = np.arange(3)  

domain_names = ["L=10 amino acids", "L=20 amino acids", "L=30 amino acids", "L=40 amino acids", "L=60 amino acids"]

fig, axs = plt.subplots(2, 3, figsize=(14, 8)) # Made slightly larger to fit the legend comfortably
axs = axs.flatten()  

for j in range(5):  
    ax = axs[j]
    # Add a horizontal grid BEHIND the bars
    ax.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
    for i in range(num_features):
        offset = (i - num_features / 2) * w + w / 2
        bar_positions = x + offset
        heights = [b[j][i], e[j][i], a[j][i]]
        stds    = [b_std[j][i], e_std[j][i], a_std[j][i]]
        # Only assign labels on the very first subplot to avoid duplicates
        bar_label = features[i] if j == 0 else None
        ax.bar(bar_positions, heights, width=w, color=colors[i],
               edgecolor="black", linewidth=1, yerr=stds, capsize=3, 
               label=bar_label, zorder=3) # zorder=3 puts bars on top of the grid
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    # Only add the Y-label to the leftmost plots to reduce text clutter
    if j == 0 or j == 3:
        ax.set_ylabel("Number of Amino Acids", fontsize=12)
    ax.set_title(domain_names[j], fontsize=14, pad=10)
    ax.set_ylim(0, 35)

# Utilize the empty 6th subplot purely for the legend!
axs[5].axis('off')  
# Extract the handles and labels from the first subplot to build the legend
handles, leg_labels = axs[0].get_legend_handles_labels()
axs[5].legend(handles, leg_labels, loc='center', fontsize=14, frameon=False)

plt.tight_layout()
plt.savefig("domain_composition_comparison.svg")
plt.show()

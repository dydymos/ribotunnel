import numpy as np
import matplotlib.pyplot as plt

# Load length data safely
end = 91.71631884486219
length = np.loadtxt("length.dat")[:-1]
cut = np.argmin(np.abs(length - end)) # Crash-proof closest match

# Load ratio and std data
b_ratio = np.loadtxt("bacteria_tunnel_composition_ratio.dat")
a_ratio = np.loadtxt("archaea_tunnel_composition_ratio.dat")
e_ratio = np.loadtxt("eukaryota_core_tunnel_composition_ratio.dat")

b_ratio_std = np.loadtxt("bacteria_tunnel_composition_ratio_std.dat")
a_ratio_std = np.loadtxt("archaea_tunnel_composition_ratio_std.dat")
e_ratio_std = np.loadtxt("eukaryota_core_tunnel_composition_ratio_std.dat")

# Load Mitochondria data (assuming no std file based on your script)
#mito_ratio = np.loadtxt("mitochondria/crei/mitochondria_tunnel_composition_ratio.dat")

# Softer, distinct colors
colors = {
    "Bacteria": "#D36B8B",    # light burgundy rose
    "Eukaryota": "#66B2B2",   # light teal
    "Archaea": "#D4AF37",     # light gold
#    "Mitochondria": "#8E44AD" # soft purple (for the loaded mito data)
}

# Pack data into a list for clean plotting [Label, Mean, Std, LineStyle]
plot_data = [
    ("Bacteria", b_ratio, b_ratio_std, "-"),
    ("Eukaryota", e_ratio, e_ratio_std, "-"),
    ("Archaea", a_ratio, a_ratio_std, "-")
]

plt.figure(figsize=(8, 6))

# Plot the 50/50 equilibrium line in the background
plt.axhline(0.5, color="black", ls='--', lw=2, zorder=1)

# X-axis values starting from the cut point
x_vals = length[cut:]

# Loop through the main domains and plot them automatically
for label, mean_val, std_val, ls in plot_data:
    plt.plot(x_vals[:-1], mean_val[cut:][:-1], label=label, color=colors[label], linewidth=3, linestyle=ls, zorder=3)
    plt.fill_between(x_vals[:-1], mean_val[cut:][:-1] - std_val[cut:][:-1], mean_val[cut:][:-1] + std_val[cut:][:-1], 
                     color=colors[label], alpha=0.3, zorder=2)

# Plot Mitochondria (no standard deviation shading since it's one organism)
# Uncomment this if you want to see the mito data you loaded!
# plt.plot(x_vals, mito_ratio[cut:], label="Mitochondria (crei)", color=colors["Mitochondria"], linewidth=3, linestyle=':', zorder=4)

# Formatting
plt.xlabel("Distance from the PTC (Å)", fontsize=16) 
plt.ylabel("RNA vs Protein ratio", fontsize=16)

# Force the legend to match the size of the larger axis labels
plt.legend(fontsize=12, loc="upper right") 

plt.ylim(0, 1.1)
plt.xlim(0, 100)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig("rna_protein_ratio.svg")
plt.show()
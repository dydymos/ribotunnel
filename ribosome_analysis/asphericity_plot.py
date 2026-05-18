import numpy as np
import matplotlib.pyplot as plt

# 1. Load length and calculate cut safely
end = 91.71631884486219
length = np.loadtxt("length.dat")[:-1]
cut = np.argmin(np.abs(length - end))

# 2. Load path asphericity data
b_path = np.loadtxt("bacteria_avg_asph.dat")
b_std = np.loadtxt("bacteria_std_asph.dat")
a_path = np.loadtxt("archaea_avg_asph.dat")
a_std = np.loadtxt("archaea_std_asph.dat")
e_path = np.loadtxt("eukaryota_core_avg_asph.dat")
e_std = np.loadtxt("eukaryota_core_std_asph.dat")

# (Optional: Uncomment if you generated Microsporidia asphericity data!)
# m_path = np.loadtxt("micro_avg_asph.dat")
# m_std = np.loadtxt("micro_std_asph.dat")

# Formatting dictionaries
colors = {
    "Bacteria": "#D36B8B",    
    "Eukaryota": "#66B2B2",   
    "Archaea": "#D4AF37",
    # "Microsporidia": "#2C3E50"
}

titles = ["L=10 amino acids", "L=20 amino acids", "L=30 amino acids", "L=40 amino acids", "L=60 amino acids"]

# Package data to cleanly loop through the domains
domains = [
    ("Bacteria", b_path, b_std),
    ("Eukaryota", e_path, e_std),
    ("Archaea", a_path, a_std)
    # ("Microsporidia", m_path, m_std)
]

fig, axs = plt.subplots(2, 3, figsize=(14, 8))
axs = axs.flatten()
x_vals = length[cut:]

# --- Subplots 1 to 5 (Asphericity) ---
for i in range(5):
    ax = axs[i]
    for label, path_data, std_data in domains:
        y_mean = path_data[i][cut:]
        y_std = std_data[i][cut:]
        
        ax.plot(x_vals, y_mean, label=label, color=colors[label], linewidth=3, zorder=3)
        ax.fill_between(x_vals, y_mean - y_std, y_mean + y_std, color=colors[label], alpha=0.2, zorder=2)
        
    ax.set_xlabel("Distance from the PTC [Å]")
    ax.set_ylabel("Tunnel Asphericity")
    ax.set_title(titles[i], fontsize=14)
    ax.set_ylim(0, 1)

# --- Subplot 6 (The Legend Canvas) ---
axs[5].axis('off') # Hide the empty axes
handles, labels = axs[0].get_legend_handles_labels()
axs[5].legend(handles, labels, loc='center', fontsize=14, frameon=False)

# Final Rendering
plt.tight_layout()
plt.savefig("asphericity_analysis.svg")
plt.show()
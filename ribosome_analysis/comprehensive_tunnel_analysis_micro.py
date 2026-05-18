import numpy as np
import matplotlib.pyplot as plt

# 1. Load length and calculate cut safely
end = 91.71631884486219
length = np.loadtxt("length.dat")[:-1]
cut = np.argmin(np.abs(length - end)) 

# 2. Load path cross-section data
b_path = np.loadtxt("bacteria_avg_path.dat")
b_std = np.loadtxt("bacteria_std_path.dat")
a_path = np.loadtxt("archaea_avg_path.dat")
a_std = np.loadtxt("archaea_std_path.dat")
e_path = np.loadtxt("eukaryota_core_avg_path.dat")
e_std = np.loadtxt("eukaryota_core_std_path.dat")
m_path = np.loadtxt("micro_avg_path.dat")
m_std = np.loadtxt("micro_std_path.dat")

# 3. Load ratio data 
b_ratio = np.loadtxt("bacteria_tunnel_composition_ratio.dat")
b_ratio_std = np.loadtxt("bacteria_tunnel_composition_ratio_std.dat")
a_ratio = np.loadtxt("archaea_tunnel_composition_ratio.dat")
a_ratio_std = np.loadtxt("archaea_tunnel_composition_ratio_std.dat")
e_ratio = np.loadtxt("eukaryota_core_tunnel_composition_ratio.dat")
e_ratio_std = np.loadtxt("eukaryota_core_tunnel_composition_ratio_std.dat")
m_ratio = np.loadtxt("micro_tunnel_composition_ratio.dat")
m_ratio_std = np.loadtxt("micro_tunnel_composition_ratio_std.dat")

# Formatting dictionaries (Added Dark Slate for Microsporidia)
colors = {
    "Bacteria": "#D36B8B",    
    "Eukaryota": "#66B2B2",   
    "Archaea": "#D4AF37",  
    "Microsporidia": "#2C3E50" # Dark, sharp contrast for the outlier
}
titles = ["L=10 amino acids", "L=20 amino acids", "L=30 amino acids", "L=40 amino acids", "L=60 amino acids"]

# Package data to cleanly loop through ALL domains
domains = [
    ("Bacteria", b_path, b_std, b_ratio, b_ratio_std),
    ("Eukaryota", e_path, e_std, e_ratio, e_ratio_std),
    ("Archaea", a_path, a_std, a_ratio, a_ratio_std),
    ("Microsporidia", m_path, m_std, m_ratio, m_ratio_std) # Simply added to the list!
]

plt.figure(figsize=(16, 9))
x_vals = length[cut:]

# --- Subplots 1 to 5 (Cross-sections) ---
for i in range(5):
    plt.subplot(2, 3, i + 1)
    for label, path_data, std_data, _, _ in domains:
        y_mean = path_data[i][cut:]
        y_std = std_data[i][cut:]
        
        # We drop the alpha of the fill_between slightly more so 4 overlapping standard deviations don't block the plot
        plt.plot(x_vals, y_mean, label=label, color=colors[label], linewidth=3, zorder=3)
        plt.fill_between(x_vals, y_mean - y_std, y_mean + y_std, color=colors[label], alpha=0.2, zorder=2)
        
    plt.xlabel("Distance from the PTC [Å]")
    plt.ylabel("Cross-section [Å²]")
    plt.title(titles[i], fontsize=14)
    plt.ylim(0, 700 if i < 3 else 1500)
    if i == 0: plt.legend(loc="upper left")

# --- Subplot 6 (RNA vs Protein Ratio) ---
plt.subplot(2, 3, 6)
plt.axhline(0.5, color="black", ls='--', lw=2, zorder=1)

valid_mask = x_vals > 0.001 
x_vals_clean = x_vals[valid_mask]

for label, _, _, ratio_data, ratio_std_data in domains:
    y_mean = ratio_data[cut:][valid_mask]
    y_std = ratio_std_data[cut:][valid_mask]    
    plt.plot(x_vals_clean[:-1], y_mean[:-1], label=label, color=colors[label], linewidth=3, zorder=3)
    plt.fill_between(x_vals_clean[:-1], y_mean[:-1] - y_std[:-1], y_mean[:-1] + y_std[:-1], color=colors[label], alpha=0.3, zorder=2)

plt.xlabel("Distance from the PTC [Å]")
plt.ylabel("RNA vs Protein ratio")
plt.title("Chemical Composition", fontsize=14)
plt.ylim(0, 1.1)

# Final Rendering
plt.tight_layout()
plt.savefig("comprehensive_tunnel_analysis.svg")
plt.show()
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

################################
# Plotting 3D plot of pathways #
################################

def paths_plot3D(org_dict, name, item):
    """
    Plots a 3D visualization of individual nascent chain pathways,
    their splines, and the averaged pathway. Forces a cubic bounding box
    to prevent visual distortion of the tunnel geometry.
    """
    fig = plt.figure(figsize=(8, 8))
    # Create short local variables so we don't have to type the massive dictionary path every time
    res_paths = org_dict[name]['res occ path'][item]
    spline_paths = org_dict[name]['spline path'][item]
    avg_res_path = np.array(org_dict[name]['res occ path avg'][item])
    avg_spline = org_dict[name]['spline path avg'][item]
    # Using the 3rd path ('2') to determine the bounding box dimensions
    path_coords = np.array(res_paths['2'])
    x, y, z = path_coords.T
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
    # Use a loop to plot the 3 replicates to keep the code DRY (Don't Repeat Yourself)
    labels = {'0': '1st', '1': '2nd', '2': '3rd'}
    for i in ['0', '1', '2']:
        coords = np.array(res_paths[i])
        spline = spline_paths[i]
        # Scatter for raw positions, Line for spline
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], alpha=0.5) 
        ax.plot(spline[0], spline[1], spline[2], lw=2, label=f'The {labels[i]} path - spline')
    # Plot the Average Path in black
    ax.scatter(avg_res_path[:, 0], avg_res_path[:, 1], avg_res_path[:, 2], color="black")
    ax.plot(avg_spline[0], avg_spline[1], avg_spline[2], lw=3, color="black", label='The average path - spline')
    # Set the calculated limits for each axis
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_zlim(z_limits)
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    plt.legend()
    # Using an f-string for cleaner string concatenation
    plt.title(f"Ribosome from {org_dict[name]['name']} at {item} length")
    plt.show()

###################################################
# Function to plot total properties of the tunnel #
###################################################

def plot_tunnel_prop(org_dict, org_ids, output_name, short=False):
    """Generates violin/strip plots for tunnel physico-chemical properties."""
    chain_lengths = ["10", "20", "30", "40"] if short else ["10", "20", "30", "40", "60"]
    fig = plt.figure(figsize=(15, 10))
    categories = ["polar", "hydrophobic", "positive", "negative", "special"]
    colors = ["pink", "orange", "red", "blue", "silver"]
    # enumerate(..., 1) automatically gives us j=1, j=2, etc. for the subplots
    for j, length in enumerate(chain_lengths, 1):
        plt.subplot(2, 3, j)
        # Extract data for this chain length across all organisms into a clean 2D array
        data_matrix = np.array([org_dict[name]["proteins"][length] for name in org_ids])
        # Loop through the 5 categories to keep the plotting code DRY
        for col_idx, (cat, color) in enumerate(zip(categories, colors)):
            y_data = data_matrix[:, col_idx]
            x_labels = np.full(len(org_ids), cat)
            sns.violinplot(x=x_labels, y=y_data, color=color, edgecolor=color, inner=None, linewidth=1.5, alpha=0.2)
            sns.stripplot(x=x_labels, y=y_data, color=color, size=6, jitter=True)
        plt.ylim(-1, 45)
        plt.ylabel("Number of AA",fontsize=18)
        plt.title(f"NC length {length} AA",fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
    plt.tight_layout()  # Automatically adjusts spacing so labels don't overlap!
    plt.savefig(output_name)
    plt.close()




#######################################################
# Function to extract tunnel properties for plotting  #
#######################################################

def path_prop(org_dict, org_ids, short=False):
    """
    Extracts the physico-chemical matrices and path lengths for a given 
    set of organisms to prepare them for distance-based plotting.
    """
    polar, hydro, posit, negat, speci, rna, length_p = [], [], [], [], [], [], []
    # Dynamically select the chain length key based on the 'short' flag
    key = '40' if short else '60'
    for name in org_ids:
        # Create a local reference to the matrix to keep lines clean
        prot_matrix = org_dict[name]["proteins_matrix"][key]
        polar.append(prot_matrix[:, 0])
        hydro.append(prot_matrix[:, 1])
        posit.append(prot_matrix[:, 2])
        negat.append(prot_matrix[:, 3])
        speci.append(prot_matrix[:, 4])
        rna.append(org_dict[name]["rna_matrix"][key])
        # [:-1] drops the final zero to match the N-1 tangent vectors
        length_p.append(org_dict[name]['path length'][key][:-1])
    return polar, hydro, posit, negat, speci, rna, length_p

###########################################################
# Plot properties as a function of distance from the PTC  #
###########################################################

def plot_prop_dist(org_dict, org_ids, output_name, end=None, short=False):
    """
    Plots the physico-chemical composition and RNA/Protein ratio 
    along the length of the ribosomal exit tunnel.
    """
    # Note: Assuming path_prop is defined elsewhere in your file!
    polar, hydro, posit, negat, speci, rna, length = path_prop(org_dict, org_ids, short=short)
    # Convert lists to NumPy arrays for fast element-wise math
    polar = np.array(polar)
    hydro = np.array(hydro)
    posit = np.array(posit)
    negat = np.array(negat)
    speci = np.array(speci)
    rna = np.array(rna)
    length = np.array(length)
    # Calculate the total protein count (Summing along the property axis)
    protein = polar + hydro + posit + negat + speci
    # Safely find the cutoff index
    if short or end is None or str(end).lower() == "nan":
        cut = 0
    else:
        # np.argmin finds the index of the closest value to 'end', 
        # which is safer than exact '==' float matching.
        cut = np.argmin(np.abs(length[0] - end))
    x_dist = length[0][cut:]
    # ----------------------------
    # 1. Main Property Subplots
    # ----------------------------
    plt.figure(figsize=(15, 10))
    # Pack data into a list to loop through and avoid repetitive code
    plot_configs = [
        (polar, 'pink', 'Polar', 10, "Number of amino acids"),
        (hydro, 'orange', 'Hydrophobic', 10, "Number of amino acids"),
        (posit, 'red', 'Positively charged', 10, "Number of amino acids"),
        (negat, 'blue', 'Negatively charged', 10, "Number of amino acids"),
        (speci, 'silver', 'Special', 10, "Number of amino acids"),
        (rna, 'black', 'RNA', 27, "Number of nucleotides")
    ]
    for i, (data, color, title, y_max, y_label) in enumerate(plot_configs, 1):
        plt.subplot(2, 3, i)   
        mean_data = np.mean(data, axis=0)[cut:]
        std_data = np.std(data, axis=0)[cut:]
        plt.plot(x_dist, mean_data, color=color, lw=3)
        plt.fill_between(x_dist, mean_data - std_data, mean_data + std_data, color=color, alpha=0.5)
        plt.ylim(0, y_max)
        plt.xlim(0, 100)
        plt.xlabel("Distance from the PTC [Å]",fontsize=18)
        plt.ylabel(y_label,fontsize=18)
        plt.title(title,fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{output_name}_length.svg")
    plt.close()
    # ----------------------------
    # 2. RNA vs Protein Ratio Plot
    # ----------------------------
    plt.figure()
    mean_prot = np.mean(protein, axis=0)
    mean_rna = np.mean(rna, axis=0)
    std_prot = np.std(protein, axis=0)
    std_rna = np.std(rna, axis=0)
    total = mean_prot + mean_rna
    # np.errstate prevents terminal warnings if 'total' is temporarily zero
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = mean_rna / total
        # Propagation of uncertainty formula
        ratio_std = np.sqrt((mean_prot / total**2)**2 * std_rna**2 + (mean_rna / total**2)**2 * std_prot**2)
        # Replace any NaNs generated by division-by-zero with 0
        ratio = np.nan_to_num(ratio)
        ratio_std = np.nan_to_num(ratio_std)
    plt.plot(x_dist, ratio[cut:], color="black", lw=3)
    plt.fill_between(x_dist, ratio[cut:] - ratio_std[cut:], ratio[cut:] + ratio_std[cut:], color="black", alpha=0.2)
    plt.axhline(y=0.5, ls="--", color="gray")
    plt.xlabel("Distance from the PTC [Å]",fontsize=18)
    plt.ylabel("RNA vs Protein ratio",fontsize=18)
    plt.ylim(0, 1.1)
    plt.xlim(0, 100)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(f"{output_name}_ratio.svg")
    plt.close()
    # Save the math to files safely
    np.savetxt(f"{output_name}_ratio.dat", ratio)
    np.savetxt(f"{output_name}_ratio_std.dat", ratio_std)



###########################################
# Functions to plot tunnel cross sections #
###########################################

def plot_avg_cross_section(org_dict, org_ids, output_name, end=None, short=False):
    """
    Plots the average cross-sectional area of the tunnel as a function 
    of distance from the PTC, across different nascent chain lengths.
    """
    fig = plt.figure(figsize=(15, 10))
    colors = ['C0', 'C1', 'C2', 'C3', 'C4']    
    if not short:
        tunnel_length = org_dict[org_ids[0]]['path length']['60'][:-1]
        chain_lengths = ["10", "20", "30", "40", "60"]
    else:
        tunnel_length = org_dict[org_ids[0]]['path length']['40'][:-1]
        chain_lengths = ["10", "20", "30", "40"]
    # Safely find the cutoff index
    if short or end is None or str(end).lower() == "nan":
        cut = 0
    else:
        cut = np.argmin(np.abs(tunnel_length - end))
    x_vals = tunnel_length[cut:]
    # Dictionary to cache means and stds so we don't have to recalculate them for subplot 6
    computed_data = {}
    # ----------------------------------------------------
    # Plot individual subplots (1 to N)
    # ----------------------------------------------------
    for j, item in enumerate(chain_lengths):
        # Extract data into a 2D float array safely
        path_matrix = np.array([org_dict[org]['cross section'][item] for org in org_ids], dtype=float)
        # Calculate statistics, ignoring NaNs
        mean_vals = np.nanmean(path_matrix, axis=0)[cut:]
        std_vals = np.nanstd(path_matrix, axis=0)[cut:]
        # Cache for later
        computed_data[item] = (mean_vals, std_vals)
        plt.subplot(2, 3, j + 1)
        color = colors[j]
        plt.plot(x_vals, mean_vals, lw=4, color=color)
        plt.fill_between(x_vals, mean_vals - std_vals, mean_vals + std_vals, color=color, alpha=0.5)
        plt.title(f"N = {item} amino acids",fontsize=16)
        plt.xlim(0, 100)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # Custom Y-limits based on chain length
        if not short and item == "60":
            plt.axvline(x=92, ls='--', color='grey', alpha=0.7)
            plt.ylim(0, 2500)
        else:
            if not short and item == "40":
                plt.axvline(x=92, ls='--', color='grey', alpha=0.7)
            plt.ylim(0, 1500)
    # ----------------------------------------------------
    # Plot the final combined subplot (Subplot 6)
    # ----------------------------------------------------
    plt.subplot(2, 3, 6)
    for j, item in enumerate(chain_lengths):
        mean_vals, std_vals = computed_data[item]
        color = colors[j]
        # Add label for the legend
        plt.plot(x_vals, mean_vals, lw=4, color=color, label=f"N={item}")
        # Dropped alpha to 0.2 here so overlapping standard deviations don't completely block each other out
        plt.fill_between(x_vals, mean_vals - std_vals, mean_vals + std_vals, color=color, alpha=0.2)
    plt.title("All lengths",fontsize=16)
    if not short:
        plt.axvline(x=92, ls='--', color='grey', alpha=0.7)
    plt.xlim(0, 100)        
    plt.ylim(0, 1500)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc="upper right") # Display the legend
    # ----------------------------------------------------
    # Global Formatting
    # ----------------------------------------------------
    # Using Å instead of A for Angstroms
    fig.text(0.5, 0.02, 'Distance from the PTC [Å]', ha='center', va='center', fontsize=18)
    fig.text(0.02, 0.5, 'Cross section surface [Å²]', ha='center', va='center', rotation='vertical', fontsize=18)
    # rect leaves space for the global text labels we just added
    plt.tight_layout(rect=[0.03, 0.03, 1, 1]) 
    plt.savefig(output_name)
    plt.close()


# For nice figure separation
def find_closest_factors(X):
    """Calculates the optimal grid layout (rows, columns) for subplots."""
    cols = int(np.ceil(np.sqrt(X)))
    rows = int(np.ceil(X / cols))
    return rows, cols

def plot_cross_section(org_dict, org_ids, output_name, end=None, short=False):
    """
    Plots the individual smoothed cross-sectional areas for every organism
    in a grid layout.
    """
    fig = plt.figure(figsize=(18, 12)) # Slightly larger to accommodate the grid
    colors = ['C0', 'C1', 'C2', 'C3', 'C4']
    chain_lengths = ["10", "20", "30", "40"] if short else ["10", "20", "30", "40", "60"]
    # Get base tunnel length array
    tunnel_length = org_dict[org_ids[0]]['path length']['60' if not short else '40'][:-1]
    # Safely find the cutoff index
    if short or end is None or str(end).lower() == "nan":
        cut = 0
    else:
        cut = np.argmin(np.abs(tunnel_length - end))
    x_vals = tunnel_length[cut:]
    rows, cols = find_closest_factors(len(org_ids))
    for i, name in enumerate(org_ids, 1):
        plt.subplot(rows, cols, i)
        for j, item in enumerate(chain_lengths):
            y_raw = np.array(org_dict[name]['cross section'][item])
            # --- SAFE GAUSSIAN SMOOTHING ---
            valid_mask = ~np.isnan(y_raw)
            if np.any(valid_mask): # Only process if there is actual data!
                # 1. Temporarily interpolate over NaNs
                indices = np.arange(len(y_raw))
                y_interp = np.interp(indices, indices[valid_mask], y_raw[valid_mask])
                # 2. Apply smoothing
                y_smooth = gaussian_filter1d(y_interp, sigma=1)
                # 3. Put the NaNs back to prevent plotting fake physical data
                y_smooth[~valid_mask] = np.nan
            else:
                y_smooth = y_raw # If entirely empty, just keep the NaNs
            # Only add labels for the very first subplot to keep the legend clean
            label = f"NC length {item}" if i == 1 else None
            plt.plot(x_vals, y_smooth[cut:], color=colors[j], label=label, lw=3)
        plt.title(org_dict[name]['name'], style='italic', fontsize=16)
        plt.xlim(0, 100)
        plt.ylim(0, 1500)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
    # Adjust spacing between subplots
    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    # Add the global legend securely at the top
    fig.legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 0.98), fontsize=16)
    # Global Axis Labels
    fig.text(0.5, 0.02, 'Distance from the PTC [Å]', ha='center', va='center', fontsize=18)
    fig.text(0.02, 0.5, 'Cross section surface [Å²]', ha='center', va='center', rotation='vertical', fontsize=18)
    # Make room for the legend and labels so nothing overlaps
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.93])
    plt.savefig(output_name)
    plt.close()


###########################################################
# Plotting cross-sections grouped by nascent chain length #
###########################################################

def plot_cross_section_length(org_dict, org_ids, output_name, end=None, short=False, highlight_orgs=None):
    """
    Plots the smoothed cross-sectional areas, grouping all organisms onto a single plot 
    per chain length. Highlights specific divergent species in black.
    """
    if highlight_orgs is None:
        # Default highlights (e.g., specific microsporidia)
        highlight_orgs = ["ecun", "ploc", "vnec", "slop"]
    fig = plt.figure(figsize=(15, 10))
    colors = ['C0', 'C1', 'C2', 'C3', 'C4']
    chain_lengths = ["10", "20", "30", "40"] if short else ["10", "20", "30", "40", "60"]
    tunnel_length = org_dict[org_ids[0]]['path length']['40' if short else '60'][:-1]
    # Safely find the cutoff index
    if short or end is None or str(end).lower() == "nan":
        cut = 0
    else:
        cut = np.argmin(np.abs(tunnel_length - end))
    x_vals = tunnel_length[cut:]
    # Dictionary to cache the smoothed arrays so we don't recalculate them for Subplot 6
    smoothed_cache = {item: {} for item in chain_lengths}
    # ----------------------------------------------------
    # Plot individual subplots for each chain length
    # ----------------------------------------------------
    for j, item in enumerate(chain_lengths):
        plt.subplot(2, 3, j + 1)
        for name in org_ids:
            y_raw = np.array(org_dict[name]['cross section'][item])
            # --- SAFE GAUSSIAN SMOOTHING ---
            valid_mask = ~np.isnan(y_raw)
            if np.any(valid_mask):
                indices = np.arange(len(y_raw))
                y_interp = np.interp(indices, indices[valid_mask], y_raw[valid_mask])
                y_smooth = gaussian_filter1d(y_interp, sigma=1)
                y_smooth[~valid_mask] = np.nan
            else:
                y_smooth = y_raw
            # Cache the sliced data for later
            y_plot = y_smooth[cut:]
            smoothed_cache[item][name] = y_plot
            # Styling: Make highlighted organisms black, fully opaque, and draw them on top
            is_highlight = name in highlight_orgs
            color = "black" if is_highlight else colors[j]
            alpha = 1.0 if is_highlight else 0.8 # Fade the background organisms slightly
            zorder = 10 if is_highlight else 1   # Force black lines to the front
            plt.plot(x_vals, y_plot, color=color, lw=3, alpha=alpha, zorder=zorder)
        plt.xlim(0, 100)
        # Consistent Y-limits and PTC markers
        if not short and item == "60":
            plt.axvline(x=92, ls='--', color='grey', alpha=0.7) 
            plt.ylim(0, 2500)
        elif not short and item == "40":
            plt.axvline(x=92, ls='--', color='grey', alpha=0.7)
            plt.ylim(0, 1500)
        else: 
            plt.ylim(0, 1500)
        plt.title(f"N = {item} amino acids", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
    # ----------------------------------------------------
    # Plot Subplot 6 (All lengths, All organisms)
    # ----------------------------------------------------
    plt.subplot(2, 3, 6)
    for j, item in enumerate(chain_lengths):
        for name in org_ids:
            # Instantly fetch the pre-smoothed math
            y_plot = smoothed_cache[item][name]
            is_highlight = name in highlight_orgs
            color = "black" if is_highlight else colors[j]
            # Drop alpha heavily here so the "spaghetti" of 85 lines doesn't block out the shape
            alpha = 0.8 if is_highlight else 0.8 
            zorder = 10 if is_highlight else 1
            plt.plot(x_vals, y_plot, color=color, lw=3, alpha=alpha, zorder=zorder)
    if not short:
        plt.axvline(x=92, ls='--', color='grey', alpha=0.7) 
    plt.xlim(0, 100)
    plt.ylim(0, 1500)
    plt.title("All lengths", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # ----------------------------------------------------
    # Global Formatting
    # ----------------------------------------------------
    fig.text(0.5, 0.02, 'Distance from the PTC [Å]', ha='center', va='center', fontsize=18)
    fig.text(0.02, 0.5, 'Cross section surface [Å²]', ha='center', va='center', rotation='vertical', fontsize=18)
    plt.tight_layout(rect=[0.03, 0.03, 1, 1])
    plt.savefig(output_name)
    plt.close()
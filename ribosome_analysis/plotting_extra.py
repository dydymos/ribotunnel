#######################
# Plotting paths asph #
#######################
import numpy as np
import matplotlib.pyplot as plt

# Load data
end = 91.71631884486219
length = np.loadtxt("length.dat")[:-1]
cut = np.where(length == end)[0][0]

b_path = np.loadtxt("bacteria_avg_asph.dat")
b_std = np.loadtxt("bacteria_std_asph.dat")
a_path = np.loadtxt("archaea_avg_asph.dat")
a_std = np.loadtxt("archaea_std_asph.dat")
e_path = np.loadtxt("eukaryota_core_avg_asph.dat")
e_std = np.loadtxt("eukaryota_core_std_asph.dat")

# Softer, distinct colors
colors = {
    "Bacteria": "#D36B8B",    # light burgundy rose
    "Eukaryota": "#66B2B2",   # light teal
    "Archaea": "#D4AF37",     # light gold
}

titles ={
    0:"L=10 amino acids",
    1:"L=20 amino acids",
    2:"L=30 amino acids",
    3:"L=40 amino acids",
    4:"L=60 amino acids",
}


plt.figure(figsize=(12, 6))
for i in range(0,5):
    plt.subplot(2, 3, i + 1)
    # Bacteria
    plt.plot(length[cut:], b_path[i][cut:], label="Bacteria", color=colors["Bacteria"], linewidth=4)
    plt.fill_between(length[cut:], 
                     b_path[i][cut:] - b_std[i][cut:], 
                     b_path[i][cut:] + b_std[i][cut:], 
                     color=colors["Bacteria"], alpha=0.3)
    # Eukaryota
    plt.plot(length[cut:], e_path[i][cut:], label="Eukaryota", color=colors["Eukaryota"], linewidth=4)
    plt.fill_between(length[cut:], 
                     e_path[i][cut:] - e_std[i][cut:], 
                     e_path[i][cut:] + e_std[i][cut:], 
                     color=colors["Eukaryota"], alpha=0.3)
    # Archaea
    plt.plot(length[cut:], a_path[i][cut:], label="Archaea", color=colors["Archaea"], linewidth=4)
    plt.fill_between(length[cut:], 
                     a_path[i][cut:] - a_std[i][cut:], 
                     a_path[i][cut:] + a_std[i][cut:], 
                     color=colors["Archaea"], alpha=0.3)                   
    # Labels and legend
    plt.xlabel("Distance from the PTC [Å]")
    plt.ylabel("Tunnel Asphericity")
    plt.title(titles[i])
    plt.ylim(0,1)
    if (i==0): plt.legend()

plt.tight_layout()

plt.show()

######################
# Plotting avg paths #
######################
import numpy as np
import matplotlib.pyplot as plt

# Load data
end = 91.71631884486219
length = np.loadtxt("length.dat")[:-1]
cut = np.where(length == end)[0][0]
#length_mito = np.loadtxt("mitochondria/lmaj/tunnel_mito.dat")
#length_mito_60 = np.loadtxt("mitochondria/hsap/tunnel_mito_60.dat")

b_path = np.loadtxt("bacteria_avg_path.dat")
b_std = np.loadtxt("bacteria_std_path.dat")
a_path = np.loadtxt("archaea_avg_path.dat")
a_std = np.loadtxt("archaea_std_path.dat")
e_path = np.loadtxt("eukaryota_core_avg_path.dat")
e_std = np.loadtxt("eukaryota_core_std_path.dat")
#mito_path = np.loadtxt("mitochondria/lmaj/mitochondria_avg_path.dat")
#mito_path_60 = np.loadtxt("mitochondria/hsap/mitochondria_avg_path_60.dat")
#mito_std = np.loadtxt("mitochondria/lmaj/mitochondria_std_path.dat")
#mito_std_60 = np.loadtxt("mitochondria/hsap/mitochondria_std_path_60.dat")

b_ratio = np.loadtxt("bacteria_tunnel_composition_ratio.dat")
a_ratio = np.loadtxt("archaea_tunnel_composition_ratio.dat")
e_ratio = np.loadtxt("eukaryota_core_tunnel_composition_ratio.dat")
b_ratio_std = np.loadtxt("bacteria_tunnel_composition_ratio_std.dat")
a_ratio_std = np.loadtxt("archaea_tunnel_composition_ratio_std.dat")
e_ratio_std = np.loadtxt("eukaryota_core_tunnel_composition_ratio_std.dat")
#mito_ratio = np.loadtxt("mitochondria/lmaj/mitochondria_tunnel_composition_ratio.dat")

# Softer, distinct colors
colors = {
    "Bacteria": "#D36B8B",    # light burgundy rose
    "Eukaryota": "#66B2B2",   # light teal
    "Archaea": "#D4AF37",     # light gold
    "Mitochondria": "black"
}

titles ={
    0:"L=10 amino acids",
    1:"L=20 amino acids",
    2:"L=30 amino acids",
    3:"L=40 amino acids",
    4:"L=60 amino acids",
}


plt.figure(figsize=(12, 6))
for i in range(0,5):
    plt.subplot(2, 3, i + 1)
    # Bacteria
    plt.plot(length[cut:], b_path[i][cut:], label="Bacteria", color=colors["Bacteria"], linewidth=4)
    plt.fill_between(length[cut:], 
                     b_path[i][cut:] - b_std[i][cut:], 
                     b_path[i][cut:] + b_std[i][cut:], 
                     color=colors["Bacteria"], alpha=0.3)
    # Eukaryota
    plt.plot(length[cut:], e_path[i][cut:], label="Eukaryota", color=colors["Eukaryota"], linewidth=4)
    plt.fill_between(length[cut:], 
                     e_path[i][cut:] - e_std[i][cut:], 
                     e_path[i][cut:] + e_std[i][cut:], 
                     color=colors["Eukaryota"], alpha=0.3)
    # Archaea
    plt.plot(length[cut:], a_path[i][cut:], label="Archaea", color=colors["Archaea"], linewidth=4)
    plt.fill_between(length[cut:], 
                     a_path[i][cut:] - a_std[i][cut:], 
                     a_path[i][cut:] + a_std[i][cut:], 
                     color=colors["Archaea"], alpha=0.3)
    # Mitochondria
    #plt.plot(length_mito, mito_path[i], label="Mitochondria", color="grey", linewidth=4)
    
    # Labels and legend
    plt.xlabel("Distance from the PTC [Å]")
    plt.ylabel("Cross-section [Å$^2$]")
    plt.title(titles[i])
    if i==0: 
        plt.legend()
        plt.ylim(0,500)
    if i==1 or i==2:
        plt.ylim(0,500)
    if i==3 or i==4:
        plt.ylim(0,1200)

plt.tight_layout()

plt.subplot(2, 3, 6)
plt.axhline(0.5,color="black",ls='--',lw=2)
plt.plot(length[cut:][:-1], b_ratio[cut:][:-1], label="Bacteria", color=colors["Bacteria"], linewidth=4)
plt.fill_between(length[cut:][:-1],b_ratio[cut:][:-1] - b_ratio_std[cut:][:-1],b_ratio[cut:][:-1] + b_ratio_std[cut:][:-1],color=colors["Bacteria"], alpha=0.3)
plt.plot(length[cut:][:-1], e_ratio[cut:][:-1], label="Eukaryota", color=colors["Eukaryota"], linewidth=4)
plt.fill_between(length[cut:][:-1],e_ratio[cut:][:-1] - e_ratio_std[cut:][:-1],e_ratio[cut:][:-1] + e_ratio_std[cut:][:-1],color=colors["Eukaryota"], alpha=0.3)
plt.plot(length[cut:][:-1], a_ratio[cut:][:-1], label="Archaea", color=colors["Archaea"], linewidth=4)
plt.fill_between(length[cut:][:-1],a_ratio[cut:][:-1] - a_ratio_std[cut:][:-1],a_ratio[cut:][:-1] + a_ratio_std[cut:][:-1],color=colors["Archaea"], alpha=0.3)
#plt.plot(length_mito, mito_ratio, label="Mitochondria", color="grey", linewidth=4)

plt.xlabel("Distance from the PTC [Å]")
plt.ylabel("RNA vs Protein ratio")
plt.ylim(0,1.1)

plt.show()


######################
# Plotting avg paths #
# with microsporidia #
######################
import numpy as np
import matplotlib.pyplot as plt

# Load data
end = 91.71631884486219
length = np.loadtxt("length.dat")[:-1]
cut = np.where(length == end)[0][0]
#length_mito = np.loadtxt("mitochondria/lmaj/tunnel_mito.dat")
#length_mito_60 = np.loadtxt("mitochondria/hsap/tunnel_mito_60.dat")

b_path = np.loadtxt("bacteria_avg_path.dat")
b_std = np.loadtxt("bacteria_std_path.dat")
a_path = np.loadtxt("archaea_avg_path.dat")
a_std = np.loadtxt("archaea_std_path.dat")
e_path = np.loadtxt("eukaryota_core_avg_path.dat")
e_std = np.loadtxt("eukaryota_core_std_path.dat")
micro_path = np.loadtxt("micro_avg_path.dat")
micro_std = np.loadtxt("micro_std_path.dat")
#mito_path = np.loadtxt("mitochondria/lmaj/mitochondria_avg_path.dat")
#mito_path_60 = np.loadtxt("mitochondria/hsap/mitochondria_avg_path_60.dat")
#mito_std = np.loadtxt("mitochondria/lmaj/mitochondria_std_path.dat")
#mito_std_60 = np.loadtxt("mitochondria/hsap/mitochondria_std_path_60.dat")

b_ratio = np.loadtxt("bacteria_tunnel_composition_ratio.dat")
a_ratio = np.loadtxt("archaea_tunnel_composition_ratio.dat")
e_ratio = np.loadtxt("eukaryota_core_tunnel_composition_ratio.dat")
micro_ratio = np.loadtxt("micro_tunnel_composition_ratio.dat")
b_ratio_std = np.loadtxt("bacteria_tunnel_composition_ratio_std.dat")
a_ratio_std = np.loadtxt("archaea_tunnel_composition_ratio_std.dat")
e_ratio_std = np.loadtxt("eukaryota_core_tunnel_composition_ratio_std.dat")
micro_ratio_std = np.loadtxt("micro_tunnel_composition_ratio_std.dat")
#mito_ratio = np.loadtxt("mitochondria/lmaj/mitochondria_tunnel_composition_ratio.dat")

# Softer, distinct colors
colors = {
    "Bacteria": "#D36B8B",    # light burgundy rose
    "Eukaryota": "#66B2B2",   # light teal
    "Archaea": "#D4AF37",     # light gold
    "Microsporidia": "black"
}

titles ={
    0:"L=10 amino acids",
    1:"L=20 amino acids",
    2:"L=30 amino acids",
    3:"L=40 amino acids",
    4:"L=60 amino acids",
}


plt.figure(figsize=(12, 6))
for i in range(0,5):
    plt.subplot(2, 3, i + 1)
    # Bacteria
    plt.plot(length[cut:], b_path[i][cut:], label="Bacteria", color=colors["Bacteria"], linewidth=4)
    plt.fill_between(length[cut:], 
                     b_path[i][cut:] - b_std[i][cut:], 
                     b_path[i][cut:] + b_std[i][cut:], 
                     color=colors["Bacteria"], alpha=0.3)
    # Eukaryota
    plt.plot(length[cut:], e_path[i][cut:], label="Eukaryota", color=colors["Eukaryota"], linewidth=4)
    plt.fill_between(length[cut:], 
                     e_path[i][cut:] - e_std[i][cut:], 
                     e_path[i][cut:] + e_std[i][cut:], 
                     color=colors["Eukaryota"], alpha=0.3)
    # Archaea
    plt.plot(length[cut:], a_path[i][cut:], label="Archaea", color=colors["Archaea"], linewidth=4)
    plt.fill_between(length[cut:], 
                     a_path[i][cut:] - a_std[i][cut:], 
                     a_path[i][cut:] + a_std[i][cut:], 
                     color=colors["Archaea"], alpha=0.3)
    # Microsporidia
    plt.plot(length[cut:], micro_path[i][cut:], label="Microsporidia", color=colors["Microsporidia"], linewidth=4)
    plt.fill_between(length[cut:], 
                     micro_path[i][cut:] - micro_std[i][cut:], 
                     micro_path[i][cut:] + micro_std[i][cut:], 
                     color=colors["Microsporidia"], alpha=0.3)
    
    # Labels and legend
    plt.xlabel("Distance from the PTC [Å]")
    plt.ylabel("Cross-section [Å$^2$]")
    plt.title(titles[i])
    if (i==0): plt.legend()

plt.tight_layout()

plt.subplot(2, 3, 6)
plt.axhline(0.5,color="black",ls='--',lw=2)
plt.plot(length[cut:], b_ratio[cut:], label="Bacteria", color=colors["Bacteria"], linewidth=4)
plt.fill_between(length[cut:],b_ratio[cut:] - b_ratio_std[cut:],b_ratio[cut:] + b_ratio_std[cut:],color=colors["Bacteria"], alpha=0.3)
plt.plot(length[cut:], e_ratio[cut:], label="Eukaryota", color=colors["Eukaryota"], linewidth=4)
plt.fill_between(length[cut:],e_ratio[cut:] - e_ratio_std[cut:],e_ratio[cut:] + e_ratio_std[cut:],color=colors["Eukaryota"], alpha=0.3)
plt.plot(length[cut:], a_ratio[cut:], label="Archaea", color=colors["Archaea"], linewidth=4)
plt.fill_between(length[cut:],a_ratio[cut:] - a_ratio_std[cut:],a_ratio[cut:] + a_ratio_std[cut:],color=colors["Archaea"], alpha=0.3)
plt.plot(length[cut:], micro_ratio[cut:], label="Microsporidia", color=colors["Microsporidia"], linewidth=4)
plt.fill_between(length[cut:],micro_ratio[cut:] - micro_ratio_std[cut:],micro_ratio[cut:] + micro_ratio_std[cut:],color=colors["Microsporidia"], alpha=0.3)
#plt.plot(length_mito, mito_ratio, label="Mitochondria", color="grey", linewidth=4)

plt.xlabel("Distance from the PTC [Å]")
plt.ylabel("RNA vs Protein ratio")
plt.ylim(0,1.1)

plt.show()

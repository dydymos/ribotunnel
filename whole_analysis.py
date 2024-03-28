import numpy as np
import matplotlib.pyplot as plt
import mrcfile
import MDAnalysis as mda
import tools
from tools import *
import os


# Getting maps
keys = ['60','40','30','20','10']
maps = dict()
for i in keys:
    maps[i]=dict()
    maps[i]['0'] = get_map_param(i+"/average_final.mrc")
    maps[i]['1'] = get_map_param(i+"/ref_1/average_final.mrc")
    maps[i]['2'] = get_map_param(i+"/ref_2/average_final.mrc")


epsilon = 0.001

# The main path
# Generating the main path PDB file
for i in keys:
    get_main_path(maps[i]['0'],epsilon,i+"/main_path.pdb")
    get_main_path(maps[i]['1'],epsilon,i+"/ref_1/main_path.pdb")
    get_main_path(maps[i]['2'],epsilon,i+"/ref_2/main_path.pdb")

# Reading the main path file
main_path = dict()
for i in keys:
    main_path[i] = dict()
    main_path[i]['0'] = mda.Universe(i+'/main_path.pdb')
    main_path[i]['1'] = mda.Universe(i+'/ref_1/main_path.pdb')
    main_path[i]['2'] = mda.Universe(i+'/ref_2/main_path.pdb')

# Getting geometry
geom = dict()
for i in keys:
    geom[i] = dict()
    geom[i]['0'] = geometric_prop(maps[i]['0'],epsilon)
    geom[i]['1'] = geometric_prop(maps[i]['1'],epsilon)
    geom[i]['2'] = geometric_prop(maps[i]['2'],epsilon)


# number of atom in the main path
n = dict()
for i in keys:
    n_list = []
    for ii in main_path[i]:
        n_list.append(len(geom[i][ii]['z_diam']))
    n[i] = np.min(n_list)


# Averages and std errors
geom_avg = dict()
for i in keys:
    geom_avg[i] = dict()
    geom_avg[i]['z_diam'] = np.mean([geom[i]['0']['z_diam'][:n[i]],geom[i]['1']['z_diam'][:n[i]],geom[i]['2']['z_diam'][:n[i]]],0)
    geom_avg[i]['x_diam'] = np.mean([geom[i]['0']['x_diam'][:n[i]],geom[i]['1']['x_diam'][:n[i]],geom[i]['2']['x_diam'][:n[i]]],0)
    geom_avg[i]['vol'] = np.mean([geom[i]['0']['vol'][:n[i]],geom[i]['1']['vol'][:n[i]],geom[i]['2']['vol'][:n[i]]],0)
    geom_avg[i]['asp_rat'] = np.mean([geom[i]['0']['asp_rat'][:n[i]],geom[i]['1']['asp_rat'][:n[i]],geom[i]['2']['asp_rat'][:n[i]]],0)


geom_std = dict()
for i in keys:
    geom_std[i] = dict()
    geom_std[i]['z_diam'] = np.std([geom[i]['0']['z_diam'][:n[i]],geom[i]['1']['z_diam'][:n[i]],geom[i]['2']['z_diam'][:n[i]]],0)
    geom_std[i]['x_diam'] = np.std([geom[i]['0']['x_diam'][:n[i]],geom[i]['1']['x_diam'][:n[i]],geom[i]['2']['x_diam'][:n[i]]],0)
    geom_std[i]['vol'] = np.std([geom[i]['0']['vol'][:n[i]],geom[i]['1']['vol'][:n[i]],geom[i]['2']['vol'][:n[i]]],0)
    geom_std[i]['asp_rat'] = np.std([geom[i]['0']['asp_rat'][:n[i]],geom[i]['1']['asp_rat'][:n[i]],geom[i]['2']['asp_rat'][:n[i]]],0)


# RMSD
rmsd = dict()
for i in keys:
    rmsd[i] = []
    rmsd[i].append(np.sqrt(np.sum((main_path[i]['0'].atoms.positions[:n[i]]-main_path[i]['1'].atoms.positions[:n[i]])**2,1)))
    rmsd[i].append(np.sqrt(np.sum((main_path[i]['0'].atoms.positions[:n[i]]-main_path[i]['2'].atoms.positions[:n[i]])**2,1)))
    rmsd[i].append(np.sqrt(np.sum((main_path[i]['1'].atoms.positions[:n[i]]-main_path[i]['2'].atoms.positions[:n[i]])**2,1)))
    rmsd[i] = np.array(rmsd[i])

rmsd_avg = dict()
rmsd_std = dict()
for i in keys:
    rmsd_avg[i] = np.mean(rmsd[i],0)
    rmsd_std[i] = np.std(rmsd[i],0)

############
# Plotting #
############
kolors = ['C0','C1','C2','C3','C4']

### RMSD
j=0
for i in keys:
    plt.errorbar(np.arange(0,len(geom[i]['0']['x_diam'][:n[i]]))*maps[i]['0']['vox'],rmsd_avg[i],yerr=rmsd_std[i],lw=2,label=i+" aa")
    plt.plot(np.arange(0,len(geom[i]['0']['x_diam'][:n[i]]))*maps[i]['0']['vox'],rmsd[i][0][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['1']['x_diam'][:n[i]]))*maps[i]['1']['vox'],rmsd[i][1][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['2']['x_diam'][:n[i]]))*maps[i]['2']['vox'],rmsd[i][2][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    j+=1

plt.xlabel("Distance from the PTC along the Y-axis [A]")
plt.ylabel("RMSD between central paths [A]")
plt.title('RMSD')
plt.legend()
plt.xlim(0,120)
plt.ylim(0,15)
plt.title(title)
# Save the figure as an SVG file
plt.savefig('rmsd_plot.svg', format='svg')
plt.close()

# Volume
j=0
for i in keys:
    plt.errorbar(np.arange(0,len(geom[i]['0']['z_diam'][:n[i]]))*maps[i]['0']['vox'],geom_avg[i]['vol'],yerr = geom_std[i]['vol'],lw=2,label=i+' aa')
    plt.plot(np.arange(0,len(geom[i]['0']['z_diam'][:n[i]]))*maps[i]['0']['vox'],geom[i]['0']['vol'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['1']['z_diam'][:n[i]]))*maps[i]['1']['vox'],geom[i]['1']['vol'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['2']['z_diam'][:n[i]]))*maps[i]['2']['vox'],geom[i]['2']['vol'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    j+=1

plt.xlabel("Distance from the PTC along the Y-axis [A]")
plt.ylabel("Volume of the section of the exit tunne [$A^3$]")
plt.legend()
plt.xlim(0,120)
plt.ylim(0,2500)
plt.title(title)
# Save the figure as an SVG file
plt.savefig('volume_plot.svg', format='svg')
plt.close()

### Distance from Axis
plt.figure(figsize=(16,4))
plt.subplot(1,2,1)
j=0
for i in keys:
    plt.errorbar(np.arange(0,len(geom[i]['0']['x_diam'][:n[i]]))*maps[i]['0']['vox'],geom_avg[i]['x_diam'],yerr = geom_std[i]['x_diam'],lw=2,label=i+' aa')
    plt.plot(np.arange(0,len(geom[i]['0']['x_diam'][:n[i]]))*maps[i]['0']['vox'],geom[i]['0']['x_diam'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['1']['x_diam'][:n[i]]))*maps[i]['1']['vox'],geom[i]['1']['x_diam'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['2']['x_diam'][:n[i]]))*maps[i]['2']['vox'],geom[i]['2']['x_diam'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    j+=1

plt.legend()
plt.xlabel("Distance from the PTC along the Y-axis [A]")
plt.ylabel("Distance from the main path along the X-axis [A]")
plt.xlim(0,120)
plt.ylim(0,60)
plt.title(title+" "+"X-axis")

plt.subplot(1,2,2)
j=0
for i in keys:
    plt.errorbar(np.arange(0,len(geom[i]['0']['z_diam'][:n[i]]))*maps[i]['0']['vox'],geom_avg[i]['z_diam'],yerr = geom_std[i]['z_diam'],lw=2,label=i+' aa')
    plt.plot(np.arange(0,len(geom[i]['0']['z_diam'][:n[i]]))*maps[i]['0']['vox'],geom[i]['0']['z_diam'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['1']['z_diam'][:n[i]]))*maps[i]['1']['vox'],geom[i]['1']['z_diam'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['2']['z_diam'][:n[i]]))*maps[i]['2']['vox'],geom[i]['2']['z_diam'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    j+=1

plt.legend()
plt.xlabel("Distance from the PTC along the Y-axis [A]")
plt.ylabel("Distance from the main path along the Z-axis [A]")
plt.xlim(0,120)
plt.ylim(0,60)
plt.title(title+" "+"Z-axis")
plt.savefig('axis_plot.svg', format='svg')
plt.close()

# Aspect ratio plot
j=0
for i in keys:
    plt.errorbar(np.arange(0,len(geom[i]['0']['z_diam'][:n[i]]))*maps[i]['0']['vox'],geom_avg[i]['asp_rat'],yerr = geom_std[i]['asp_rat'],lw=2,label=i+' aa')
    plt.plot(np.arange(0,len(geom[i]['0']['z_diam'][:n[i]]))*maps[i]['0']['vox'],geom[i]['0']['asp_rat'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['1']['z_diam'][:n[i]]))*maps[i]['1']['vox'],geom[i]['1']['asp_rat'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    plt.plot(np.arange(0,len(geom[i]['2']['z_diam'][:n[i]]))*maps[i]['2']['vox'],geom[i]['2']['asp_rat'][:n[i]],lw=2,color=kolors[j],alpha=0.3)
    j+=1

plt.xlabel("Distance from the PTC [A]")
plt.ylabel("Aspect ratio X/Z dimension")
plt.axhline(y=1.0, color='black',ls='--')
plt.xlim(0,120)
plt.ylim(0.25,3.0)
plt.title(title)
plt.legend()
# Save the figure as an SVG file
plt.savefig('aspect_ratio_plot.svg', format='svg')
plt.close()

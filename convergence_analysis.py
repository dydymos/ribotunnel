import numpy as np
import matplotlib.pyplot as plt
import mrcfile
import MDAnalysis as mda
import tools
from tools import *
import os


path = "/media/didymos/Projects/ribosome_tunnels/methionine/"
os.chdir(path)

######################
# Bacteria ribosomes #
######################
bac_id = ['4ybb']

bac_name = ['E.coli']

bac_dict = dict()
bac_dict_long = dict()
i=0
for item in bac_id:
    bac_dict[item] = dict()
    bac_dict_long[item] = dict()
    bac_dict[item]['name'] = bac_name[i]
    bac_dict_long[item]['name'] = bac_name[i]
    i+=1


################
# Getting maps #
################

keys = ['60','40','30','20','10']

for item in bac_id:
    bac_dict[item]['maps'] = dict()
    bac_dict_long[item]['maps'] = dict()
    for i in keys:
        bac_dict[item]['maps'][i] = dict()
        bac_dict_long[item]['maps'][i] = dict()
        bac_dict[item]['maps'][i]['0'] = get_map_param(item+'/'+i+"/average_final.mrc")
        bac_dict[item]['maps'][i]['1'] = get_map_param(item+'/'+i+"/ref_1/average_final.mrc")
        bac_dict[item]['maps'][i]['2'] = get_map_param(item+'/'+i+"/ref_2/average_final.mrc")
        bac_dict_long[item]['maps'][i]['0'] = get_map_param(item+'/'+i+"/average_final_all.mrc")
        bac_dict_long[item]['maps'][i]['1'] = get_map_param(item+'/'+i+"/ref_1/average_final_all.mrc")
        bac_dict_long[item]['maps'][i]['2'] = get_map_param(item+'/'+i+"/ref_2/average_final_all.mrc")


#################
# The main path #
#################

epsilon = 0.001

# Generating the main path PDB file
for item in bac_id:
    for i in keys:
        get_main_path(bac_dict[item]['maps'][i]['0'],epsilon,item+'/'+i+"/main_path.pdb")
        get_main_path(bac_dict[item]['maps'][i]['1'],epsilon,item+'/'+i+"/ref_1/main_path.pdb")
        get_main_path(bac_dict[item]['maps'][i]['2'],epsilon,item+'/'+i+"/ref_2/main_path.pdb")
        get_main_path(bac_dict_long[item]['maps'][i]['0'],epsilon,item+'/'+i+"/main_path_all.pdb")
        get_main_path(bac_dict_long[item]['maps'][i]['1'],epsilon,item+'/'+i+"/ref_1/main_path_all.pdb")
        get_main_path(bac_dict_long[item]['maps'][i]['2'],epsilon,item+'/'+i+"/ref_2/main_path_all.pdb")

# Reading the main path file
for item in bac_id:
    bac_dict[item]['path'] = dict()
    bac_dict_long[item]['path'] = dict()
    for i in keys:
        bac_dict[item]['path'][i] = dict()
        bac_dict_long[item]['path'][i] = dict()
        bac_dict[item]['path'][i]['0'] = mda.Universe(item+'/'+i+'/main_path.pdb')
        bac_dict[item]['path'][i]['1'] = mda.Universe(item+'/'+i+'/ref_1/main_path.pdb')
        bac_dict[item]['path'][i]['2'] = mda.Universe(item+'/'+i+'/ref_2/main_path.pdb')
        bac_dict_long[item]['path'][i]['0'] = mda.Universe(item+'/'+i+'/main_path_all.pdb')
        bac_dict_long[item]['path'][i]['1'] = mda.Universe(item+'/'+i+'/ref_1/main_path_all.pdb')
        bac_dict_long[item]['path'][i]['2'] = mda.Universe(item+'/'+i+'/ref_2/main_path_all.pdb')


# Getting geometry
for item in bac_id:
    bac_dict[item]['geom'] = dict()
    bac_dict_long[item]['geom'] = dict()
    for i in keys:
        bac_dict[item]['geom'][i] = dict()
        bac_dict_long[item]['geom'][i] = dict()
        bac_dict[item]['geom'][i]['0'] = geometric_prop(bac_dict[item]['maps'][i]['0'],epsilon)
        bac_dict[item]['geom'][i]['1'] = geometric_prop(bac_dict[item]['maps'][i]['1'],epsilon)
        bac_dict[item]['geom'][i]['2'] = geometric_prop(bac_dict[item]['maps'][i]['2'],epsilon)
        bac_dict_long[item]['geom'][i]['0'] = geometric_prop(bac_dict_long[item]['maps'][i]['0'],epsilon)
        bac_dict_long[item]['geom'][i]['1'] = geometric_prop(bac_dict_long[item]['maps'][i]['1'],epsilon)
        bac_dict_long[item]['geom'][i]['2'] = geometric_prop(bac_dict_long[item]['maps'][i]['2'],epsilon)


# number of atom in the main path
for item in bac_id:
    bac_dict[item]['n'] = dict()
    bac_dict_long[item]['n'] = dict()
    for i in keys:
        n_list = []
        n_list_long = []
        for ii in bac_dict[item]['path'][i]:
            n_list.append(len(bac_dict[item]['geom'][i][ii]['z_diam']))
        for ii in bac_dict_long[item]['path'][i]:
            n_list_long.append(len(bac_dict_long[item]['geom'][i][ii]['z_diam']))
        bac_dict[item]['n'][i] = np.min(n_list)
        bac_dict_long[item]['n'][i] = np.min(n_list_long)


# Averages and std errors
for item in bac_id:
    bac_dict[item]['geom_avg'] = dict()
    bac_dict_long[item]['geom_avg'] = dict()
    for i in keys:
        bac_dict[item]['geom_avg'][i] = dict()
        bac_dict_long[item]['geom_avg'][i] = dict()
        bac_dict[item]['geom_avg'][i]['z_diam'] = np.mean([bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_avg'][i]['x_diam'] = np.mean([bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_avg'][i]['vol'] = np.mean([bac_dict[item]['geom'][i]['0']['vol'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['vol'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['vol'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_avg'][i]['asp_rat'] = np.mean([bac_dict[item]['geom'][i]['0']['asp_rat'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['asp_rat'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['asp_rat'][:bac_dict[item]['n'][i]]],0)
        bac_dict_long[item]['geom_avg'][i]['z_diam'] = np.mean([bac_dict_long[item]['geom'][i]['0']['z_diam'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['1']['z_diam'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['2']['z_diam'][:bac_dict_long[item]['n'][i]]],0)
        bac_dict_long[item]['geom_avg'][i]['x_diam'] = np.mean([bac_dict_long[item]['geom'][i]['0']['x_diam'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['1']['x_diam'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['2']['x_diam'][:bac_dict_long[item]['n'][i]]],0)
        bac_dict_long[item]['geom_avg'][i]['vol'] = np.mean([bac_dict_long[item]['geom'][i]['0']['vol'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['1']['vol'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['2']['vol'][:bac_dict_long[item]['n'][i]]],0)
        bac_dict_long[item]['geom_avg'][i]['asp_rat'] = np.mean([bac_dict_long[item]['geom'][i]['0']['asp_rat'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['1']['asp_rat'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['2']['asp_rat'][:bac_dict_long[item]['n'][i]]],0)



for item in bac_id:
    bac_dict[item]['geom_std'] = dict()
    bac_dict_long[item]['geom_std'] = dict()
    for i in keys:
        bac_dict[item]['geom_std'][i] = dict()
        bac_dict_long[item]['geom_std'][i] = dict()
        bac_dict[item]['geom_std'][i]['z_diam'] = np.std([bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_std'][i]['x_diam'] = np.std([bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_std'][i]['vol'] = np.std([bac_dict[item]['geom'][i]['0']['vol'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['vol'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['vol'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_std'][i]['asp_rat'] = np.std([bac_dict[item]['geom'][i]['0']['asp_rat'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['asp_rat'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['asp_rat'][:bac_dict[item]['n'][i]]],0)
        bac_dict_long[item]['geom_std'][i]['z_diam'] = np.std([bac_dict_long[item]['geom'][i]['0']['z_diam'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['1']['z_diam'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['2']['z_diam'][:bac_dict_long[item]['n'][i]]],0)
        bac_dict_long[item]['geom_std'][i]['x_diam'] = np.std([bac_dict_long[item]['geom'][i]['0']['x_diam'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['1']['x_diam'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['2']['x_diam'][:bac_dict_long[item]['n'][i]]],0)
        bac_dict_long[item]['geom_std'][i]['vol'] = np.std([bac_dict_long[item]['geom'][i]['0']['vol'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['1']['vol'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['2']['vol'][:bac_dict_long[item]['n'][i]]],0)
        bac_dict_long[item]['geom_std'][i]['asp_rat'] = np.std([bac_dict_long[item]['geom'][i]['0']['asp_rat'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['1']['asp_rat'][:bac_dict_long[item]['n'][i]],bac_dict_long[item]['geom'][i]['2']['asp_rat'][:bac_dict_long[item]['n'][i]]],0)



# RMSD
for item in bac_id:
    bac_dict[item]['rmsd_avg'] = dict()
    bac_dict[item]['rmsd_std'] = dict()
    bac_dict[item]['rmsd'] = dict()
    bac_dict_long[item]['rmsd_avg'] = dict()
    bac_dict_long[item]['rmsd_std'] = dict()
    bac_dict_long[item]['rmsd'] = dict()
    for i in keys:
        bac_dict[item]['rmsd'][i] = []
        bac_dict_long[item]['rmsd'][i] = []
        bac_dict[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict[item]['path'][i]['0'].atoms.positions[:bac_dict[item]['n'][i]]-bac_dict[item]['path'][i]['1'].atoms.positions[:bac_dict[item]['n'][i]])**2,1)))
        bac_dict[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict[item]['path'][i]['0'].atoms.positions[:bac_dict[item]['n'][i]]-bac_dict[item]['path'][i]['2'].atoms.positions[:bac_dict[item]['n'][i]])**2,1)))
        bac_dict[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict[item]['path'][i]['1'].atoms.positions[:bac_dict[item]['n'][i]]-bac_dict[item]['path'][i]['2'].atoms.positions[:bac_dict[item]['n'][i]])**2,1)))
        bac_dict[item]['rmsd_avg'][i] = np.mean(np.array(bac_dict[item]['rmsd'][i]),0)
        bac_dict[item]['rmsd_std'][i] = np.std(np.array(bac_dict[item]['rmsd'][i]),0)
        bac_dict_long[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict_long[item]['path'][i]['0'].atoms.positions[:bac_dict_long[item]['n'][i]]-bac_dict_long[item]['path'][i]['1'].atoms.positions[:bac_dict_long[item]['n'][i]])**2,1)))
        bac_dict_long[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict_long[item]['path'][i]['0'].atoms.positions[:bac_dict_long[item]['n'][i]]-bac_dict_long[item]['path'][i]['2'].atoms.positions[:bac_dict_long[item]['n'][i]])**2,1)))
        bac_dict_long[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict_long[item]['path'][i]['1'].atoms.positions[:bac_dict_long[item]['n'][i]]-bac_dict_long[item]['path'][i]['2'].atoms.positions[:bac_dict_long[item]['n'][i]])**2,1)))
        bac_dict_long[item]['rmsd_avg'][i] = np.mean(np.array(bac_dict_long[item]['rmsd'][i]),0)
        bac_dict_long[item]['rmsd_std'][i] = np.std(np.array(bac_dict_long[item]['rmsd'][i]),0)


############
# Plotting #
############
kolors = ['C0','C1','C2','C3','C4']

##################
### Pathway RMSD #
##################
plt.figure(figsize=(12,12))
plt.tight_layout()
for item in bac_id:
    plt.subplot(1,2,1)
    j=0
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['rmsd_avg'][i],yerr=bac_dict[item]['rmsd_std'][i],lw=2,label=i+" aa")
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['rmsd'][i][0][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['rmsd'][i][1][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['rmsd'][i][2][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.xlabel("Distance from the PTC along the Y-axis [A]")
    plt.ylabel("RMSD between central paths [A]")
    plt.title('RMSD')
    plt.legend()
    plt.xlim(0,120)
    plt.ylim(0,15)
    plt.title(bac_dict[item]['name'])
    plt.subplot(1,2,2)
    j=0
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict_long[item]['geom'][i]['0']['x_diam'][:bac_dict_long[item]['n'][i]]))*bac_dict_long[item]['maps'][i]['0']['vox'],bac_dict_long[item]['rmsd_avg'][i],yerr=bac_dict_long[item]['rmsd_std'][i],lw=2,label=i+" aa")
        plt.plot(np.arange(0,len(bac_dict_long[item]['geom'][i]['0']['x_diam'][:bac_dict_long[item]['n'][i]]))*bac_dict_long[item]['maps'][i]['0']['vox'],bac_dict_long[item]['rmsd'][i][0][:bac_dict_long[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict_long[item]['geom'][i]['1']['x_diam'][:bac_dict_long[item]['n'][i]]))*bac_dict_long[item]['maps'][i]['1']['vox'],bac_dict_long[item]['rmsd'][i][1][:bac_dict_long[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict_long[item]['geom'][i]['2']['x_diam'][:bac_dict_long[item]['n'][i]]))*bac_dict_long[item]['maps'][i]['2']['vox'],bac_dict_long[item]['rmsd'][i][2][:bac_dict_long[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.xlabel("Distance from the PTC along the Y-axis [A]")
    plt.ylabel("RMSD between central paths [A]")
    plt.title('RMSD')
    plt.legend()
    plt.xlim(0,120)
    plt.ylim(0,15)
    plt.title(bac_dict_long[item]['name'])
    # Save the figure as an SVG file
    plt.savefig(item+'/rmsd_plot_conv.svg', format='svg')
    plt.close()



##########
# Volume #
##########
plt.figure(figsize=(12,12))
plt.tight_layout()
for item in bac_id:
    plt.subplot(1,2,1)
    j=0
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom_avg'][i]['vol'],yerr = bac_dict[item]['geom_std'][i]['vol'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.xlabel("Distance from the PTC along the Y-axis [A]")
    plt.ylabel("Volume of the section of the exit tunne [$A^3$]")
    plt.legend()
    plt.xlim(0,120)
    plt.ylim(0,2500)
    plt.title(bac_dict[item]['name'])
    plt.subplot(1,2,2)
    j=0
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict_long[item]['geom'][i]['0']['z_diam'][:bac_dict_long[item]['n'][i]]))*bac_dict_long[item]['maps'][i]['0']['vox'],bac_dict_long[item]['geom_avg'][i]['vol'],yerr = bac_dict_long[item]['geom_std'][i]['vol'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict_long[item]['geom'][i]['0']['z_diam'][:bac_dict_long[item]['n'][i]]))*bac_dict_long[item]['maps'][i]['0']['vox'],bac_dict_long[item]['geom'][i]['0']['vol'][:bac_dict_long[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict_long[item]['geom'][i]['1']['z_diam'][:bac_dict_long[item]['n'][i]]))*bac_dict_long[item]['maps'][i]['1']['vox'],bac_dict_long[item]['geom'][i]['1']['vol'][:bac_dict_long[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict_long[item]['geom'][i]['2']['z_diam'][:bac_dict_long[item]['n'][i]]))*bac_dict_long[item]['maps'][i]['2']['vox'],bac_dict_long[item]['geom'][i]['2']['vol'][:bac_dict_long[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.xlabel("Distance from the PTC along the Y-axis [A]")
    plt.ylabel("Volume of the section of the exit tunne [$A^3$]")
    plt.legend()
    plt.xlim(0,120)
    plt.ylim(0,2500)
    plt.title(bac_dict_long[item]['name'])
    # Save the figure as an SVG file
    plt.savefig(item+'/volume_plot_conv.svg', format='svg')
    plt.close()

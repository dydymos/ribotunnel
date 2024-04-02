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
bac_id = ['4ybb','abau','bbur','bsub','cacn','drad','ecoli','efae','fjoh','linn','llac','lmon','mpne','msme','mtub','paer','pura','saur','tthe']

bac_name = ['E.coli','A.baumannii','B.burgdorferi','B.subtilis','C.acnes','D.radiodurans','E.coli','E.faecalis','F.johnsoniae','L.innocua','L.lactis','L.monocytogenes','M.pneumoniae','M.smegmatis','M.tuberculosis','P.aeruginosa','P.urativorans','S.aureus','T.thermophilus']

bac_dict = dict()
i=0
for item in bac_id:
    bac_dict[item] = dict()
    bac_dict[item]['name'] = bac_name[i]
    i+=1

################
# Getting maps #
################

keys = ['60','40','30','20','10']

for item in bac_id:
    bac_dict[item]['maps'] = dict()
    for i in keys:
        bac_dict[item]['maps'][i] = dict()
        bac_dict[item]['maps'][i]['0'] = get_map_param(item+'/'+i+"/average_final.mrc")
        bac_dict[item]['maps'][i]['1'] = get_map_param(item+'/'+i+"/ref_1/average_final.mrc")
        bac_dict[item]['maps'][i]['2'] = get_map_param(item+'/'+i+"/ref_2/average_final.mrc")


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


# Reading the main path file
for item in bac_id:
    bac_dict[item]['path'] = dict()
    for i in keys:
        bac_dict[item]['path'][i] = dict()
        bac_dict[item]['path'][i]['0'] = mda.Universe(item+'/'+i+'/main_path.pdb')
        bac_dict[item]['path'][i]['1'] = mda.Universe(item+'/'+i+'/ref_1/main_path.pdb')
        bac_dict[item]['path'][i]['2'] = mda.Universe(item+'/'+i+'/ref_2/main_path.pdb')



# Getting geometry
for item in bac_id:
    bac_dict[item]['geom'] = dict()
    for i in keys:
        bac_dict[item]['geom'][i] = dict()
        bac_dict[item]['geom'][i]['0'] = geometric_prop(bac_dict[item]['maps'][i]['0'],epsilon)
        bac_dict[item]['geom'][i]['1'] = geometric_prop(bac_dict[item]['maps'][i]['1'],epsilon)
        bac_dict[item]['geom'][i]['2'] = geometric_prop(bac_dict[item]['maps'][i]['2'],epsilon)


# number of atom in the main path
for item in bac_id:
    bac_dict[item]['n'] = dict()
    for i in keys:
        n_list = []
        for ii in bac_dict[item]['path'][i]:
            n_list.append(len(bac_dict[item]['geom'][i][ii]['z_diam']))
        bac_dict[item]['n'][i] = np.min(n_list)


# Averages and std errors
for item in bac_id:
    bac_dict[item]['geom_avg'] = dict()
    for i in keys:
        bac_dict[item]['geom_avg'][i] = dict()
        bac_dict[item]['geom_avg'][i]['z_diam'] = np.mean([bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_avg'][i]['x_diam'] = np.mean([bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_avg'][i]['vol'] = np.mean([bac_dict[item]['geom'][i]['0']['vol'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['vol'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['vol'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_avg'][i]['asp_rat'] = np.mean([bac_dict[item]['geom'][i]['0']['asp_rat'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['asp_rat'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['asp_rat'][:bac_dict[item]['n'][i]]],0)


for item in bac_id:
    bac_dict[item]['geom_std'] = dict()
    for i in keys:
        bac_dict[item]['geom_std'][i] = dict()
        bac_dict[item]['geom_std'][i]['z_diam'] = np.std([bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_std'][i]['x_diam'] = np.std([bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_std'][i]['vol'] = np.std([bac_dict[item]['geom'][i]['0']['vol'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['vol'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['vol'][:bac_dict[item]['n'][i]]],0)
        bac_dict[item]['geom_std'][i]['asp_rat'] = np.std([bac_dict[item]['geom'][i]['0']['asp_rat'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['1']['asp_rat'][:bac_dict[item]['n'][i]],bac_dict[item]['geom'][i]['2']['asp_rat'][:bac_dict[item]['n'][i]]],0)


# RMSD
for item in bac_id:
    bac_dict[item]['rmsd_avg'] = dict()
    bac_dict[item]['rmsd_std'] = dict()
    bac_dict[item]['rmsd'] = dict()
    for i in keys:
        bac_dict[item]['rmsd'][i] = []
        bac_dict[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict[item]['path'][i]['0'].atoms.positions[:bac_dict[item]['n'][i]]-bac_dict[item]['path'][i]['1'].atoms.positions[:bac_dict[item]['n'][i]])**2,1)))
        bac_dict[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict[item]['path'][i]['0'].atoms.positions[:bac_dict[item]['n'][i]]-bac_dict[item]['path'][i]['2'].atoms.positions[:bac_dict[item]['n'][i]])**2,1)))
        bac_dict[item]['rmsd'][i].append(np.sqrt(np.sum((bac_dict[item]['path'][i]['1'].atoms.positions[:bac_dict[item]['n'][i]]-bac_dict[item]['path'][i]['2'].atoms.positions[:bac_dict[item]['n'][i]])**2,1)))
        bac_dict[item]['rmsd_avg'][i] = np.mean(np.array(bac_dict[item]['rmsd'][i]),0)
        bac_dict[item]['rmsd_std'][i] = np.std(np.array(bac_dict[item]['rmsd'][i]),0)


############
# Plotting #
############
kolors = ['C0','C1','C2','C3','C4']

##################
### Pathway RMSD #
##################
for item in bac_id:
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
    # Save the figure as an SVG file
    plt.savefig(item+'/rmsd_plot.svg', format='svg')
    plt.close()

############################
# Pathways RMSD - ONE PLOT #
############################
plt.figure(figsize=(28,28))
plt.tight_layout()
k=1
for item in bac_id:
    j=0
    plt.subplot(5,4,k)
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['rmsd_avg'][i],yerr=bac_dict[item]['rmsd_std'][i],lw=2,label=i+" aa")
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['rmsd'][i][0][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['rmsd'][i][1][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['rmsd'][i][2][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.xlabel("Distance from the PTC along the Y-axis [A]",fontsize=14)
    plt.ylabel("RMSD between central paths [A]",fontsize=14)
    if (item=='4ybb'):plt.legend(fontsize=14)
    plt.xlim(0,120)
    plt.ylim(0,15)
    if (item=='4ybb'):
        plt.title('E.coli - X-ray',fontsize=18)
    else:
        plt.title(bac_dict[item]['name'],fontsize=18)
    # Save the figure as an SVG file
    k+=1

plt.savefig('rmsd_plot.svg', format='svg')
plt.close()

##########
# Volume #
##########
for item in bac_id:
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
    # Save the figure as an SVG file
    plt.savefig(item+'/volume_plot.svg', format='svg')
    plt.close()



######################
# Volume - ONE PLOT #
######################
plt.figure(figsize=(28,28))
plt.tight_layout()
k=1
for item in bac_id:
    j=0
    plt.subplot(5,4,k)
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['geom'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom_avg'][i]['vol'],yerr = bac_dict[item]['geom_std'][i]['vol'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.xlabel("Distance from the PTC along the Y-axis [A]",fontsize=14)
    plt.ylabel("Volume of the section of the exit tunne [$A^3$]",fontsize=14)
    if (item=='4ybb'):plt.legend(fontsize=14)
    plt.xlim(0,120)
    plt.ylim(0,2500)
    if (item=='4ybb'):
        plt.title('E.coli - X-ray',fontsize=18)
    else:
        plt.title(bac_dict[item]['name'],fontsize=18)
    k+=1

# Save the figure as an SVG file
plt.savefig('volume_plot.svg', format='svg')
plt.close()


###################
# Volume - LENGTH #
###################
plt.figure(figsize=(28,28))
plt.tight_layout()
k=5
j=0
for i in keys:
    plt.subplot(3,2,k)
    for item in bac_id:
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['vol'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
    plt.xlabel("Distance from the PTC along the Y-axis [A]",fontsize=14)
    plt.ylabel("Volume of the section of the exit tunne [$A^3$]",fontsize=14)
    plt.xlim(0,120)
    plt.ylim(0,2500)
    k-=1
    j+=1


# Save the figure as an SVG file
plt.savefig('volume_plot.svg', format='svg')
plt.close()

############################
# Distance from the X-Axis #
############################
for item in bac_id:
    j=0
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom_avg'][i]['x_diam'],yerr = bac_dict[item]['geom_std'][i]['x_diam'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.legend()
    plt.xlabel("Distance from the PTC along the Y-axis [A]")
    plt.ylabel("Distance from the main path along the X-axis [A]")
    plt.xlim(0,120)
    plt.ylim(0,60)
    plt.title(bac_dict[item]['name'])
    plt.savefig(item+'/x-plot.svg', format='svg')
    plt.close()


#######################################
# Distance from the X-Axis - ONE PLOT #
#######################################

plt.figure(figsize=(28,28))
plt.tight_layout()
k=1
for item in bac_id:
    j=0
    plt.subplot(5,4,k)
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom_avg'][i]['x_diam'],yerr = bac_dict[item]['geom_std'][i]['x_diam'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['x_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['x_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['x_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    if (item=='4ybb'):plt.legend(fontsize=14)
    plt.xlabel("Distance from the PTC along the Y-axis [A]",fontsize=14)
    plt.ylabel("Distance from the main path along the X-axis [A]",fontsize=14)
    plt.xlim(0,120)
    plt.ylim(0,65)
    if (item=='4ybb'):
        plt.title('E.coli - X-ray',fontsize=18)
    else:
        plt.title(bac_dict[item]['name'],fontsize=18)
    k+=1

# Save the figure as an SVG file
plt.savefig('x-plot.svg', format='svg')
plt.close()





############################
# Distance from the Z-Axis #
############################
for item in bac_id:
    j=0
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom_avg'][i]['z_diam'],yerr = bac_dict[item]['geom_std'][i]['z_diam'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.legend()
    plt.xlabel("Distance from the PTC along the Y-axis [A]")
    plt.ylabel("Distance from the main path along the Z-axis [A]")
    plt.xlim(0,120)
    plt.ylim(0,60)
    plt.title(bac_dict[item]['name'])
    plt.savefig(item+'/z-plot.svg', format='svg')
    plt.close()


#######################################
# Distance from the Z-Axis - ONE PLOT #
#######################################

plt.figure(figsize=(28,28))
plt.tight_layout()
k=1
for item in bac_id:
    j=0
    plt.subplot(5,4,k)
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom_avg'][i]['z_diam'],yerr = bac_dict[item]['geom_std'][i]['z_diam'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    if (item=='4ybb'):plt.legend(fontsize=14)
    plt.xlabel("Distance from the PTC along the Y-axis [A]",fontsize=14)
    plt.ylabel("Distance from the main path along the Z-axis [A]",fontsize=14)
    plt.xlim(0,120)
    plt.ylim(0,65)
    if (item=='4ybb'):
        plt.title('E.coli - X-ray',fontsize=18)
    else:
        plt.title(bac_dict[item]['name'],fontsize=18)
    k+=1

# Save the figure as an SVG file
plt.savefig('z-plot.svg', format='svg')
plt.close()

#####################
# Aspect ratio plot #
#####################
for item in bac_id:
    j=0
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom_avg'][i]['asp_rat'],yerr = bac_dict[item]['geom_std'][i]['asp_rat'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['asp_rat'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['asp_rat'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['asp_rat'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Aspect ratio X/Z dimension")
    plt.axhline(y=1.0, color='black',ls='--')
    plt.xlim(0,120)
    plt.ylim(0.25,3.0)
    plt.title(bac_dict[item]['name'])
    plt.legend()
    # Save the figure as an SVG file
    plt.savefig(item+'/aspect_ratio_plot.svg', format='svg')
    plt.close()


################################
# Aspect ratio plot - ONE PLOT #
################################
plt.figure(figsize=(28,28))
plt.tight_layout()
k=1
for item in bac_id:
    plt.subplot(5,4,k)
    j=0
    for i in keys:
        plt.errorbar(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom_avg'][i]['asp_rat'],yerr = bac_dict[item]['geom_std'][i]['asp_rat'],lw=2,label=i+' aa')
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['0']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['0']['vox'],bac_dict[item]['geom'][i]['0']['asp_rat'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['1']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['1']['vox'],bac_dict[item]['geom'][i]['1']['asp_rat'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        plt.plot(np.arange(0,len(bac_dict[item]['geom'][i]['2']['z_diam'][:bac_dict[item]['n'][i]]))*bac_dict[item]['maps'][i]['2']['vox'],bac_dict[item]['geom'][i]['2']['asp_rat'][:bac_dict[item]['n'][i]],lw=2,color=kolors[j],alpha=0.3)
        j+=1
    plt.xlabel("Distance from the PTC [A]",fontsize=14)
    plt.ylabel("Aspect ratio X/Z dimension",fontsize=14)
    plt.axhline(y=1.0, color='black',ls='--')
    plt.xlim(0,120)
    plt.ylim(0.25,3.0)
    if (item=='4ybb'):
        plt.title('E.coli - X-ray',fontsize=18)
    else:
        plt.title(bac_dict[item]['name'],fontsize=18)
    if (item=='4ybb'):plt.legend(fontsize=14)
    k+=1

# Save the figure as an SVG file
plt.savefig('aspect_ratio_plot.svg', format='svg')
plt.close()



# Correlation between maps
plt.tight_layout()
k=5
for i in keys:
    plt.subplot(3,2,k)
    id0 = np.where(bac_dict['tthe']['maps'][i]['0']['data']>epsilon)
    id1 = np.where(bac_dict['tthe']['maps'][i]['1']['data']>epsilon)
    id2 = np.where(bac_dict['tthe']['maps'][i]['2']['data']>epsilon)
    set_id0 = set(zip(id0[0], id0[1], id0[2]))
    set_id1 = set(zip(id1[0], id1[1], id1[2]))
    set_id2 = set(zip(id2[0], id2[1], id2[2]))
    common_id0_id1 = len(set_id0.intersection(set_id1))
    common_id0_id2 = len(set_id0.intersection(set_id2))
    common_id1_id2 = len(set_id1.intersection(set_id2))
    cc_matrix=np.ones([3,3])
    cc_matrix[0,1] = common_id0_id1/len(set_id0)
    cc_matrix[0,2] = common_id0_id2/len(set_id0)
    cc_matrix[1,0] = common_id0_id1/len(set_id1)
    cc_matrix[2,0] = common_id0_id2/len(set_id2)
    cc_matrix[1,2] = common_id1_id2/len(set_id1)
    cc_matrix[2,1] = common_id1_id2/len(set_id2)
    plt.imshow(cc_matrix,cmap='coolwarm_r',vmin=0.5, vmax=1.0)
    plt.title(i)
    plt.colorbar()
    k-=1

import numpy as np
import matplotlib.pyplot as plt
import mrcfile
import MDAnalysis as mda
import os


title = "A.baumannii"

# Function to derive all important parameters from the density map
def get_map_param(map_file):
    param = dict()
    map = mrcfile.open(map_file)
    param['data'] = map.data
    param['dim'] = np.shape(map.data)
    param['origin'] = np.array(map.header['origin'].tolist())
    param['vox'] = map.voxel_size.tolist()[0]
    return param

# Function to find a main path in the density that corresponds to the most dense region of the tunnel density along the Y-axis
def get_main_path(map_param,epsilon,output):
    indx = 1
    # PDB format template
    pdb_format = "ATOM  {:5d}  {:<4s}MOL     1    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           {}\n"
    atom_type = "C"
    # Output
    file = open(output,'w')
    # Y-axis cross-sections
    for i in range(0,map_param['dim'][1]):
        section = map_param['data'][:,i,:]
        if np.max(section) > epsilon:
            x_pos = np.where(section==np.max(section))[1][0]*map_param['vox'] + map_param['origin'][0]
            z_pos = np.where(section==np.max(section))[0][0]*map_param['vox'] + map_param['origin'][2]
            y_pos = i*map_param['vox'] + map_param['origin'][1]
            file.write(pdb_format.format(indx, atom_type, x_pos, y_pos, z_pos, atom_type))
            indx+=1
    file.close()

# Function that derives the most basic geometric properties of the exit tunnel
def geometric_prop(map_param,epsilon):
    geometry = dict()
    geometry["z_diam"] = []
    geometry["x_diam"] = []
    geometry["asp_rat"] = []
    geometry["vol"] = []
    for i in range(0,map_param['dim'][1]):
        section = map_param['data'][:,i,:]
        if np.max(section) > epsilon:
            geometry['vol'].append(len(np.where(section>epsilon)[0])*map_param['vox']*map_param['vox'])
            x_id = np.where(section==np.max(section))[0][0]
            z_id = np.where(section==np.max(section))[1][0]
            # Z-dim
            idx = z_id
            while idx>=0:
                if section[x_id,idx]>epsilon:idx-=1
                else:
                    z_len=z_id - idx
                    break
            idx = z_id+1
            while idx<=map_param['dim'][2]:
                if section[x_id,idx]>epsilon:idx+=1
                else:
                    z_len+= idx - z_id - 1
                    break
            geometry['z_diam'].append(z_len*map_param['vox'])
            # X-dim
            idx = x_id
            while idx>=0:
                if section[idx,z_id]>epsilon:idx-=1
                else:
                    x_len=x_id - idx
                    break
            idx = x_id+1
            while idx<=map_param['dim'][0]:
                if section[idx,z_id]>epsilon:idx+=1
                else:
                    x_len+= idx - x_id - 1
                    break
            geometry['x_diam'].append(x_len*map_param['vox'])
            geometry['asp_rat'].append(x_len/z_len)
    geometry["z_diam"] = np.array(geometry["z_diam"])
    geometry["x_diam"] = np.array(geometry["x_diam"])
    geometry['asp_rat'] = np.array(geometry['asp_rat'])
    geometry['vol'] = np.array(geometry['vol'])
    return geometry


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

import numpy as np
import matplotlib.pyplot as plt
import mrcfile
import MDAnalysis as mda
import os

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

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.ndimage import gaussian_filter1d

################################
# Plotting 3D plot of pathways #
################################
def paths_plot3D(org_dict,name,item):
    fig = plt.figure()
    path_coords = org_dict[name]['res occ path'][item]['2']
    x = np.array(path_coords).T[0]
    y = np.array(path_coords).T[1]
    z = np.array(path_coords).T[2]
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
    ax.scatter(np.array(org_dict[name]['res occ path'][item]['0'])[:,0], np.array(org_dict[name]['res occ path'][item]['0'])[:,1], np.array(org_dict[name]['res occ path'][item]['0'])[:,2])
    ax.scatter(np.array(org_dict[name]['res occ path'][item]['1'])[:,0], np.array(org_dict[name]['res occ path'][item]['1'])[:,1], np.array(org_dict[name]['res occ path'][item]['1'])[:,2])
    ax.scatter(np.array(org_dict[name]['res occ path'][item]['2'])[:,0], np.array(org_dict[name]['res occ path'][item]['2'])[:,1], np.array(org_dict[name]['res occ path'][item]['2'])[:,2])
    ax.scatter(np.array(org_dict[name]['res occ path avg'][item])[:,0], np.array(org_dict[name]['res occ path avg'][item])[:,1], np.array(org_dict[name]['res occ path avg'][item])[:,2],color="black")
    ax.plot(org_dict[name]['spline path'][item]['0'][0], org_dict[name]['spline path'][item]['0'][1], org_dict[name]['spline path'][item]['0'][2], lw=2, label='The 1st path - spline')
    ax.plot(org_dict[name]['spline path'][item]['1'][0], org_dict[name]['spline path'][item]['1'][1], org_dict[name]['spline path'][item]['1'][2], lw=2, label='The 2nd path - spline')
    ax.plot(org_dict[name]['spline path'][item]['2'][0], org_dict[name]['spline path'][item]['2'][1], org_dict[name]['spline path'][item]['2'][2], lw=2, label='The 3rd path - spline')
    ax.plot(org_dict[name]['spline path avg'][item][0], org_dict[name]['spline path avg'][item][1], org_dict[name]['spline path avg'][item][2], lw=2, color="black", label='The average path - spline')
    # Set the calculated limits for each axis
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_zlim(z_limits)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()
    plt.title("Ribosome from "+org_dict[name]['name']+" at "+item+" length")
    plt.show()


###################################################
# Function to plot total properties of the tunnel #
###################################################
def plot_tunnel_prop(org_dict,org_id,output_name,short=False):
    fig = plt.figure(figsize=(15, 10))
    proteins_prop = dict()
    if (short==False):
        lista = ["10", "20", "30", "40", "60"]
    if (short==True):
        lista = ["10", "20", "30", "40"]
    for item in lista:
        proteins_prop[item] = np.zeros([len(org_id),5])
    i=0
    max_proteins = 0
    for name in org_id:
        for item in lista:
            proteins_prop[item][i] = org_dict[name]["proteins"][item]
            if np.max(org_dict[name]["proteins"][item]) > max_proteins: 
                max_proteins = np.max(org_dict[name]["proteins"][item])
        i+=1    
    j=1
    for item in lista:
        plt.subplot(2,3,j)
        sns.violinplot(x=np.full(len(org_id), "polar"),y=proteins_prop[item][:,0],color="pink",inner=None,linewidth=1.5,edgecolor="pink",alpha = 0.2)
        sns.stripplot(x=np.full(len(org_id), "polar"),y=proteins_prop[item][:,0],color="pink",size=6)
        sns.violinplot(x=np.full(len(org_id), "hydrophobic"),y=proteins_prop[item][:,1],color="orange",inner=None,linewidth=1.5,edgecolor="orange",alpha = 0.2)
        sns.stripplot(x=np.full(len(org_id), "hydrophobic"),y=proteins_prop[item][:,1],color="orange",size=6)
        sns.violinplot(x=np.full(len(org_id), "positive"),y=proteins_prop[item][:,2],color="red",inner=None,linewidth=1.5,edgecolor="red",alpha = 0.2)
        sns.stripplot(x=np.full(len(org_id), "positive"),y=proteins_prop[item][:,2],color="red",size=6)
        sns.violinplot(x=np.full(len(org_id), "negative"),y=proteins_prop[item][:,3],color="blue",inner=None,linewidth=1.5,edgecolor="blue",alpha = 0.2)
        sns.stripplot(x=np.full(len(org_id), "negative"),y=proteins_prop[item][:,3],color="blue",size=6)
        sns.violinplot(x=np.full(len(org_id), "special"),y=proteins_prop[item][:,4],color="silver",inner=None,linewidth=1.5,edgecolor="silver",alpha = 0.2)
        sns.stripplot(x=np.full(len(org_id), "special"),y=proteins_prop[item][:,4],color="silver",size=6)
        plt.ylim(-1,45)
        plt.ylabel("Number of AA")
        plt.title("NC length "+item+" AA")
        j+=1
    plt.savefig(output_name)
    plt.close()

#######################################################
# Function to plot tunnel properties along its length #
#######################################################
def path_prop(org_dict,org_id,short=False):
    polar = []
    hydro = [] 
    negat = []
    posit = []
    speci = []
    rna = []
    length_p = []
    if (short==False):
        for name in org_id:
            polar.append(org_dict[name]["proteins_matrix"]['60'][:,0])
            hydro.append(org_dict[name]["proteins_matrix"]['60'][:,1])
            posit.append(org_dict[name]["proteins_matrix"]['60'][:,2])
            negat.append(org_dict[name]["proteins_matrix"]['60'][:,3])
            speci.append(org_dict[name]["proteins_matrix"]['60'][:,4])
            rna.append(org_dict[name]["rna_matrix"]['60'])
            length_p.append(org_dict[name]['path length']['60'][:-1])
    if (short==True):
        for name in org_id:
            polar.append(org_dict[name]["proteins_matrix"]['40'][:,0])
            hydro.append(org_dict[name]["proteins_matrix"]['40'][:,1])
            posit.append(org_dict[name]["proteins_matrix"]['40'][:,2])
            negat.append(org_dict[name]["proteins_matrix"]['40'][:,3])
            speci.append(org_dict[name]["proteins_matrix"]['40'][:,4])
            rna.append(org_dict[name]["rna_matrix"]['40'])
            length_p.append(org_dict[name]['path length']['40'][:-1])
    return polar,hydro,posit,negat,speci,rna,length_p


def plot_prop_dist(org_dict,org_id,output_name,end,short=False):
    plt.figure(figsize=(15, 10))
    if (short==False):
        polar,hydro,posit,negat,speci,rna,length = path_prop(org_dict,org_id,short=False)
        if (end == "NaN"): cut = 0
        else:
            cut = np.where(length[0]==end)[0][0]
    if (short==True):
        polar,hydro,posit,negat,speci,rna,length = path_prop(org_dict,org_id,short=True)
        cut = 0
    protein = []
    for i in range(0,len(org_id)):
        protein.append(polar[i]+hydro[i]+posit[i]+negat[i]+speci[i])
    plt.subplot(2,3,1)
    plt.plot(length[0][cut:],np.mean(polar,0)[cut:],color='pink',lw = 3)
    plt.fill_between(length[0][cut:],np.mean(polar,0)[cut:]-np.std(polar,0)[cut:],np.mean(polar,0)[cut:]+np.std(polar,0)[cut:],color='pink',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")
    plt.title("Polar")
    plt.subplot(2,3,2)
    plt.plot(length[0][cut:],np.mean(hydro,0)[cut:],color='orange',lw = 3)
    plt.fill_between(length[0][cut:],np.mean(hydro,0)[cut:]-np.std(hydro,0)[cut:],np.mean(hydro,0)[cut:]+np.std(hydro,0)[cut:],color='orange',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")
    plt.title("Hydrophobic")
    plt.subplot(2,3,3)
    plt.plot(length[0][cut:],np.mean(posit,0)[cut:],color='red',lw = 3)
    plt.fill_between(length[0][cut:],np.mean(posit,0)[cut:]-np.std(posit,0)[cut:],np.mean(posit,0)[cut:]+np.std(posit,0)[cut:],color='red',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")
    plt.title("Positively charged")
    plt.subplot(2,3,4)
    plt.plot(length[0][cut:],np.mean(negat,0)[cut:],color='blue',lw = 3)
    plt.fill_between(length[0][cut:],np.mean(negat,0)[cut:]-np.std(negat,0)[cut:],np.mean(negat,0)[cut:]+np.std(negat,0)[cut:],color='blue',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")  
    plt.title("Negatively charged")
    plt.subplot(2,3,5)
    plt.plot(length[0][cut:],np.mean(speci,0)[cut:],color='silver',lw = 3)
    plt.fill_between(length[0][cut:],np.mean(speci,0)[cut:]-np.std(speci,0)[cut:],np.mean(speci,0)[cut:]+np.std(speci,0)[cut:],color='silver',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")
    plt.title("Special")
    plt.subplot(2,3,6)
    plt.plot(length[0][cut:],np.mean(rna,0)[cut:],color='black',lw = 3)
    plt.fill_between(length[0][cut:],np.mean(rna,0)[cut:]-np.std(rna,0)[cut:],np.mean(rna,0)[cut:]+np.std(rna,0)[cut:],color='black',alpha = 0.5)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of nucleotides")
    plt.title("RNA")
    plt.xlim(0,100)
    plt.ylim(0,27)
    plt.savefig(output_name+"_length.svg")
    plt.close()
    plt.figure()
    total = np.mean(protein,0)+np.mean(rna,0)
    ratio = np.mean(rna,0)/total
    ratio_std = np.sqrt((np.mean(protein,0) / total**2)**2 * np.std(rna,0)**2 + (np.mean(rna,0) / total**2)**2 * np.std(protein,0)**2)
    plt.plot(length[0][cut:],ratio[cut:],color="black",lw = 3)
    plt.axhline(y=0.5,ls="--")
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("RNA vs Protein ratio")
    plt.ylim(0,1.1)
    plt.xlim(0,100)
    plt.savefig(output_name+"_ratio.svg")
    np.savetxt(output_name+"_ratio.dat",ratio)
    np.savetxt(output_name+"_ratio_std.dat",ratio_std)
    plt.close()

def plot_prop_dist_old(org_dict,org_id,output_name):
    plt.figure(figsize=(15, 10))
    polar,hydro,posit,negat,speci,rna,length = path_prop(org_dict,org_id)
    protein = []
    for i in range(0,len(org_id)):
        protein.append(polar[i]+hydro[i]+posit[i]+negat[i]+speci[i])
    plt.subplot(2,3,1)
    plt.plot(length[0],np.mean(polar,0),color='pink',lw = 3)
    plt.fill_between(length[0],np.mean(polar,0)-np.std(polar,0),np.mean(polar,0)+np.std(polar,0),color='pink',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")
    plt.title("Polar")
    plt.subplot(2,3,2)
    plt.plot(length[0],np.mean(hydro,0),color='orange',lw = 3)
    plt.fill_between(length[0],np.mean(hydro,0)-np.std(hydro,0),np.mean(hydro,0)+np.std(hydro,0),color='orange',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")
    plt.title("Hydrophobic")
    plt.subplot(2,3,3)
    plt.plot(length[0],np.mean(posit,0),color='red',lw = 3)
    plt.fill_between(length[0],np.mean(posit,0)-np.std(posit,0),np.mean(posit,0)+np.std(posit,0),color='red',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")
    plt.title("Positively charged")
    plt.subplot(2,3,4)
    plt.plot(length[0],np.mean(negat,0),color='blue',lw = 3)
    plt.fill_between(length[0],np.mean(negat,0)-np.std(negat,0),np.mean(negat,0)+np.std(negat,0),color='blue',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")  
    plt.title("Negatively charged")
    plt.subplot(2,3,5)
    plt.plot(length[0],np.mean(speci,0),color='silver',lw = 3)
    plt.fill_between(length[0],np.mean(speci,0)-np.std(speci,0),np.mean(speci,0)+np.std(speci,0),color='silver',alpha = 0.5)
    plt.ylim(0,10)
    plt.xlim(0,100)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of amino acids")
    plt.title("Special")
    plt.subplot(2,3,6)
    plt.plot(length[0],np.mean(rna,0),color='black',lw = 3)
    plt.fill_between(length[0],np.mean(rna,0)-np.std(rna,0),np.mean(rna,0)+np.std(rna,0),color='black',alpha = 0.5)
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("Number of nucleotides")
    plt.title("RNA")
    plt.xlim(0,100)
    plt.ylim(0,27)
    plt.savefig(output_name+"_length.svg")
    plt.close()
    plt.figure()
    ratio = np.mean(rna,0)/(np.mean(protein,0)+np.mean(rna,0))
    plt.plot(length[0],ratio,color="black",lw = 3)
    plt.axhline(y=0.5,ls="--")
    plt.xlabel("Distance from the PTC [A]")
    plt.ylabel("RNA vs Protein ratio")
    plt.ylim(0,1.1)
    plt.xlim(0,100)
    plt.savefig(output_name+"_ratio.svg")
    plt.close()

###########################################
# Functions to plot tunnel cross sections #
###########################################
def plot_avg_cross_section(org_dict,org_id,output_name,end,short=False):
    fig = plt.figure(figsize=(15, 10))
    kolors = ['C0','C1','C2','C3','C4']
    j=0
    if (short==False):
        tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
        lista = ["10", "20", "30", "40","60"]
        if (end=="NaN"): cut = 0
        else:
            cut = np.where(tunnel_length==end)[0][0]
    if (short==True):
        tunnel_length = org_dict[org_id[0]]['path length']['40'][:-1]
        lista = ["10", "20", "30", "40"]
        cut=0
    for item in lista:
        path_matrix = np.zeros([len(org_id),tunnel_length.size,])
        for i,org in enumerate(org_id):
            path_matrix[i] = org_dict[org]['cross section'][item]
        plt.subplot(2,3,j+1)
        plt.plot(tunnel_length[cut:],np.nanmean(path_matrix,0)[cut:],lw=3,color=kolors[j])
        plt.fill_between(tunnel_length[cut:],np.nanmean(path_matrix,0)[cut:]-np.nanstd(path_matrix,0)[cut:],np.nanmean(path_matrix,0)[cut:]+np.nanstd(path_matrix,0)[cut:],color=kolors[j],alpha=0.5)
        plt.title("N = "+item +" amino acids")
        plt.xlim(0,100)
        if (short==False and item == "60"):
            plt.axvline(x=92,ls='--',color='grey')
            plt.ylim(0,2500)
        elif (short==False and item == "40"):
            plt.axvline(x=92,ls='--',color='grey')
            plt.ylim(0,1500)
        else:
            plt.ylim(0,1500)
        j+=1
    j=0
    plt.subplot(2,3,6)
    for item in lista:
        path_matrix = np.zeros([len(org_id),tunnel_length.size,])
        for i,org in enumerate(org_id):
            path_matrix[i] = org_dict[org]['cross section'][item]
        plt.plot(tunnel_length[cut:],np.nanmean(path_matrix,0)[cut:],lw=3,color=kolors[j])
        plt.fill_between(tunnel_length[cut:],np.nanmean(path_matrix,0)[cut:]-np.nanstd(path_matrix,0)[cut:],np.nanmean(path_matrix,0)[cut:]+np.nanstd(path_matrix,0)[cut:],color=kolors[j],alpha=0.5)
        plt.title("All lengths")
        if (short==False):
            plt.axvline(x=92,ls='--',color='grey')
        plt.xlim(0,100)        
        plt.ylim(0,1500)
        j+=1
    fig.text(0.5, 0.05, 'Distance from the PTC [A]', ha='center', va='center', fontsize=14)
    fig.text(0.05, 0.5, r'Cross section surface $[A^2]$', ha='center', va='center', rotation='vertical', fontsize=14)
    plt.savefig(output_name)
    plt.close()


def plot_avg_cross_section_old(org_dict,org_id,output_name):
    fig = plt.figure(figsize=(15, 10))
    kolors = ['C0','C1','C2','C3','C4']
    j=0
    tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
    for item in ["10", "20", "30", "40","60"]:
        path_matrix = np.zeros([len(org_id),tunnel_length.size,])
        for i,org in enumerate(org_id):
            path_matrix[i] = org_dict[org]['cross section'][item]
        plt.subplot(2,3,j+1)
        plt.plot(tunnel_length,np.nanmean(path_matrix,0),lw=3,color=kolors[j])
        plt.fill_between(tunnel_length,np.nanmean(path_matrix,0)-np.nanstd(path_matrix,0),np.nanmean(path_matrix,0)+np.nanstd(path_matrix,0),color=kolors[j],alpha=0.5)
        plt.title("N = "+item +" amino acids")
        plt.xlim(0,100)
        if item == "60":
            plt.axvline(x=92,ls='--',color='grey')
            plt.ylim(0,2500)
        elif item == "40":
            plt.axvline(x=92,ls='--',color='grey')
            plt.ylim(0,1500)
        else:
            plt.ylim(0,1500)
        j+=1
    j=0
    plt.subplot(2,3,6)
    for item in ["10", "20", "30", "40","60"]:
        path_matrix = np.zeros([len(org_id),tunnel_length.size,])
        for i,org in enumerate(org_id):
            path_matrix[i] = org_dict[org]['cross section'][item]
        plt.plot(tunnel_length,np.nanmean(path_matrix,0),lw=3,color=kolors[j])
        plt.fill_between(tunnel_length,np.nanmean(path_matrix,0)-np.nanstd(path_matrix,0),np.nanmean(path_matrix,0)+np.nanstd(path_matrix,0),color=kolors[j],alpha=0.5)
        plt.title("All lengths")
        plt.axvline(x=92,ls='--',color='grey')
        plt.xlim(0,100)        
        plt.ylim(0,1500)
        j+=1
    fig.text(0.5, 0.05, 'Distance from the PTC [A]', ha='center', va='center', fontsize=14)
    fig.text(0.05, 0.5, r'Cross section surface $[A^2]$', ha='center', va='center', rotation='vertical', fontsize=14)
    plt.savefig(output_name)
    plt.close()

# For nice figure separation
def find_closest_factors(X):
    # Start with the square root of X to get n and m close to each other
    n = int(np.sqrt(X))
    # Increase n until n * m >= X
    while n * np.ceil(X / n) < X:
        n += 1
    m = int(np.ceil(X / n))
    return n, m

def plot_cross_section(org_dict,org_id,output_name,end):
    fig = plt.figure(figsize=(15, 10))
    kolors = ['C0','C1','C2','C3','C4']
    # Plot cross sections surface
    i=1
    n,m = find_closest_factors(len(org_id))
    tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
    if (end=="NaN"): cut = 0
    else:
        cut = np.where(tunnel_length==end)[0][0]
    for name in org_id:
        plt.subplot(n,m,i)
        j = 0
        for item in ["10", "20", "30", "40","60"]:
            # Extract data
            x = tunnel_length[cut:]
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            if i == 1:
                plt.plot(x,y_smooth[cut:],color=kolors[j],label="NC length "+item,lw=2)
            else:
                plt.plot(x,y_smooth[cut:],color=kolors[j],lw=2)
            j+=1
        plt.title(org_dict[name]['name'])
        plt.xlim(0,100)
        plt.ylim(0,1500)      
        i+=1
    plt.subplots_adjust(hspace=0.5, wspace=0.4)
    # Add the global legend
    fig.legend(loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0.92), fontsize=12)
    fig.text(0.5, 0.05, 'Distance from the PTC [A]', ha='center', va='center', fontsize=14)
    fig.text(0.05, 0.5, r'Cross section surface $[A^2]$', ha='center', va='center', rotation='vertical', fontsize=14)
    plt.savefig(output_name)
    plt.close()

def plot_cross_section_old(org_dict,org_id,output_name):
    fig = plt.figure(figsize=(15, 10))
    kolors = ['C0','C1','C2','C3','C4']
    # Plot cross sections surface
    i=1
    n,m = find_closest_factors(len(org_id))
    tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
    for name in org_id:
        plt.subplot(n,m,i)
        j = 0
        for item in ["10", "20", "30", "40","60"]:
            # Extract data
            x = tunnel_length
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            if i == 1:
                plt.plot(x,y_smooth,color=kolors[j],label="NC length "+item,lw=2)
            else:
                plt.plot(x,y_smooth,color=kolors[j],lw=2)
            j+=1
        plt.title(org_dict[name]['name'])
        plt.xlim(0,100)
        plt.ylim(0,1500)      
        i+=1
    plt.subplots_adjust(hspace=0.5, wspace=0.4)
    # Add the global legend
    fig.legend(loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0.92), fontsize=12)
    fig.text(0.5, 0.05, 'Distance from the PTC [A]', ha='center', va='center', fontsize=14)
    fig.text(0.05, 0.5, r'Cross section surface $[A^2]$', ha='center', va='center', rotation='vertical', fontsize=14)
    plt.savefig(output_name)
    plt.close()

def plot_cross_section_length(org_dict,org_id,output_name,end):
    fig = plt.figure(figsize=(15, 10))
    kolors = ['C0','C1','C2','C3','C4']
    tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
    if (end=="NaN"): cut = 0
    else:
        cut = np.where(tunnel_length==end)[0][0]
    j = 0   
    for item in ["10", "20", "30", "40","60"]:
        plt.subplot(2,3,j+1)
        for name in org_id:
            # Extract data
            x = tunnel_length[cut:]
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            if name in ["ecun",'ploc','vnec','slop']:
                plt.plot(x,y_smooth[cut:],color="black",lw=2)
            else:
                plt.plot(x,y_smooth[cut:],color=kolors[j],lw=2)
            plt.xlim(0,100)
            if item == "60":
                plt.axvline(x=92,ls='--',color='grey') 
                plt.ylim(0,2500)
            elif item == "40":
                plt.axvline(x=92,ls='--',color='grey')
                plt.ylim(0,1500)
            else: plt.ylim(0,1500)
            plt.title("N = "+item +" amino acids")
        j+=1
    j=0    
    plt.subplot(2,3,6)
    for item in ["10", "20", "30", "40","60"]:
        for name in org_id:
            # Extract data
            x = tunnel_length[cut:]
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            if name in ["ecun",'ploc','vnec','slop']:
                plt.plot(x,y_smooth[cut:],color="black",lw=2)
            else:
                plt.plot(x,y_smooth[cut:],color=kolors[j],lw=2)
            plt.axvline(x=92,ls='--',color='grey') 
            plt.xlim(0,100)
            plt.ylim(0,1500)
            plt.title("All lengths")
        j+=1
    fig.text(0.5, 0.05, 'Distance from the PTC [A]', ha='center', va='center', fontsize=14)
    fig.text(0.05, 0.5, r'Cross section surface $[A^2]$', ha='center', va='center', rotation='vertical', fontsize=14)
    plt.savefig(output_name)
    plt.close()

def plot_cross_section_length_old(org_dict,org_id,output_name):
    fig = plt.figure(figsize=(15, 10))
    kolors = ['C0','C1','C2','C3','C4']
    tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
    j = 0   
    for item in ["10", "20", "30", "40","60"]:
        plt.subplot(2,3,j+1)
        for name in org_id:
            # Extract data
            x = tunnel_length
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            if name in ["ecun",'ploc','vnec','slop']:
                plt.plot(x,y_smooth,color="black",lw=2)
            else:
                plt.plot(x,y_smooth,color=kolors[j],lw=2)
            plt.xlim(0,100)
            if item == "60":
                plt.axvline(x=92,ls='--',color='grey') 
                plt.ylim(0,2500)
            elif item == "40":
                plt.axvline(x=92,ls='--',color='grey')
            else: plt.ylim(0,1500)
            plt.title("N = "+item +" amino acids")
        j+=1
    j=0    
    plt.subplot(2,3,6)
    for item in ["10", "20", "30", "40","60"]:
        for name in org_id:
            # Extract data
            x = tunnel_length
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            if name in ["ecun",'ploc','vnec','slop']:
                plt.plot(x,y_smooth,color="black",lw=2)
            else:
                plt.plot(x,y_smooth,color=kolors[j],lw=2)
            plt.axvline(x=92,ls='--',color='grey') 
            plt.xlim(0,100)
            plt.ylim(0,1500)
            plt.title("All lengths")
        j+=1
    fig.text(0.5, 0.05, 'Distance from the PTC [A]', ha='center', va='center', fontsize=14)
    fig.text(0.05, 0.5, r'Cross section surface $[A^2]$', ha='center', va='center', rotation='vertical', fontsize=14)
    plt.savefig(output_name)
    plt.close()


def plot_cross_section_length_special(org_dict,org_id,output_name,end):
    fig = plt.figure(figsize=(15, 10))
    kolors = ['C0','C1','C2','C3','C4']
    tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
    if (end=="NaN"): cut = 0
    else: 
        cut = np.where(tunnel_length==end)[0][0]
    j = 0   
    for item in ["10", "20", "30", "40","60"]:
        plt.subplot(2,3,j+1)
        if item == "20": 
            cluster_1 = ["tthe","msme","cacn","efae"]
            cluster_2 = ["saur","linn","paer","abau","ecoli"]
            cluster_3 = ["mtub","llac","lmon","mpne","bbur","drad","bsub","pura","fjoh"]
        else: 
            cluster_1 = []
            cluster_2 = []
            cluster_3 = org_id
        for name in cluster_3:
            # Extract data
            x = tunnel_length[cut:]
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            plt.plot(x,y_smooth[cut:],color=kolors[j],lw=2)
            plt.xlim(0,100)
            if item == "60":
                plt.axvline(x=92,ls='--',color='grey') 
                plt.ylim(0,2500)
            elif item == "40":
                plt.axvline(x=92,ls='--',color='grey')
                plt.ylim(0,1500)
            else: plt.ylim(0,1500)
            plt.title("N = "+item +" amino acids")            
        for name in cluster_1:
            # Extract data
            x = tunnel_length[cut:]
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            plt.plot(x,y_smooth[cut:],color="black",lw=2,label=name)           
            plt.xlim(0,100)
            if item == "60":
                plt.axvline(x=92,ls='--',color='grey') 
                plt.ylim(0,2500)
            elif item == "40":
                plt.axvline(x=92,ls='--',color='grey')
                plt.ylim(0,1500)
            else: plt.ylim(0,1500)
            plt.title("N = "+item +" amino acids")
            plt.legend()
        for name in cluster_2:
            # Extract data
            x = tunnel_length[cut:]
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            plt.plot(x,y_smooth[cut:],color="red",lw=2,label=name)           
            plt.xlim(0,100)
            if item == "60":
                plt.axvline(x=92,ls='--',color='grey') 
                plt.ylim(0,2500)
            elif item == "40":
                plt.axvline(x=92,ls='--',color='grey')
                plt.ylim(0,1500)
            else: plt.ylim(0,1500)
            plt.title("N = "+item +" amino acids")
        j+=1
    j=0    
    plt.subplot(2,3,6)
    for item in ["10", "20", "30", "40","60"]:
        for name in org_id:
            # Extract data
            x = tunnel_length[cut:]
            y = np.array(org_dict[name]['cross section'][item])
            # Apply Gaussian filter
            sigma = 1  # Standard deviation for Gaussian kernel; adjust as needed
            y_smooth = gaussian_filter1d(y, sigma)
            if name in ["ecun",'ploc','vnec','slop']:
                plt.plot(x,y_smooth[cut:],color="black",lw=2)
            else:
                plt.plot(x,y_smooth[cut:],color=kolors[j],lw=2)
            plt.axvline(x=92,ls='--',color='grey') 
            plt.xlim(0,100)
            plt.ylim(0,1500)
            plt.title("All lengths")
        j+=1
    fig.text(0.5, 0.05, 'Distance from the PTC [A]', ha='center', va='center', fontsize=14)
    fig.text(0.05, 0.5, r'Cross section surface $[A^2]$', ha='center', va='center', rotation='vertical', fontsize=14)
    plt.savefig(output_name)
    plt.close()

##### ALL PLOTS ######
#fig = plt.figure(figsize=(15, 10))
#tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
#cut = np.where(tunnel_length==end)[0][0]
#j = 0
#for item in ["10","20","30","40","60"]:
#    plt.subplot(2,3,j+1)
#    path_matrix_b = np.zeros([len(bac_id),tunnel_length.size,])
#    for i,org in enumerate(bac_id):
#        path_matrix_b[i] = bac_dict[org]['cross section'][item]
    #path_matrix_e = np.zeros([len(euk_core_id),tunnel_length.size,])
    #for i,org in enumerate(euk_core_id):
    #    path_matrix_e[i] = euk_core_dict[org]['cross section'][item]
    #path_matrix_a = np.zeros([len(arc_id),tunnel_length.size,])
    #for i,org in enumerate(arc_id):
    #    path_matrix_a[i] = arc_dict[org]['cross section'][item]
#    path_matrix_c = np.zeros([len(chloro_id),tunnel_length.size,])
#    for i,org in enumerate(chloro_id):
#        path_matrix_c[i] = chloro_dict[org]['cross section'][item]
#    plt.subplot(2,3,j+1)
#    plt.plot(tunnel_length[cut:],np.nanmean(path_matrix_b,0)[cut:],lw=3,color="#77DD77",label="Bacteria")
    #plt.plot(tunnel_length[cut:],np.nanmean(path_matrix_e,0)[cut:],lw=3,color="#4A90E2",label="Eukaryota")
    #plt.plot(tunnel_length[cut:],np.nanmean(path_matrix_a,0)[cut:],lw=3,color="#FF6B6B",label="Archaea")
#    plt.plot(tunnel_length[cut:],np.nanmean(path_matrix_c,0)[cut:],lw=3,color="orange",label="Chloroplasts")
#    plt.fill_between(tunnel_length[cut:],np.nanmean(path_matrix_b,0)[cut:]-np.nanstd(path_matrix_b,0)[cut:],np.nanmean(path_matrix_b,0)[cut:]+np.nanstd(path_matrix_b,0)[cut:],color="#77DD77",alpha=0.5)
    #plt.fill_between(tunnel_length[cut:],np.nanmean(path_matrix_e,0)[cut:]-np.nanstd(path_matrix_e,0)[cut:],np.nanmean(path_matrix_e,0)[cut:]+np.nanstd(path_matrix_e,0)[cut:],color="#4A90E2",alpha=0.5)
    #plt.fill_between(tunnel_length[cut:],np.nanmean(path_matrix_a,0)[cut:]-np.nanstd(path_matrix_a,0)[cut:],np.nanmean(path_matrix_a,0)[cut:]+np.nanstd(path_matrix_a,0)[cut:],color="#FF6B6B",alpha=0.5)
#    plt.title("N = "+item +" amino acids")
#    plt.xlim(0,100)
#    if item == "60":
#        plt.axvline(x=92,ls='--',color='grey')
#        plt.legend()
#        plt.ylim(0,2300)
#    elif item == "40":
#        plt.axvline(x=92,ls='--',color='grey')
#        plt.ylim(0,1300)
#    else:
#        plt.ylim(0,1300)
#    j+=1

#fig.text(0.5, 0.05, 'Distance from the PTC [A]', ha='center', va='center', fontsize=14)
#fig.text(0.05, 0.5, r'Cross section surface $[A^2]$', ha='center', va='center', rotation='vertical', fontsize=14)
#plt.savefig("bac_chloro.svg")

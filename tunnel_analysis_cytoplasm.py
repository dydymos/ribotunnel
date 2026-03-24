import os
import numpy as np
import functions
import plotting
from functions import get_dictionary, read_grid, read_occupancy_path, calculate_spline, write_paths, write_tunnel_path, write_tunnel_asph, get_tunnel, get_edge_from_occupancy, write_tunnel_prop, get_rings, get_discs#, path_interp
from plotting import paths_plot3D, plot_tunnel_prop, plot_prop_dist, plot_avg_cross_section, plot_cross_section, plot_cross_section_length

################
# Ribosome IDs #
################

# "4ybb" - X-ray E.coli
bac_id = ['abau','bbur','bsub','drad','ecoli','efae','fjoh','linn','llac','lmon','mpne','msme','mtub','paer','pura','saur','tthe']
euk_id = ['calb','cele','cthe','dmel','drer','ecun','egra','ggal','glam','hsap','klac','ldon','lmaj','mmus','ncra','ntab','ocun','pfal','ploc','rnor','scer','slop','slyc','spom','sscr','taes','tbru','tcru','tgon','tthe','tvag','vnec','xlae']
euk_core_id = ['calb','cele','cthe','dmel','drer','egra','ggal','glam','hsap','klac','ldon','lmaj','mmus','ncra','ntab','ocun','pfal','rnor','scer','slyc','spom','sscr','taes','tbru','tcru','tgon','tthe','tvag','xlae']
arc_id = ['hmar','pfur','saci','tkod','pcal']
micro_id = ['ecun','ploc','slop','vnec']

# ribosome names
bac_name = ['A.baumannii','B.burgdorferi','B.subtilis','D.radiodurans','E.coli','E.faecalis','F.johnsoniae','L.innocua','L.lactis','L.monocytogenes','M.pneumoniae','M.smegmatis','M.tuberculosis','P.aeruginosa','P.urativorans','S.aureus','T.thermophilus']
arc_name = ["H.marismorturi","P.furiosus","S.acidocaldarius","T.kodakarensis",'P.calidifontis']
euk_name = ["C.albicans", "C.elegans", "C.thermophilia", "D.melanogaster", "D.rerio", "E.cuniculi", "E.gracilis", "G.gallus","G.lamblia", "H.sapiens", "K.lactis", "L.donovani", "L.major", "M.musculus", "N.crasa", "N.tabacum", "O.cuniculus", "P.falciparum", "P.locustae", "R.norvegicus", "S.cerevisiae", "S.lophii", "S.lycoperscium", "S.pombe", "S.scrofa", "T.aestivum", "T.brucei", "T.cruzi", "T.gondii", "T.thermophila","T.vaginalis","V.necatrix","X.laevis"]
euk_core_name = ["C.albicans", "C.elegans", "C.thermophilia", "D.melanogaster", "D.rerio", "E.gracilis", "G.gallus","G.lamblia", "H.sapiens", "K.lactis", "L.donovani", "L.major", "M.musculus", "N.crasa", "N.tabacum", "O.cuniculus", "P.falciparum", "R.norvegicus", "S.cerevisiae", "S.lycoperscium", "S.pombe", "S.scrofa", "T.aestivum", "T.brucei", "T.cruzi", "T.gondii", "T.thermophila","T.vaginalis","X.laevis"]
micro_name = ["E.cuniculi","P.locustae","S.lophii","V.necatrix"]

###############
## Reference ##
###############

# Run reference path
ref_id = ["ecun"]
ref_name = ["E.cuniculi"]
path = "/home/twlodarski/Projects/ribosome_tunnels/eukaryota/"
ref_dict = get_dictionary(ref_id,ref_name)

# Grid for calcualtions
grid_resolution = 1.0  # Grid resolution in Angstroms
# getting grid for each trajectory and finding the common one with the use of get_grid.py script
# When grid is ready we can calculate occupancy map for each trajectory using -> occupancy_map.py script
grid_file = "/home/twlodarski/Projects/ribosome_tunnels/grid_data.dat"
grid_dim_min,grid_dim_max = read_grid(grid_file)

# Read occupancy max path
read_occupancy_path(path,ref_dict,ref_id)

# Calculations of the spline function
calculate_spline(ref_dict,ref_id,s=60,NN=25,short=False)

# write occ path avg, its interpolation and spline into PDB files:
write_paths(path,ref_dict,ref_id,NN=25)


# Reference tunnel path
ref_tunnel_path = np.array(ref_dict[ref_id[0]]['spline path avg']['60'])

#######################################
# Selection of the group of organismsFit #
#######################################
theme = "bacteria"

# Env setup
if theme == "bacteria":
    path = "/home/twlodarski/Projects/ribosome_tunnels/bacteria/"
elif theme == "eukaryota":
    path = "/home/twlodarski/Projects/ribosome_tunnels/eukaryota/"
elif theme == "eukaryota_core":
    path = "/home/twlodarski/Projects/ribosome_tunnels/eukaryota/"  
elif theme == "micro":
    path = "/home/twlodarski/Projects/ribosome_tunnels/eukaryota/"          
elif theme == "archaea":
    path = "/home/twlodarski/Projects/ribosome_tunnels/archaea/"


##############################
# Generating main dictionary #
##############################

bac_dict = get_dictionary(bac_id,bac_name)
arc_dict = get_dictionary(arc_id,arc_name)
euk_dict = get_dictionary(euk_id,euk_name)
euk_core_dict = get_dictionary(euk_core_id,euk_core_name)
micro_dict = get_dictionary(micro_id,micro_name)

if theme == "bacteria":
    org_dict = bac_dict
    org_id = bac_id
elif theme == "eukaryota":
    org_dict = euk_dict
    org_id = euk_id
elif theme == "eukaryota_core":
    org_dict = euk_core_dict
    org_id = euk_core_id   
elif theme == "micro":
    org_dict = micro_dict
    org_id = micro_id    
elif theme == "archaea":
    org_dict = arc_dict
    org_id = arc_id


# Read occupancy max path
read_occupancy_path(path,org_dict,org_id,short=False)

# Interpolation with the spline functionget_rings(org_dict,dot_threshold):
if org_id == bac_id or org_id == chloro_id or org_id == mito_id: nn = 20
else: nn = 25

calculate_spline(org_dict,org_id,s=60,NN=nn,short=False)

# write occ path avg, its interpolation and spline into PDB files:
write_paths(path,org_dict,org_id,NN=nn)

#### 3D PLOTs
paths_plot3D(org_dict,name="ecoli",item="60")

#################################################################################
# Getting edges and interior of the occupancy maps at certain map threshold     #
# Generating average occupancy map and corresponding edges.                     #
#################################################################################

get_edge_from_occupancy(path,grid_dim_min,grid_resolution,org_dict,org_id,0.0001,short=False)

####################################################################
# Getting the structure of the ribosome tunnel that is interacting #
# with the growing nascent chain along with its physico-chemical   #
# properties                                                       #
####################################################################

get_tunnel(path,org_dict,org_id,dist_threshold=5.0,short=False)

output_name = theme+"_tunnel_composition_new.svg"
plot_tunnel_prop(org_dict,org_id,output_name,short=False)

write_tunnel_prop(org_dict,org_id,theme,short=False)

####################################################################
# Calculating properties along the average tunnel                  #
# We caclulate properties along the spline and perpendicular to it #
####################################################################

# Define a small threshold of the dot product use to finding perpendicular vectors to the spline
dot_threshold = 0.5 # Adjust this value based on your specific requirement

# Getting rings
#get_rings(path,org_dict,ref_dict,org_id,dot_threshold,ref_org="ecun")
get_rings(path,org_dict,org_id,dot_threshold,ref_tunnel_path,short=False)


# Getting discs and cross-sections
get_discs(path,org_dict,org_id,dot_threshold,ref_tunnel_path,short=False)
#get_discs_without_clusters(path,org_dict,org_id,dot_threshold,ref_org="linn")

tunnel_length = org_dict['lmaj']['path length']['60'][:-1]
np.savetxt("tunnel_mito.dat",tunnel_length) 

# Plot properties in the function of the distance from the PTC:
output_name = theme+"_tunnel_composition"
plot_prop_dist(org_dict,org_id,output_name,end="NaN",short=False)

# Getting interpolation of the cross-sections on the common path
#interp_M[theme],new_length[theme] = path_interp(org_dict,org_id)

# Plotting avg cross-section end = 91.71631884486219
output_name = theme+"_avg_cross_sections_new.svg"
plot_avg_cross_section(org_dict,org_id,output_name,end="NaN",short=False)


item="20"
tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
path_matrix = np.zeros([len(org_id),tunnel_length.size,])
for i,org in enumerate(org_id):
    path_matrix[i] = org_dict[org]['cross section'][item]

path_avg = np.nanmean(path_matrix,0)
plt.plot(tunnel_length,path_avg,'o-')



plt.axvline(32)
plt.axvline(50)
plt.axvline(67)
#plt.axvline(70)
plt.show()

#
write_tunnel_path(org_dict,org_id,theme,short=False)
write_tunnel_asph(org_dict,org_id,theme)
# Bacteria features: 29, 45, 60 and 70.
# Eukaryota features: 29, 47, 65
# Archaea feature: 32, 50, 67

# Plotting cross-sections 
output_name = theme+"_cross_sections_new.svg"
plot_cross_section(org_dict,org_id,output_name,end="NaN")


# PCA
end = 91.71631884486219
tunnel_length = org_dict[org_id[0]]['path length']['60'][:-1]
cut = np.where(tunnel_length==end)[0][0]


for item in ["10", "20", "30", "40","60"]:

M_10 = []
start_new = 0
end_new = np.inf

for name in org_id:
    item = "10"
    # Extract and smooth
    y = np.array(org_dict[name]['cross section'][item])
    y_smooth = gaussian_filter1d(y, sigma=1) 
    # Find valid (non-NaN) range for this profile
    non_nan_indices = np.where(~np.isnan(y_smooth))[0]
    if non_nan_indices.size > 0:
        start = non_nan_indices[0]
        end = non_nan_indices[-1]
        start_new = max(start_new, start)
        end_new = min(end_new, end)


# Plotting cross-sections for each length 
output_name = theme+"_cross_sections_length_new.svg"
plot_cross_section_length(org_dict,org_id,output_name,end="NaN")

################
# Mitochondria #
################


# Calculating RMSD matrices:
rmsd_M = get_rmsd_m(org_id,interp_M)

# Ploting RMSD matrices and clustering
treshold = [40,40,58,131,220]
plot_rmsd(org_id,rmsd_M)


import os
import MDAnalysis as mda
import numpy as np
import argparse
import mrcfile
from tqdm import tqdm

# ==============================================================================
# CORE FUNCTIONS
# ==============================================================================
def map_traj(u, grid_shape, origin_min, grid_resolution):
    """Maps atom positions to a 3D grid to calculate occupancy using fast NumPy arrays."""
    occupancy_grid = np.zeros(grid_shape)
    for ts in u.trajectory:
        positions = (u.atoms.positions - origin_min) / grid_resolution
        indices = np.round(positions).astype(int)
        indices = np.clip(indices, 0, grid_shape - 1)
        unique_indices = np.unique(indices, axis=0)
        occupancy_grid[tuple(unique_indices.T)] += 1
    occupancy_grid /= len(u.trajectory)
    return occupancy_grid

def write_map(filepath, occupancy_grid, grid_resolution, origin_coords):
    """Safely writes the MRC file using explicit arguments."""
    with mrcfile.new(filepath, overwrite=True) as mrc:
        mrc.set_data(occupancy_grid.T.astype(np.float32))
        mrc.voxel_size = grid_resolution
        mrc.header.origin.x = float(origin_coords[0])
        mrc.header.origin.y = float(origin_coords[1])
        mrc.header.origin.z = float(origin_coords[2])

# ==============================================================================
# CONFIGURATION & SETUP
# ==============================================================================
parser = argparse.ArgumentParser(description="Convert MD trajectories into 3D occupancy maps.")
# Use 'theme' instead of 'organism' to perfectly match your notebook variables
parser.add_argument('-o','--theme', type=str, required=True,
                    choices=['bacteria', 'eukaryota', 'archaea'])
args = parser.parse_args()

BASE_PATH = "./data_new/" # Pointing to the newly organized directory
GRID_RES = 1.0  

# 1. Dictionaries (Synced perfectly with your notebook)
bacteria = {
    'abau': 'A.baumannii', 'bbur': 'B.burgdorferi', 'bsub': 'B.subtilis',
    'drad': 'D.radiodurans', 'ecoli': 'E.coli', 'efae': 'E.faecalis',
    'fjoh': 'F.johnsoniae', 'linn': 'L.innocua', 'llac': 'L.lactis',
    'lmon': 'L.monocytogenes', 'mpne': 'M.pneumoniae', 'msme': 'M.smegmatis',
    'mtub': 'M.tuberculosis', 'paer': 'P.aeruginosa', 'pgin': 'P.gingivalis',
    'pura': 'P.urativorans', 'saur': 'S.aureus', 'sfra': 'S.fradiae',
    'tthe': 'T.thermophilus', 'vnat': 'V.natrigens'
}

archaea = {
    'hmar': 'H.marismorturi', 'mace': 'M.acetivorans', 'pfur': 'P.furiosus', 
    'saci': 'S.acidocaldarius', 'tkod': 'T.kodakarensis', 'pcal': 'P.calidifontis'
}

eukaryotes = {
    'atha': 'A.thaliana', 'calb': 'C.albicans', 'cele': 'C.elegans',
    'cpor': 'C.porcellus', 'cser': 'C.sericeus', 'cthe': 'C.thermophilia',
    'dmel': 'D.melanogaster', 'drer': 'D.rerio', 'ecun': 'E.cuniculi',
    'egra': 'E.gracilis', 'ehis': 'E.histolytica',
    'ggal': 'G.gallus', 'glam': 'G.lamblia', 'hgla': 'H.glaber',
    'hsap': 'H.sapiens', 'klac': 'K.lactis', 'ldon': 'L.donovani',
    'lmaj': 'L.major', 'mmus': 'M.musculus', 'ncra': 'N.crasa',
    'ntab': 'N.tabacum', 'ocun': 'O.cuniculus', 'pfal': 'P.falciparum',
    'ploc': 'P.locustae', 'rnor': 'R.norvegicus', 'scer': 'S.cerevisiae',
    'slop': 'S.lophii', 'slyc': 'S.lycoperscium', 'spom': 'S.pombe',
    'sscr': 'S.scrofa', 'taes': 'T.aestivum', 'tbru': 'T.brucei',
    'tcru': 'T.cruzi', 'tgon': 'T.gondii', 'tthe': 'T.thermophila',
    'tvag': 'T.vaginalis', 'vnec': 'V.necatrix', 'xlae': 'X.laevis'
}


# 2. Master Domain Configuration
DOMAINS = {
    "bacteria":       {"folder": "bacteria",  "ids": list(bacteria.keys()),       "top": "NC_FME"},
    "eukaryota":      {"folder": "eukaryota", "ids": list(eukaryotes.keys()),     "top": "NC_MET"},
    "archaea":        {"folder": "archaea",   "ids": list(archaea.keys()),        "top": "NC_MET"}
}

ITEMS = ["10", "20", "30", "40", "60"]
TRAJECTORIES = {
    "fitted.xtc": "occupancy_map.mrc",
    "ref_1/fitted.xtc": "ref_1/occupancy_map.mrc",
    "ref_2/fitted.xtc": "ref_2/occupancy_map.mrc"
}

# Load the correct grid boundary file
grid_filename = "grid_data.dat"

# Safety fallback: Check new base first, then fall back to old data folder
grid_path = os.path.join(BASE_PATH, grid_filename)

grid_data = np.loadtxt(grid_path)

grid_dim_min = {str(int(grid_data[i, 0])): grid_data[i, 1:] for i in range(0, len(grid_data), 2)}
grid_dim_max = {str(int(grid_data[i+1, 0])): grid_data[i+1, 1:] for i in range(0, len(grid_data), 2)}
# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
config = DOMAINS[args.theme]
print(f"\n--- Generating MRC Maps for Theme: {args.theme.upper()} ---")

for item in tqdm(ITEMS, desc="Lengths"):
    origin = grid_dim_min[item]
    grid_shape = np.ceil((grid_dim_max[item] - origin) / GRID_RES).astype(int)
    for org_id in config["ids"]:
        work_dir = os.path.join(BASE_PATH, config["folder"], org_id, item)
        topology_file = os.path.join(BASE_PATH, 'top/',f"{config['top']}_{item}.pdb")
        for traj_name, out_mrc in TRAJECTORIES.items():
            traj_file = os.path.join(work_dir, traj_name)
            mrc_file = os.path.join(work_dir, out_mrc)
            if os.path.exists(traj_file):
                try:
                    u = mda.Universe(topology_file, traj_file)
                    occupancy_grid = map_traj(u, grid_shape, origin, GRID_RES)
                    write_map(mrc_file, occupancy_grid, GRID_RES, origin)
                except Exception as e:
                    tqdm.write(f"❌ Error processing {traj_file}: {e}")

print("✅ Occupancy mapping complete.\n")
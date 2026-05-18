import os
import MDAnalysis as mda
import numpy as np
import argparse
from tqdm import tqdm

# ==============================================================================
# CORE FUNCTIONS
# ==============================================================================
PDB_FORMAT = "ATOM  {:5d}  {:<4s}MOL     1    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           {}\n"
            

def residue_occ_path_coord(u, output_name, grid_shape, grid_dim_min, grid_resolution):
    """Calculates the max occupancy coordinate for each residue's CA atom and writes a PDB."""
    atom_type = "C"
    # Pre-select atoms for each residue to avoid redundant selections
    residue_atoms = {r: u.select_atoms(f"resid {r} and name CA") for r in u.residues.resids}
    occupancy_grid = {r: np.zeros(grid_shape) for r in u.residues.resids}
    # Generate occupancy grid for each residue
    for ts in u.trajectory:
        for r, atoms in residue_atoms.items():
            if len(atoms) == 0: 
                continue # Safety skip if a residue lacks a CA atom
            positions = (atoms.positions - grid_dim_min) / grid_resolution
            indices = np.round(positions).astype(int)
            # Safety clip (prevents crash if an atom wiggles slightly outside global bounds)
            indices = np.clip(indices, 0, grid_shape - 1) 
            occupancy_grid[r][tuple(indices[0])] += 1
    path_coords = []
    # Safely open file using a context manager
    with open(output_name, 'w') as file:
        for r in occupancy_grid:
            occupancy_grid[r] /= len(u.trajectory)
            # Find the max occupancy value and its grid indices for each residue
            max_index = np.unravel_index(np.argmax(occupancy_grid[r]), occupancy_grid[r].shape)
            # Convert grid indices back to real-space coordinates
            max_coords = np.array(max_index) * grid_resolution + grid_dim_min
            path_coords.append(max_coords)
            file.write(PDB_FORMAT.format(r, atom_type, max_coords[0], max_coords[1], max_coords[2], atom_type))
    return path_coords

# ==============================================================================
# CONFIGURATION & SETUP
# ==============================================================================
parser = argparse.ArgumentParser(description="Extract max occupancy pathways for each residue.")
parser.add_argument('-o','--theme', type=str, required=True, choices=['bacteria', 'eukaryota', 'archaea'])
args = parser.parse_args()

BASE_PATH = "./data_new/"
GRID_RES = 1.0  

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

DOMAINS = {
    "bacteria":  {"folder": "bacteria",  "ids": list(bacteria.keys()),   "top": "NC_FME"},
    "eukaryota": {"folder": "eukaryota", "ids": list(eukaryotes.keys()), "top": "NC_MET"},
    "archaea":   {"folder": "archaea",   "ids": list(archaea.keys()),    "top": "NC_MET"}
}

ITEMS = ["10", "20", "30", "40", "60"]

# Map the input XTC to the desired output PDB
TRAJECTORIES = {
    "fitted.xtc": "occ_max_res.pdb",
    "ref_1/fitted.xtc": "ref_1/occ_max_res.pdb",
    "ref_2/fitted.xtc": "ref_2/occ_max_res.pdb"
}

grid_data = np.loadtxt(os.path.join(BASE_PATH, "grid_data.dat"))
grid_dim_min = {str(int(grid_data[i, 0])): grid_data[i, 1:] for i in range(0, len(grid_data), 2)}
grid_dim_max = {str(int(grid_data[i+1, 0])): grid_data[i+1, 1:] for i in range(0, len(grid_data), 2)}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
config = DOMAINS[args.theme]
print(f"\n--- Generating Max Occupancy Pathways for Theme: {args.theme.upper()} ---")

for item in tqdm(ITEMS, desc="Lengths"):
    origin = grid_dim_min[item]
    grid_shape = np.ceil((grid_dim_max[item] - origin) / GRID_RES).astype(int)
    for org_id in config["ids"]:
        work_dir = os.path.join(BASE_PATH, config["folder"], org_id, item)
        topology_file = os.path.join(BASE_PATH, 'top/',f"{config['top']}_{item}.pdb")
        if not os.path.exists(topology_file):
            continue 
        for traj_name, out_pdb in TRAJECTORIES.items():
            traj_file = os.path.join(work_dir, traj_name)
            output_pdb = os.path.join(work_dir, out_pdb)
            if os.path.exists(traj_file):
                try:
                    u = mda.Universe(topology_file, traj_file)
                    residue_occ_path_coord(u, output_pdb, grid_shape, origin, GRID_RES)
                except Exception as e:
                    tqdm.write(f"❌ Error processing {traj_file}: {e}")

print("✅ Occupancy pathway generation complete.\n")
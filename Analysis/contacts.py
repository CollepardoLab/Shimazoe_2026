import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import seaborn as sns


dump_file = "equil_T0.dump"  
contact_cutoff = 7.5
csv_output = "linker_histone_dna_segment_contacts.csv"

linker_segments = {
        'NTD':(0,38),
        'Globular':(38,108),
        'CTD':(108,213)
        }


DNA_SEGMENT_SIZE = 195 * 3  # 585 beads
SINGLE_BLOCK_TOTAL = 21264 
BLOCK_PROTEIN = 14244
BLOCK_DNA = 7020

dna_region_boundaries = [
    ("linker", 0, 72),
    ("side", 72, 138),
    ("back", 138, 204),
    ("side", 204, 270),
    ("dyad", 270, 315),
    ("side", 315, 381),
    ("back", 381, 447),
    ("side", 447, 513),
    ("linker", 513, 585),
]


def get_linker_segment(res_idx):
    for name, (start, end) in linker_segments.items():
        if start <= res_idx < end:
            return name
    return None

def get_dna_segment_index(global_dna_index):
    return global_dna_index // DNA_SEGMENT_SIZE

def get_dna_subregion(relative_index):
    for region, start, end in dna_region_boundaries:
        if start <= relative_index < end:
            return region
    return None

contacts = []

with open(dump_file, 'r') as f:
    frame_idx = -1
    atom_data = []

    while True:
        line = f.readline()
        if not line:
            break

        if line.startswith("ITEM: TIMESTEP"):
            frame_idx += 1
            atom_data = []
            continue

        if line.startswith("ITEM: NUMBER OF ATOMS"):
            n_atoms = int(f.readline().strip())
            continue

        if line.startswith("ITEM: BOX BOUNDS"):
            for _ in range(3):
                f.readline()
            continue

        if line.startswith("ITEM: ATOMS"):
            headers = line.strip().split()[2:]
            id_col = headers.index('id')
            x_col = headers.index('x')
            y_col = headers.index('y')
            z_col = headers.index('z')
            mol_col = headers.index('mol')
            type_col = headers.index('type')

            for _ in range(n_atoms):
                parts = f.readline().split()
                aid = int(parts[id_col])
                mol = int(parts[mol_col])
                atype = int(parts[type_col])
                pos = np.array([float(parts[x_col]), float(parts[y_col]), float(parts[z_col])])
                atom_data.append((aid, mol, atype, pos))

            linker_atoms_by_segment = {seg: [] for seg in linker_segments}
            for aid, mol, atype, pos in atom_data:
                if atype not in [41, 42] and mol % 2 == 1:  # linker histone
                    res_idx = (aid - 1) % 213
                    segment = get_linker_segment(res_idx)
                    if segment:
                        linker_atoms_by_segment[segment].append((aid, pos))

            dna_atoms = [(aid, pos) for aid, mol, atype, pos in atom_data if atype in [41, 42]]
            if not dna_atoms:
                continue

            dna_coords = np.array([pos for _, pos in dna_atoms])
            dna_ids = [aid for aid, _ in dna_atoms]
            tree = cKDTree(dna_coords)

            for segment_name, linker_atoms in linker_atoms_by_segment.items():
                for lid, lpos in linker_atoms:
                    indices = tree.query_ball_point(lpos, r=contact_cutoff)
                    for idx in indices:
                        did = dna_ids[idx]
                        dist = np.linalg.norm(lpos - dna_coords[idx])
                        if n_atoms == SINGLE_BLOCK_TOTAL:
                            relative_idx = idx
                        else:
                            dna_base = 0
                            for i in range(9):
                                if idx < (i + 1) * BLOCK_DNA:
                                    break
                                dna_base += BLOCK_DNA
                            relative_idx = idx - dna_base

                        dna_segment_idx = get_dna_segment_index(relative_idx)
                        dna_subregion = get_dna_subregion(relative_idx % DNA_SEGMENT_SIZE)
                        contacts.append((frame_idx, lid, segment_name, did, dna_segment_idx, dna_subregion, dist))


df = pd.DataFrame(contacts, columns=["frame", "linker_atom_id", "linker_segment", "dna_atom_id", "dna_segment", "dna_subregion", "distance"])
df.to_csv(csv_output, index=False)
print(f"Saved {len(df)} contacts to {csv_output}")



#df_filtered = df[df["frame"] >= 300]

interaction_counts = df.groupby(["linker_segment", "dna_subregion"]).size().unstack(fill_value=0)

linker_order = ["NTD", "Globular", "CTD"]
dna_order = ["linker", "side", "dyad", "back"]
interaction_counts = interaction_counts.reindex(index=linker_order, columns=dna_order, fill_value=0)

total_contacts = interaction_counts.values.sum()
interaction_percentages = (interaction_counts / total_contacts) * 100

plt.figure(figsize=(10, 8))
sns.heatmap(interaction_percentages, cmap="Reds", annot_kws={"size": 10}, cbar_kws={'label': 'Interaction %'})
plt.title("Linker Histone – DNA Region Contact Frequency")
plt.xlabel("DNA Region")
plt.ylabel("Linker Histone Segment")
plt.tight_layout()
plt.savefig("linker_histone_dna_contact_heatmap.svg", dpi=300)
plt.show()


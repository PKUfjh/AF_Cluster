import Bio.PDB
from tmtools import tm_align
from biopandas.pdb import PandasPdb
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import os
from matplotlib import pyplot as plt

# TODO: export this to a general file, not in evaluation scripts; and import here
# knowledge base
PROT_ONE2THREE = {
    "R": "ARG",
    "H": "HIS",
    "K": "LYS",
    "D": "ASP",
    "E": "GLU",
    "S": "SER",
    "T": "THR",
    "N": "ASN",
    "Q": "GLN",
    "C": "CYS",
    "G": "GLY",
    "P": "PRO",
    "A": "ALA",
    "V": "VAL",
    "I": "ILE",
    "L": "LEU",
    "M": "MET",
    "F": "PHE",
    "Y": "TYR",
    "W": "TRP",
    "X": "UNK",
}
PROT_THREE2ONE = {v: k for k, v in PROT_ONE2THREE.items()}

# Function to extract the sequence and coordinates from a PDB file
def extract_seq_and_coords_from_model(pandas_pdb: PandasPdb, model_id: int = 1) -> Tuple[str, np.ndarray]:
    """Reads the CA coords and FASTA sequence for a chain, excluding residues marked as 'UNK'."""
    model = pandas_pdb.get_model(model_id)

    # Create a mask for CA atoms where the residue is not 'UNK'
    ca_mask = (model.df["ATOM"].atom_name == "CA") & (
        model.df["ATOM"].residue_name != "UNK")

    # Apply the mask to extract coordinates and residue names
    ca_coords = model.df["ATOM"][ca_mask][["x_coord",
                                           "y_coord", "z_coord"]].values.astype(np.float64)
    three_seq = model.df["ATOM"][ca_mask].residue_name.tolist()

    # Convert the three-letter codes to one-letter codes, excluding 'UNK'
    fasta_seq = ''.join([PROT_THREE2ONE.get(aa, 'X')
                        for aa in three_seq if aa != "UNK"])

    return fasta_seq, ca_coords

# Function calculate the TM-score between two structures using tmtools package
# Input: 2 pdb files
# Output: TM-score
def calc_tmscore(pdb1, pdb2):
    # Read the PDB files
    ppdb1 = PandasPdb().read_pdb(pdb1)
    ppdb2 = PandasPdb().read_pdb(pdb2)

    # Extract the sequence and coordinates for each model
    seq1, coords1 = extract_seq_and_coords_from_model(ppdb1)
    seq2, coords2 = extract_seq_and_coords_from_model(ppdb2)

    # Calculate the TM-score
    tmscore = tm_align(coords1, coords2, seq1, seq2)

    return tmscore.tm_norm_chain2

# Example
pdb1 = '/lustre/home/2001110077/software/AF_Cluster/KaiB_inference/KaiB_000.pdb'
pdb2 = '/lustre/home/2001110077/software/AF_Cluster/KaiB_inputs/2QKE_A.pdb'
tmscore = calc_tmscore(pdb1, pdb2)
print(tmscore)

# Function to calculate the TM-score between between the pdb files in one directory w.r.t two reference pdb files
# Output is a dictionary with the TM-scores, where the key is the pdb file name
# Output a 2D scatter plot of the TM-scores, save the plot as a png file, not show the plot
# Output a txt file which contains the TM-scores
# The dic key and plot axis are the reference pdb filenames (only the filename, not the full path)
# The dict keys are sorted in the txt file
def calc_tmscore_all(pdb_dir, ref1, ref2):
    # Get all the pdb files in the directory
    pdb_files = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f.endswith('.pdb')]

    # Calculate the TM-score for each pdb file
    tmscores = {}
    for pdb in pdb_files:
        tmscore1 = calc_tmscore(pdb, ref1)
        tmscore2 = calc_tmscore(pdb, ref2)
        tmscores[os.path.basename(pdb)] = [tmscore1, tmscore2]

    # Save the TM-scores to a txt file
    with open(os.path.join(pdb_dir, 'tmscores.txt'), 'w') as f:
        for pdb, scores in sorted(tmscores.items()):
            f.write(f"{pdb}\t{scores[0]}\t{scores[1]}\n")

    # Create a 2D scatter plot
    tmscores = np.array(list(tmscores.values()))

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(tmscores[:, 0], tmscores[:, 1], c=tmscores[:, 0], cmap='viridis', alpha=0.7, edgecolors='w', s=50)

    # Add labels and title, strip the .pdb postfix
    ax.set_xlabel(f"TM-score to {os.path.basename(ref1)[:-4]}", fontsize=12)
    ax.set_ylabel(f"TM-score to {os.path.basename(ref2)[:-4]}", fontsize=12)
    ax.set_title("2D Scatter Plot of TM-scores", fontsize=14)

    # Customize grid and axes
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_xlim([min(tmscores[:, 0]) - 0.1, max(tmscores[:, 0]) + 0.1])
    ax.set_ylim([min(tmscores[:, 1]) - 0.1, max(tmscores[:, 1]) + 0.1])

    # Add colorbar (if needed, to show how colors correspond to values)
    cbar = plt.colorbar(scatter)
    cbar.set_label('TM-score (Color scale)', rotation=270, labelpad=15)

    # Save the figure
    plt.tight_layout()
    plt.savefig(os.path.join(pdb_dir, 'tmscores_refined.png'),dpi=300)
    plt.close()


ref_pdb_1 = "/lustre/home/2001110077/software/AF_Cluster/KaiB_inputs/2QKE_A.pdb"
ref_pdb_2 = "/lustre/home/2001110077/software/AF_Cluster/KaiB_inputs/5JYT_A.pdb"
calc_tmscore_all('/lustre/home/2001110077/software/AF_Cluster/KaiB_inference', ref_pdb_1, ref_pdb_2)
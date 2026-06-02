from Bio.PDB import PDBParser
from collections import Counter
import pandas as pd
import sys

if len(sys.argv) < 2:
    print("Usage: python analyzer.py protein.pdb")
    sys.exit()

pdb_file = sys.argv[1]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)

chains = 0
residues = 0
atoms = 0

chain_ids = []
aa_counter = Counter()

for model in structure:
    for chain in model:

        chains += 1
        chain_ids.append(chain.id)

        for residue in chain:

            if residue.id[0] == " ":
                residues += 1
                aa_counter[residue.resname] += 1

            for atom in residue:
                atoms += 1

print("=" * 50)
print("PROTEIN STRUCTURE ANALYSIS REPORT")
print("=" * 50)

print(f"Chains          : {chains}")
print(f"Chain IDs       : {', '.join(chain_ids)}")
print(f"Residues        : {residues}")
print(f"Atoms           : {atoms}")
print(f"AA Types        : {len(aa_counter)}")

most_common = aa_counter.most_common(1)[0]
print(f"Most Common AA  : {most_common[0]} ({most_common[1]})")

print("\nAmino Acid Composition")
print("-" * 30)

for aa, count in aa_counter.items():
    print(f"{aa}: {count}")

df = pd.DataFrame(
    aa_counter.items(),
    columns=["Amino_Acid", "Count"]
)

df.to_csv("amino_acid_composition.csv", index=False)

with open("report.txt", "w") as f:
    f.write("PROTEIN STRUCTURE ANALYSIS REPORT\n")
    f.write(f"Chains: {chains}\n")
    f.write(f"Chain IDs: {', '.join(chain_ids)}\n")
    f.write(f"Residues: {residues}\n")
    f.write(f"Atoms: {atoms}\n")

print("\nCSV saved as amino_acid_composition.csv")
print("Report saved as report.txt")

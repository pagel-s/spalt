#!/usr/bin/env python3
import os
import re
import csv
import argparse
from rdkit import Chem

def extract_number_from_filename(filename):
    """Extracts the trailing number before '_aligned.sdf'"""
    match = re.search(r"(\d+)_aligned\.sdf$", filename)
    return int(match.group(1)) if match else float("inf")

def parse_sdf_file(filepath):
    suppl = Chem.SDMolSupplier(filepath, removeHs=False)
    best = None
    best_fit = float("-inf")

    for mol in suppl:
        if mol is None:
            continue
        try:
            fit = float(mol.GetProp("AlignmentFitness"))
        except Exception:
            continue  # skip molecules without property

        if fit > best_fit:
            best_fit = fit
            best = mol

    if best is None:
        return None

    row = {
        "File": os.path.basename(filepath),
        "AlignmentFitness": best.GetProp("AlignmentFitness"),
        "AlignmentRMSE": best.GetProp("AlignmentRMSE") if best.HasProp("AlignmentRMSE") else None,
        "SMILES": best.GetProp("SMILES") if best.HasProp("SMILES") else Chem.MolToSmiles(best),
    }
    return row

def main():
    parser = argparse.ArgumentParser(description="Extract best molecule per SDF file (highest AlignmentFitness)")
    parser.add_argument("folder", help="Folder containing *_aligned.sdf files")
    args = parser.parse_args()

    folder = args.folder
    rows = []

    # Only files ending in _aligned.sdf
    files = [f for f in os.listdir(folder) if f.endswith("_aligned.sdf")]
    files.sort(key=extract_number_from_filename)

    for filename in files:
        filepath = os.path.join(folder, filename)
        row = parse_sdf_file(filepath)
        if row:
            rows.append(row)

    if rows:
        out_file = os.path.join(folder, "summary.csv")
        with open(out_file, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=["File", "AlignmentFitness", "AlignmentRMSE", "SMILES"])
            writer.writeheader()
            writer.writerows(rows)
        print(f"[+] Summary saved to {out_file}")
    else:
        print("No valid molecules found.")

if __name__ == "__main__":
    main()

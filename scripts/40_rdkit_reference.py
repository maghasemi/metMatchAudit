"""Reference InChIKeys computed with RDKit, for the v0.3 feasibility test.

Takes the five rows discussed during the manual review plus a random
sample of 50 LEOCC rows and 10 Broad rows that have a SMILES, and
computes InChI and InChIKey from the SMILES with RDKit. The R script
scripts/40_v03_feasibility.R compares these with what PubChem returns
for the same SMILES.

The CXSMILES extension some Broad rows carry (the part after a space,
such as "|&1:1,4|") is dropped before parsing, as the R script does.

Run from the project root: python3 scripts/40_rdkit_reference.py
Needs: pandas, rdkit (pip install rdkit)
"""
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")
base = "data-raw/MetLinkR/inst/extdata/Example_output/"
leocc = pd.read_csv(base + "LEOCC_Metabolon_Annotations_metLinkR.csv", dtype=str, encoding="latin-1")
broad = pd.read_csv(base + "Broad_2022Aug_annotations_metLinkR.csv", dtype=str, encoding="latin-1")

leocc["row_id"] = ["LEOCC_Metabolon:%d" % (i + 1) for i in range(len(leocc))]
broad["row_id"] = ["Broad_2022Aug:%d" % (i + 1) for i in range(len(broad))]
# Keep only the columns we need before renaming, since the Broad file
# already has a column called "name".
leocc = leocc[["row_id", "CHEMICAL_NAME", "SMILES", "INCHIKEY"]].set_axis(
    ["row_id", "name", "smiles", "stored_key"], axis=1)
broad = broad[["row_id", "Metabolite", "smiles", "inChIKey"]].set_axis(
    ["row_id", "name", "smiles", "stored_key"], axis=1)

named = ["dehydroepiandrosterone sulfate (DHEA-S)", "glycochenodeoxycholate",
         "glycolithocholate", "taurocholate"]
picked = [leocc[leocc["name"] == n].assign(set="named") for n in named]
picked.append(broad[broad["name"] == "Quinic acid"].assign(set="named"))

has_smiles_l = leocc[leocc["smiles"].notna() & ~leocc["name"].isin(named)]
has_smiles_b = broad[broad["smiles"].notna() & (broad["name"] != "Quinic acid")]
picked.append(has_smiles_l.sample(50, random_state=20261001).assign(set="sample"))
picked.append(has_smiles_b.sample(10, random_state=20261001).assign(set="sample"))

ref = pd.concat(picked)[["set", "row_id", "name", "smiles", "stored_key"]].copy()

def compute(smiles):
    s = smiles.split(" ")[0]
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return pd.Series({"rdkit_inchi": None, "rdkit_key": None})
    return pd.Series({"rdkit_inchi": inchi.MolToInchi(mol), "rdkit_key": inchi.MolToInchiKey(mol)})

ref = pd.concat([ref.reset_index(drop=True), ref["smiles"].reset_index(drop=True).apply(compute)], axis=1)
ref.to_csv("scripts/40_rdkit_reference.csv", index=False)
print(ref["set"].value_counts().to_dict(), "rows;", ref["rdkit_key"].isna().sum(), "failed to parse")

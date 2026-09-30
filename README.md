# metMatchAudit

Check whether metabolites matched across datasets are really the same compound.

Tools such as [MetLinkR](https://github.com/ncats/MetLinkR) link metabolomics datasets by giving every row a standardized name. Rows that end up with the same name are treated as one metabolite. metMatchAudit runs after that step. For every pair of rows that share a name, it looks at the structural evidence and reports how the two actually relate, with the rule that decided it.

Status: v0.3. Tested on the MetLinkR example data and evaluated against a blind manual review (see [Evaluation](#evaluation)). v0.2 added lipid matches across vendors; v0.3 adds stereo comparison from full structures and a check of each row's structure against its own InChIKey.

## Why

A shared name does not always mean a shared structure. The MetLinkR authors compared their output with manual curation and found groupings of compounds reported at different structural resolution, of stereoisomers, and in a few cases of compounds that are not even isomers (Patt et al., J. Proteome Res. 2025). Their advice is to check the output by hand.

Checking by hand is what I did first, on the MetLinkR example output. Of 53 groups whose rows carried different InChIKeys:

| what the rows turned out to be | groups |
|---|---|
| one annotation broader than the other (for example SM 36:1;O2 against SM d18:1/18:0) | 27 |
| different compounds merged under one name | 13 |
| an error in a source file | 10 |
| stereoisomers | 2 |
| unclear | 1 |

Some of the wrong merges matter biologically: omega-3 and omega-6 docosapentaenoic acid, the endocannabinoid 2-arachidonoylglycerol and its isomerisation product 1-arachidonoylglycerol, 1,2- and 1,3-diacylglycerols, and asparagine with acetylspermidine. Running the audit on the rest of the example data showed internal standards merged with the metabolites they stand in for (inosine-15N4, glycocholate-d4, thymine-d4).

The source errors were a surprise. Five steroid structures in the Metabolon example file have inverted stereocentres in the ring system, so the file's own InChIKey and SMILES describe a different stereoisomer than the name (dehydroepiandrosterone sulfate, glycochenodeoxycholate, glycolithocholate, taurocholate, taurodeoxycholate). Two rows in the Broad file carry the name of one bile acid and the structure of another.

## What it reports

For each pair of rows in a group, one of:

| relationship | meaning |
|---|---|
| `exact` | same compound, or the same annotation at the same level |
| `broader_narrower` | same compound family, one row reported with less detail |
| `stereo_variant` | same connectivity, different stereochemistry, stated in a name |
| `distinct` | different compounds |
| `unresolved` | the evidence does not settle it |

When the evidence is not enough, the answer is `unresolved`. The tool does not guess. Every result carries a rule ID and a short evidence text, for example:

| row A | row B | relationship | rule | evidence |
|---|---|---|---|---|
| SM 36:1;O2 (Broad) | stearoyl sphingomyelin (d18:1/18:0) (Metabolon) | broader_narrower | R05 | SM 36:1;O2: level species vs structure |
| 2-arachidonoylglycerol (20:4) | 1-arachidonylglycerol (20:4) | distinct | R12 | different connectivity |
| Inosine | Inosine-15N4 [iSTD] | distinct | R02 | one row is an internal standard |

## How it decides

Four kinds of evidence, none of them typed in by hand:

- the InChIKey in the source file, split into its connectivity, stereo and protonation parts
- fixed text checks on the name: group names such as "(2 or 3)-decenoate", tentative marks, isomer labels, internal standards, stated stereochemistry
- for lipids, a profile at species, molecular species and structure level. Rows with an InChIKey are placed in the [SwissLipids](https://www.swisslipids.org) hierarchy. Names in shorthand notation are parsed with [rgoslin](https://bioconductor.org/packages/rgoslin/). Both routes end in the same notation, so a Metabolon name like "stearoyl sphingomyelin (d18:1/18:0)" and a Broad name like "SM 36:1;O2" can be compared.
- optionally, the full structure (the file's InChI, or its SMILES read with [RDKit](https://www.rdkit.org) through reticulate). Two InChIKeys can only say that stereo differs; two InChIs show whether one row simply leaves some centres undefined (cholesterol with one defined centre against cholesterol with all eight) or whether they disagree. A structure is used only when its own key matches the row's stored key.

The rules run in a fixed order and the first one that applies decides. The full list, with the real case behind each rule, is in [notes/design.md](notes/design.md). Every change to it is logged with its reason in [notes/decisions.md](notes/decisions.md).

## Installation

rgoslin comes from Bioconductor, so install it first:

```r
install.packages("BiocManager")
BiocManager::install("rgoslin")
install.packages("remotes")
remotes::install_github("maghasemi/metMatchAudit")
```

The lipid evidence needs the SwissLipids table. Download `lipids.tsv` from <https://www.swisslipids.org/#/downloads> (about 690 MB, CC BY 4.0).

The structure checks are optional and need RDKit through reticulate. reticulate can set it up by itself:

```r
install.packages("reticulate")
reticulate::py_require("rdkit")
```

Without RDKit the audit runs as before and skips those checks.

## Usage

MetLinkR does not record which column holds the metabolite name or the InChIKey in each file, so you say it:

```r
library(metMatchAudit)

files <- data.frame(
  path         = c("metLinkR_output/Metabolon_metLinkR.csv",
                   "metLinkR_output/Broad_metLinkR.csv"),
  name_col     = c("CHEMICAL_NAME", "Metabolite"),
  inchikey_col = c("INCHIKEY", "inChIKey"),
  smiles_col   = c("SMILES", "smiles"),   # optional
  inchi_col    = c(NA, "inChI")           # optional
)

audit <- audit_metlinkr(files, swisslipids = "lipids.tsv")
audit
write_audit(audit, "audit_output")
```

`write_audit()` writes `pairs.csv`, `groups.csv` (the most severe relationship per group), `row_flags.csv`, `rows.csv` (all evidence per row), `lipid_matches.csv` (lipid matches across files, see below) and `run_info.txt` (package version, date, md5 of the SwissLipids file).

Loading SwissLipids takes about half a minute. For several audits, load it once with `load_swisslipids()` and pass the result.

## Evaluation

The manual review of the 53 groups was finished and committed before the tool's output on those groups was opened. Each version was compared with it group by group:

| outcome | v0.1 | v0.2 | v0.3 |
|---|---|---|---|
| agrees with the review | 28 | 29 | 34 |
| unresolved where the review has an answer | 14 | 13 | 8 |
| a different definite answer | 1 | 1 | 1 |
| source-error groups called exact | 0 of 10 | 0 of 10 | 0 of 10 |

On the 42 groups where the review gives a definite relationship, v0.3 agrees on 34 and leaves 8 unresolved. It makes no wrong definite claim on them. Of the 13 wrong merges, it finds 12. The one wrong answer is xylose, where the ring and open-chain forms have different InChIKeys; this was listed as a known limit before the first evaluation.

What changed between versions: v0.2 resolved an SM group row through its chain sum; v0.3 resolved five groups (cholesterol, cholic, deoxycholic, chenodeoxycholic and tauro-beta-muricholic acid) where one row defines only part of the other's stereo.

The 8 unresolved cases:

- 5 groups with a group-name row that has no lipid profile, or that is compared with a specific structure. The tool does not yet check which structures a group lists.
- 3 groups where the stereo really differs but no name states it (androsterone against etiocholanolone glucuronide), or where a row's structure does not match its own key (quinic acid, 2-hydroxystearic acid in the Broad file). By design these stay unresolved.

The stereo comparison was also tested on 317 structures from the example data by removing or inverting stereo marks with RDKit. Removing marks never produced a false conflict except in three symmetric molecules such as inositol, and inverting a centre never produced a false "subset". Only "subset" leads to a claim.

Everything needed to repeat this is in the repository:

```
bash scripts/00_get_data.sh              # MetLinkR example data, fixed commit
Rscript scripts/20_audit_example.R       # run the audit
Rscript scripts/21_evaluate_v0.1.R       # compare with notes/review.csv
Rscript scripts/31_v02_check_proposals.R # compare lipid matches with notes/v02_proposal_review.csv
Rscript scripts/41_v03_rdkit_feasibility.R  # RDKit checks and stereo comparison on R11 pairs
```

The manual labels, with the PubChem records behind each decision, are in [notes/review.csv](notes/review.csv).

## Lipid matches across vendors

`propose_lipid_matches()` looks for lipid rows in different files that describe the same Species, whether or not a matching step linked them. On the example data it proposes 60 matches between the Metabolon and Broad files; MetLinkR had linked 12 of them. All 60 were checked by hand against the names, and all are correct ([notes/v02_proposal_review.csv](notes/v02_proposal_review.csv)). For example:

| Metabolon | Broad | relationship |
|---|---|---|
| 1-palmitoyl-2-oleoyl-GPC (16:0/18:1) | PC 34:1 | broader_narrower |
| 1-(1-enyl-palmitoyl)-2-oleoyl-GPC (P-16:0/18:1)* | PC O-34:2 | broader_narrower |
| sphingomyelin (d18:1/20:0, d16:1/22:0)* | SM 38:1;O2 | broader_narrower |

A Broad Species often matches several Metabolon structures, such as PC 34:4 with both 14:0/20:4 and 16:1/18:3. These are marked one_to_many. A Species-level value usually covers all of them, so it should not be compared with each structure on its own.

Group-name rows like the last example get a profile only when every listed member adds up to the same chain sum and that sum matches the Species SwissLipids gives for the row's InChIKey. In the example data, 19 of 28 such rows qualify.

It works without MetLinkR too:

```r
db <- load_swisslipids("lipids.tsv")
rows <- row_evidence(my_rows)          # row_id, file, original_name, inchikey
matches <- propose_lipid_matches(lipid_evidence(rows, db), group_col = NULL)
```

## Source file checks

With RDKit, each row's own structure is turned into an InChIKey and compared with the key stored in the same row (check C04). In the Broad example file, 21 rows disagree: 2 in connectivity (xylose, L-urobilin) and 19 in the stereo block. For N-formylmethionine, the Broad InChI describes the natural L form while the stored key describes the D form.

The Metabolon file has a different kind of error that this check cannot see: five steroid rows whose SMILES and InChIKey agree with each other but describe a stereoisomer other than the one named. They were found during the manual review and are listed in [notes/findings.md](notes/findings.md).

## Known limits

- Sugars in ring and open-chain form get different InChIKeys and are called distinct.
- The stereo comparison can report a false conflict for symmetric molecules such as inositol. A conflict never leads to a claim on its own, so this can only leave a pair unresolved.
- A structure whose SMILES and key agree but which does not match its name (the Metabolon steroids above) is not detected.
- The name checks are tuned on Metabolon and Broad naming.
- Some lipids have no usable SwissLipids entry (d18:2 and d17:1 sphingoid bases in the test data). The rgoslin grammar for trivial fatty acid names is not used, because on the test data it read small organic acids as wax esters.
- Glycosylceramide pairs that differ only in headgroup naming (HexCer against GlcCer) stay unresolved.
- Pair checks cover only rows that MetLinkR grouped. Lipid matches cover all lipid rows with a profile; other unmatched rows are out of scope.
- A group-name row compared with a specific structure of the same Species stays unresolved, since the listed members are not checked yet.

## Next

- Check a row's structure against its name, to catch cases like the Metabolon steroids. This needs a reliable name-to-structure source; PubChem synonyms turned out to be wrong in at least one case during the review, so the design is still open.
- Check which structures a group-name row lists, to resolve the remaining group rows.

## Repository layout

```
R/             package code
tests/         unit tests, built from real cases in the review
scripts/       analysis scripts, numbered by phase
notes/         design, decisions, findings, manual review
results/       script output
```

## References

- Patt A. et al. metLinkR: Facilitating Metaanalysis of Human Metabolomics Data through Automated Linking of Metabolite Identifiers. J. Proteome Res. 2025, 24, 2403-2407.
- Aimo L. et al. The SwissLipids knowledgebase for lipid biology. Bioinformatics 2015, 31, 2860-2866.
- Kopczynski D. et al. Goslin: A Grammar of Succinct Lipid Nomenclature. Anal. Chem. 2020, 92, 10957-10960.
- RDKit: Open-source cheminformatics. https://www.rdkit.org
- Liebisch G. et al. Update on LIPID MAPS classification, nomenclature, and shorthand notation for MS-derived lipid structures. J. Lipid Res. 2020, 61, 1539-1555.

## License

MIT. SwissLipids data is not included and keeps its own license (CC BY 4.0).

Amin Ghasemi
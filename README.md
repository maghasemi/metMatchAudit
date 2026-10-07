# metMatchAudit

Check whether metabolites matched across datasets are really the same compound.

Tools such as [MetLinkR](https://github.com/ncats/MetLinkR) link metabolomics datasets by giving every row a standardized name. Rows that end up with the same name are treated as one metabolite. metMatchAudit runs after that step. For every pair of rows that share a name, it looks at the structural evidence and reports how the two actually relate, with the rule that decided it. It also proposes lipid matches across vendors that the matching step missed.

Status: v0.7. Developed on the MetLinkR example data and evaluated against a blind manual review, then tested blind on five COMETS datasets with expert mappings and in two further blind tests on files it had not seen (see [Evaluation](#evaluation), [External test on COMETS](#external-test-on-comets), [Second blind test](#second-blind-test) and [Third blind test](#third-blind-test)).

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
| scyllo-inositol | myo-inositol | unresolved | R03c | identical InChIKey, but the names state different isomers |

## How it decides

Every piece of evidence comes from a source file, a fixed parsing rule or a versioned database file. Nothing is typed in by hand.

- The InChIKey in the source file, split into its connectivity, stereo and protonation parts.
- For rows without a key, a key from the row's PubChem CID or HMDB ID, through PubChem's parent compound. Every key records where it came from.
- Fixed text checks on the name: group names such as "(2 or 3)-decenoate", tentative marks, isomer labels, internal standards, stated stereochemistry, omega positions.
- For lipids, a profile at species, molecular species and structure level. Rows with an InChIKey are placed in the [SwissLipids](https://www.swisslipids.org) hierarchy, or found in the [LIPID MAPS Structure Database](https://www.lipidmaps.org/databases/lmsd) when SwissLipids has nothing. Names in shorthand notation are parsed with [rgoslin](https://bioconductor.org/packages/rgoslin/). All routes end in the same notation, so a Metabolon name like "stearoyl sphingomyelin (d18:1/18:0)" and a Broad name like "SM 36:1;O2" can be compared.
- For files in Biocrates p180 nomenclature, the vendor's own old-to-new nomenclature table, so "PC aa C36:4" becomes PC 36:4 and "C14:2" becomes CAR 14:2. These values are sums measured by flow injection, so they always count as species level.
- Optionally, the full structure (the file's InChI, or its SMILES read with [RDKit](https://www.rdkit.org) through reticulate). Two InChIKeys can only say that stereo differs; two InChIs show whether one row simply leaves some centres undefined or whether they disagree.

The rules run in a fixed order and the first one that applies decides. The full list, with the real case behind each rule, is in [notes/design.md](notes/design.md). Every change is logged with its reason in [notes/decisions.md](notes/decisions.md), including the ideas that failed.

## Installation

rgoslin comes from Bioconductor, so install it first:

```r
install.packages("BiocManager")
BiocManager::install("rgoslin")
install.packages("remotes")
remotes::install_github("maghasemi/metMatchAudit")
```

Data files, not included in the package:

- SwissLipids `lipids.tsv` from <https://www.swisslipids.org/#/downloads> (about 690 MB, CC BY 4.0). Needed.
- LMSD standard SDF from <https://www.lipidmaps.org/databases/lmsd/download> (about 290 MB unzipped, CC BY 4.0). Optional.

Identifier resolution needs httr2 and jsonlite and internet access on the first run; answers are cached. The structure checks need RDKit through reticulate:

```r
install.packages(c("httr2", "jsonlite", "reticulate"))
reticulate::py_require("rdkit")
```

Without these the audit still runs and skips the parts that need them.

## Usage

MetLinkR does not record which column holds the metabolite name or the identifiers in each file, so you say it:

```r
library(metMatchAudit)

files <- data.frame(
  path         = c("out/Metabolon_metLinkR.csv", "out/Broad_metLinkR.csv", "out/Biocrates_metLinkR.csv"),
  name_col     = c("CHEMICAL_NAME", "Metabolite", "metabolite_name"),
  inchikey_col = c("INCHIKEY", "inChIKey", NA),
  pubchem_col  = c("PUBCHEM", NA, NA),            # optional
  hmdb_col     = c("HMDB", "hmdbId", NA),         # optional
  vendor       = c(NA, NA, "biocrates_p180")      # optional
)

sl   <- load_swisslipids("lipids.tsv")
lmsd <- load_lmsd("structures.sdf")               # optional
audit <- audit_metlinkr(files, swisslipids = sl, lmsd = lmsd,
                        cache_file = "pubchem_cache.rds")
audit
write_audit(audit, "audit_output")
```

`write_audit()` writes `pairs.csv`, `groups.csv` (the most severe relationship per group), `row_flags.csv`, `rows.csv` (all evidence per row), `lipid_matches.csv` and `run_info.txt` (package version, date, md5 of every database file).

Loading SwissLipids takes about half a minute and LMSD about a minute. Both objects can be kept with `saveRDS()`.

## Evaluation

The manual review of the 53 groups was finished and committed before the tool's output on those groups was opened. Each version was compared with it group by group:

| outcome | v0.1 | v0.2 | v0.3 | v0.4 | v0.4.2 and v0.6 |
|---|---|---|---|---|---|
| agrees with the review | 28 | 29 | 34 | 34 | 35 |
| unresolved where the review has an answer | 14 | 13 | 8 | 8 | 8 |
| a different definite answer | 1 | 1 | 1 | 1 | 0 |
| source-error groups called exact | 0 of 10 | 0 of 10 | 0 of 10 | 0 of 10 | 0 of 10 |

The one wrong answer up to v0.4 was xylose, whose ring and open-chain forms have different InChIKeys. Since v0.4.2, two rows with the same name and conflicting keys stay unresolved, which is what the review says for that group.

The 8 unresolved cases are 5 groups with a group-name row the tool cannot check yet, and 3 where stereo really differs but no name states it, or where a row's structure does not match its own key. By design these stay unresolved.

Lipid matches on the example data: 110 proposals between the Metabolon and Broad files, all checked by hand and all correct ([notes/v02_proposal_review.csv](notes/v02_proposal_review.csv), [notes/v06_dev_lipid_review.csv](notes/v06_dev_lipid_review.csv)).

## External test on COMETS

To see how the tool does on data it was not built on, I ran MetLinkR on five COMETS datasets and compared the result with the COMETS expert mappings from the MetLinkR supplement. The criteria were written and committed before the run, and the expert labels were not opened until the tool's output was final. Our MetLinkR run found 183 of 214 multi-file expert groups whole (85.5%; the paper reports 85.3%).

Blind result, v0.4.1:

| criterion | result | |
|---|---|---|
| C1: definite claims (exact, distinct) agree with the expert, at least 90% | 121 of 125 (96.8%) | met |
| C2: with the audit, precision and F1 not lower than MetLinkR alone | F1 0.893 against 0.895 | not met |

The audit's definite claims were mostly right, but it made few of them. 259 of 401 pairs stayed unresolved, mostly because three of the five files carry names only. Its two `distinct` calls were both wrong, and it separated none of MetLinkR's real wrong merges.

Reading the 42 pairs the expert does not group, after the result was sealed: 7 are different compounds, 4 stereoisomers, 17 a broader against a narrower annotation or a class total, 9 the same compound measured on a different platform that the expert keeps apart, and 5 look like inconsistencies in the expert labels (Nightingale alanine is grouped with beta-alanine, for example). The expert answers "same measurement"; this tool answers "same compound". Details are in [notes/findings.md](notes/findings.md) and [notes/comets_review.csv](notes/comets_review.csv).

Changes made after that test, and their results on COMETS. These numbers are after adaptation, not blind:

- v0.4.2: identical keys are not enough for `exact` when the names state different isomers (R03c), and the same name with conflicting keys stays unresolved (R12n). The four wrong claims of the blind run became unresolved. On the example data, R03c also caught a wrong `exact` call no review had seen, (R)- against (S)-3-hydroxybutyrylcarnitine.
- v0.6: LMSD as a second lipid source and the Biocrates p180 table. 23 COMETS pairs moved from unresolved to a definite answer, 22 of them across files. Two became `distinct`, valerylcarnitine against glutarylcarnitine and C14:1 against C14:2, both MetLinkR merges of different compounds; one of the two is a pair across files. Across files, definite claims now agree with the expert in 122 of 122, precision goes from 0.895 to 0.897 and F1 from 0.895 to 0.896 (recounted by `scripts/97_comets_recount.R`; an earlier version of this README gave 123, 0.900 and 0.897, worked out by hand instead of counted). 36 new lipid matches across files; the expert groups 18 of them, and the other 18 are Biocrates sums matched to single structures, which the expert links one to one or not at all.

## Second blind test

v0.6.0 ran on two files it had not seen. One is my own Metabolon annotation of a mouse dataset (494 rows, private; only counts are published). The other is the annotation from the aging mouse brain atlas of Ding et al. 2021 (1547 rows, MS-DIAL-style names). Before the run I labeled 86 pairs in 60 randomly drawn MetLinkR groups and 40 random B lipid rows, then committed the criteria and the md5 of the label files.

| criterion | result | |
|---|---|---|
| C1: definite claims (exact, distinct) agree with the label, at least 90% | 48 of 65 (73.8%) | not met |
| C2: groups labeled A to D with a wrong definite answer, at most 5% | 12 of 58 (20.7%) | not met |
| C3: lipid proposals involving a sampled B row are in the labels, at least 90% | 1 of 1 | met, on one proposal |
| C3 recall | 1 of 5 | |

11 of the 17 wrong definite claims come from one habit of the B file. Features named "Adenosine isomer", "Carnosine isomer A" and the like carry the parent compound's identifiers, so the tool called them `exact` with the parent. Most of the rest are rows whose name and key disagree: "ACar 17:0" names a sum but carries a straight-chain key. The tool follows the key, and the labels follow the name. Lipid recall was low because the tool could not read MS-DIAL names such as "PE 38:6; PE 18:2_20:4;". 375 of 960 B lipid rows got no profile, and the rest got a structure-level profile from their key that their names do not support.

The tool held up on the earlier test sets but did not carry over to an unseen naming style. Running v0.4.1 on the same files gave 48 of 67 and 13 of 58, so the later versions did not cause the drop. Details: [notes/findings.md](notes/findings.md).

v0.7, after adaptation: a new rule R03i keeps a pair unresolved when a name marks a feature of unknown identity, either by the word "isomer" or, in MS-DIAL lipid names, by a lone letter after the last semicolon ("PC 36:2; A"). On the same files, C1 goes to 48 of 53 (90.6%) and C2 to 8 of 58. The development data and COMETS do not change. These numbers come from the files the rule was written for, so they show the fix works, not that it generalizes. Reading MS-DIAL lipid names was planned for v0.7 and dropped: only 74.9% of the names parsed, short of the 80% set beforehand.

## Third blind test

v0.7.0 ran on my Metabolon file and on two published lists from Lai et al. 2021 (Nat. Commun.; mouse feces, 533 rows, and serum, 231 rows, each with PubChem CIDs and HMDB-style names). As before, I drew 60 MetLinkR groups and 40 lipid and carnitine rows at random, labeled them, and committed the criteria and the label md5s before the run. The rule for stereo in names was also fixed in advance.

| criterion | result | |
|---|---|---|
| C1: definite claims (exact, distinct) agree with the label, at least 90% | 78 of 78 (100%) | met |
| C2: groups labeled A to D with a wrong definite answer, at most 5% | 7 of 60 (11.7%) | not met |
| C3: lipid proposals involving a sampled row are in the labels, at least 90% | 13 of 13 | met |
| C3 recall | 13 of 20 | |

Every `exact` call was right, and no lipid proposal was wrong. All seven wrong groups have one cause. Two rows name the same compound, but one lab's PubChem record has stereo and the other's does not, and rule R09 then calls the pair broader_narrower. In five of the seven, neither name says anything about stereo. The same cause appeared in the second test. The sample was mostly same-compound pairs, so it tests `exact` calls far more than `distinct` ones. Details: [notes/findings.md](notes/findings.md).

## What did not work

Two ideas were tried and set aside. Both are logged in [notes/decisions.md](notes/decisions.md).

- Giving name-only rows a key from RefMet. MetLinkR's standardized names are RefMet names, so a key from RefMet would agree with MetLinkR by construction and audit nothing.
- Checking each row against the RefMet entry of its group. The check found real problems, among them a Metabolon row for 2-hydroxycarbamazepine filed under norcodeine and two RefMet records with wrong formulas (CAR 20:4 stored as C31H54NO4). But only 17 of 40 flags were wrong assignments, short of the 80% set beforehand, and the check could not change any pair result. The code is kept on branch `shelved/refmet-check`.

## Known limits

- Rows with only a name, and no identifier the tool can use, get no evidence. This is most of the COMETS Biocrates amino acids, the LC/MS file and the Nightingale file.
- Nightingale class totals (Omega-6, total fatty acids) and other aggregate measures have no relationship of their own yet.
- Sugars in ring and open-chain form get different InChIKeys. With the same name they stay unresolved; with different names they are called distinct (N-acetylneuraminic acid in the second blind test).
- Identifiers are trusted where the name contradicts them only for stereo terms, omega positions and isomer markers. A sum name with a single-structure key ("ACar 17:0" with a straight-chain key) gets the key's answer.
- MS-DIAL-style lipid names ("PE 40:6;", several annotations separated by ";", class names such as Cer-NS) are not read. Such rows get a lipid profile only from their key.
- The stereo comparison can report a false conflict for symmetric molecules such as inositol. A conflict never leads to a claim on its own.
- A structure whose SMILES and key agree but which does not match its name (the Metabolon steroids above) is not detected.
- The name checks are tuned on Metabolon, Broad and Biocrates naming. The Biocrates table covers the p180 kit as published in 2023; analytes it does not list get no profile.
- Pair checks cover only rows the matching step grouped. Lipid matches cover all lipid rows with a profile.

## Next

- A stereo-only key difference between two names that state no stereo should not give broader_narrower. This was the only cause of wrong answers in the third test and one of the causes in the second.
- A further blind test with more pairs of different compounds, since the third test had almost none.
- Reading MS-DIAL lipid names, with a class map for names such as Cer-NS and HexCer-NDS.
- A name against key check for lipids. In the second test, B names PC 35:3 for a row whose key is PE 38:3.
- Nightingale biomarker definitions and an aggregate relationship for class totals.
- Checking which structures a group-name row lists.

## Repository layout

```
R/             package code
inst/extdata/  the Biocrates p180 nomenclature table, extracted by scripts/80
tests/         unit tests, built from real cases in the reviews
scripts/       analysis scripts, numbered by phase
notes/         design, decisions, findings, manual reviews
results/       script output
```

## Data sources and references

- Patt A. et al. metLinkR: Facilitating Metaanalysis of Human Metabolomics Data through Automated Linking of Metabolite Identifiers. J. Proteome Res. 2025, 24, 2403-2407.
- Aimo L. et al. The SwissLipids knowledgebase for lipid biology. Bioinformatics 2015, 31, 2860-2866. CC BY 4.0.
- LIPID MAPS Structure Database, <https://www.lipidmaps.org>. Sud M. et al. LMSD: LIPID MAPS structure database. Nucleic Acids Res. 2007, 35, D527-D532. CC BY 4.0.
- Liebisch G. et al. Update on LIPID MAPS classification, nomenclature, and shorthand notation for MS-derived lipid structures. J. Lipid Res. 2020, 61, 1539-1555.
- Kopczynski D. et al. Goslin: A Grammar of Succinct Lipid Nomenclature. Anal. Chem. 2020, 92, 10957-10960.
- biocrates life sciences ag. AbsoluteIDQ p180 kit, old vs. new nomenclature (v1-2023).
- PubChem, <https://pubchem.ncbi.nlm.nih.gov>, through PUG REST.
- RDKit: Open-source cheminformatics, <https://www.rdkit.org>.

## License

MIT. SwissLipids and LMSD data are not included and keep their own license (CC BY 4.0). The Biocrates nomenclature table in `inst/extdata` lists analyte names from the vendor's published document.

Amin Ghasemi

# Decisions

One entry per decision, with the reason. Add new entries at the bottom and never rewrite old ones. If a decision changes, add a new entry that says so.

## 2026-09-23

Language is R. MetLinkR, RaMP and rgoslin are all R packages, so the comparison work happens there anyway. A Python version can come after the R one is finished.

renv tracks package versions. `renv.lock` is committed and the local package library is not.

Phase 1 is exploratory. There is no package structure until the phase 1 findings are written up in `notes/findings.md`.

Phase 1a uses the example output that ships with MetLinkR, at commit `a6aa145ef8cf76fb762172c6f6dd466c0b30c101`. We use it as published and do not edit it.

These example outputs do not fully match the current manifest (`inst/extdata/HarmInputFiles.csv`). The LEOCC output has 1548 rows while the closest input file has 1546, and the DDLPS output has no manifest entry at all. They probably come from an older run. Phase 1b will rerun MetLinkR ourselves with the RaMP database version written down.

Phase 1a uses only two of the output files, LEOCC (Metabolon) and Broad_2022Aug, because they are the ones with an InChIKey column. COMETS and DDLPS have PubChem IDs but no InChIKey. They come in later, once we add a PubChem lookup.

Invalid InChIKeys are counted and left out of the grouping. They are not corrected.

The `pattern` column in `results/inchikey_groups.csv` describes how the InChIKeys differ. It is not a relationship label. A block2 difference is often stereochemistry, but that has to be confirmed in review before any row gets called a stereo variant.

Review categories, fixed before looking at any review results:

- A: same compound, differs only in charge or protonation
- B: stereoisomers
- C: positional isomers or different compounds, so the MetLinkR grouping is wrong
- D: one row is a broader or narrower annotation of the other
- E: data error in the source file
- F: unclear after checking the evidence

Review selection: every `different_connectivity` group. For same-connectivity groups, all of them if there are 40 or fewer, otherwise a random sample of 40 using seed 20260923.

Manual review is used for evaluation only. Nothing typed by hand goes into the matching logic.

## 2026-09-25

## 2026-09-25: Lipid feasibility test, criteria set before running

Question: can each lipid row get a level-aware annotation from evidence alone (SwissLipids by InChIKey, rgoslin by name), with nothing typed in by hand?

SwissLipids file: lipids.tsv, md5 b6ddd5bc3003967dce95e905e5748953, 687497215 bytes, downloaded 2026-09-25.

Scope comes from the source files' own labels. LEOCC: SUPER_PATHWAY "Lipid". Broad: superClass is one of the eight LIPID MAPS categories. Results are reported per sub-pathway or class.

Complex lipids, LEOCC SUB_PATHWAY: Phosphatidylcholine (PC), Phosphatidylethanolamine (PE), Phosphatidylinositol (PI), Phosphatidylserine (PS), Lysophospholipid, Plasmalogen, Lysoplasmalogen, Diacylglycerol, Monoacylglycerol, Sphingomyelins, Dihydrosphingomyelins, Ceramides, Dihydroceramides, Hexosylceramides (HCER), Lactosylceramides (LCER), Ceramide PEs.

Complex lipids, Broad mainClass: Glycerophosphocholines, Glycerophosphoethanolamines, Monoradylglycerols, Diradylglycerols, Triradylglycerols, Sphingomyelins, Ceramides.

Criteria:
- LEOCC complex, non-group rows: at least 80% have exactly one Species ancestor in SwissLipids.
- Broad complex rows: at least 80% parsed by rgoslin to Species level or finer.

Other lipid rows are reported without a threshold. Group-name rows are measured separately and do not count toward the criteria.

Qualitative check: for lipid rows in groups already labeled in review.csv, compare SwissLipids and rgoslin output with the manual label. One unexplained disagreement stops the plan.

InChIKey matching runs in three tiers: full key, first 25 characters (ignores protonation, since SwissLipids stores pH 7.3 forms), first 14 characters. The tier is recorded per row.

If the criteria are met, SwissLipids and rgoslin go into design.md as lipid evidence sources. If not, lipid rules stay limited to InChIKey comparison, and design.md says so.

## 2026-09-29: Lipid feasibility, first run and one change after it

First run results: both criteria met. LEOCC complex non-group rows with one Species ancestor: 114 of 124 (91.9%). Broad complex rows parsed to Species or finer: 169 of 169 (100%). No disagreement with the manual review labels among the 15 labeled groups that contain lipid rows.

The 4 rows where the hierarchy and the hit abbreviation disagreed are all glycosylceramides. The hierarchy gives HexCer or Hex2Cer, the hit gives GlcCer or LacCer, with the same chain sums. This is a headgroup resolution difference, not a conflict. For cross-dataset comparison, the Species name from the SwissLipids hierarchy is the one to use, since it is written at one consistent level.

Change made after seeing the first run: 5 rows matched SwissLipids entries with no path to a Species (three 1-monoacylglycerols whose entries have no Level, two d18:2 ceramides whose parent chain stops early). For these rows only, the script now collects all SwissLipids entries with the same first 14 InChIKey characters and accepts their Species if they all agree on exactly one. Same connectivity implies same chain composition, so they cannot belong to different Species. The result goes into new columns. The pre-registered criteria are still read from the original columns, and both numbers are reported.

SwissLipids has 1292 entries with an empty Level column, 797 of them with an InChIKey. Not investigated yet.

## 2026-09-30: Name checks tested on all example names

Each name check was run on all 2514 unique names in the five MetLinkR example output files, and every match was read by hand. Four problems were found and fixed, each now covered by a test:

- The isotope pattern matched "18O2" inside the formula of "glucuronide of C10H18O2". Isotope labels must now not follow a letter or digit.
- The isomer label pattern matched residue positions in peptide names such as "bradykinin, des-arg(9)". A space before the label is now required.
- alpha and beta are left out of the stereo check, because in names such as "alpha-ketoglutarate" they mark a position. Some real stereo statements will be missed. That leads to "unresolved", which is the safe direction.
- The COMETS file joins synonyms with " or ", as in "Sarcosine or Sarcosine (N-Methylglycine)". These were first read as group names. The text cannot tell synonyms from real alternatives (the same data also has "branched-chain, straight-chain, or cyclopropyl 12:1 fatty acid"). A new check, has_alternative_names(), reports the form only. New rule R08b makes such pairs unresolved, and these rows get no lipid profile.

## 2026-09-30: Lipid rules run before the identical-key rule

The Broad row "SM 36:1;O2" carries the InChIKey LKQLRGMMMAHREN-YJFXYUILSA-N, the same key as the LEOCC row "stearoyl sphingomyelin (d18:1/18:0)". Broad reports a Species and attached one full structure as its key. With R03 (identical keys, exact) first, this pair would be called exact, while the manual review labeled it broader/narrower (D). The lipid rules R04 to R06 now run before R03, so a pair of lipid rows is judged by the level each name reports before the keys are compared. Rule IDs stay the same; only the order changed.

Related choice in lipid_evidence(): when a row's name parses with rgoslin and SwissLipids also places its key, the profile takes its level from the name, since the name states what was measured.

## 2026-09-30: Only three rgoslin grammars are accepted

On the real data, 245 rows got a lipid profile from their name. Checked by grammar:

| grammar | rows | result |
|---|---|---|
| Shorthand2020 | 182 | all correct |
| Goslin | 12 | all correct |
| LipidMaps | 7 | all correct |
| FattyAcids | 44 | not reliable |

FattyAcids read small organic acids and amino acids as wax esters ("3-hydroxyisobutyrate" as WE 4:0;O, "2-aminobutyrate" as WE 4:0;N) and miscounted double bonds ("Pentadecenoic acid" as FA 15:0, "Eicosenoic acid" as FA 20:0). Some of its results were correct, but the good ones cannot be told apart from the bad ones, so the whole grammar is left out. parse_lipid_names() now accepts Shorthand2020, Goslin and LipidMaps only, and reports the grammar for every name.

The two assumptions in the step 3 tests were checked with real rgoslin: "SM(d36:1)" gives SM 36:1;O2 and "DG(18:1_18:2)" gives molecular species DG 18:1_18:2, both through the LipidMaps grammar.
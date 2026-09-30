# Findings

Numbers here come from the scripts, not from memory. When a number is copied in, note which results file it came from.

## Phase 1a: InChIKey check on MetLinkR example output

## v0.2 results

Package output checked against the hand review (`scripts/31_v02_check_proposals.R`): all 53 reviewed proposals are reproduced, and 6 new ones come from group-name rows that now have a lipid profile. The 6 were checked by adding up each listed member's chains; all are correct. Total: 59 proposals, 59 correct, 48 of them not linked by MetLinkR.

Group-name rows with a profile: 19 of 28, all in the Metabolon file.

Effect on the v0.1 evaluation (53 reviewed groups): agree 28 to 29, safe_miss 14 to 13, wrong unchanged at 1 (xylose). The SM 36:1 group (C) now comes out distinct through R04. Of the 13 wrong merges by MetLinkR, the tool now finds 12.

## v0.3 feasibility: RDKit and stereo comparison

Source: `scripts/41_v03_rdkit_feasibility.R`, output in `results/v03_rdkit/`.

All three criteria met on the development machine, with the same results as the earlier run elsewhere: 65 of 65 reference keys reproduced; no B pair classified as subset; Broad xylose flagged.

R11 pairs by review label: D pairs 6 subset and 1 without a usable structure (quinic acid); B pairs 1 conflict (androsterone against etiocholanolone) and 1 without a usable structure (2-hydroxystearate); E pairs 6 conflict and 1 subset (tauroursodeoxycholate, where the structures are the same compound with partial stereo and the E label comes from the Broad row's name).

Broad rows whose own InChI disagrees with their stored InChIKey: 21. Two differ in connectivity (xylose, L-urobilin); the other 19 in the stereo block, often one side with stereo and the other without. For N-formylmethionine the Broad InChI gives the L form and the stored key the D form, which confirms that the key, not the compound, is wrong.

## v0.3 results

Evaluation against the manual review (53 groups), with structure checks on:

| outcome | v0.1 | v0.2 | v0.3 |
|---|---|---|---|
| agree | 28 | 29 | 34 |
| safe_miss | 14 | 13 | 8 |
| wrong | 1 | 1 | 1 |

Five D groups (cholic acid, cholesterol, deoxycholic acid, chenodeoxycholic acid, tauro-beta-muricholic acid) moved from R11 to R09b and now agree. No B group was called broader/narrower. The one wrong answer is still xylose.

Lipid matches: 60. The new one is cholesterol, whose two rows share Species ST 27:1;O at the same level; R09b now resolves it as broader_narrower, matching the review label D. Added to notes/v02_proposal_review.csv as correct: 60 of 60.

Row check C04 flags 21 Broad rows whose own InChI disagrees with their stored key, 2 in connectivity (xylose, L-urobilin).
### Input rows

Source: `results/rows_summary.csv`

| file | rows | with InChIKey | invalid InChIKey | no MetLinkR name |
|---|---|---|---|---|
| LEOCC_Metabolon | 1548 | 1228 | 0 | 427 |
| Broad_2022Aug | 511 | 210 | 1 | 9 |

### Groups sharing a MetLinkR name

Source: `results/pattern_counts.csv`

| pattern | groups |
|---|---|
| identical_inchikey | 124 |
| same_connectivity_different_block2 | 37 |
| protonation_only | 0 |
| different_connectivity | 16 |

Total groups with two or more rows: 177.
`protonation_only` does not appear in `pattern_counts.csv` because no group matched it, so it is recorded here as 0.

### Manual review

Source: `notes/review.csv`

| category | groups |
|---|---|
| A | 0 |
| B | 2 |
| C | 13 |
| D | 27 |
| E | 10 |
| F | 1 |

Examples (three per category where possible, with the PubChem record that settled it):

## Lipid feasibility test

Criteria and scope: `notes/decisions.md`, 2026-09-25. Scripts: `scripts/10_swisslipids_slim.R`, `scripts/11_lipid_feasibility.R`. Output: `results/lipid_feasibility/`.

| measure | result |
|---|---|
| LEOCC complex non-group rows with one Species ancestor (criterion) | 114 of 124, 91.9% |
| Same, with the connectivity fallback added after the first run | 117 of 124, 94.4% |
| Broad complex rows parsed by rgoslin to Species or finer (criterion) | 169 of 169, 100% |
| Disagreements with manual review labels (15 groups with lipid rows) | 0 |

Both criteria are met. SwissLipids and rgoslin go into the design as lipid evidence sources.

The Metabolon row "stearoyl sphingomyelin (d18:1/18:0)" and the Broad row "SM 36:1;O2" both resolve to the string "SM 36:1;O2", one through SwissLipids and one through rgoslin. Cross-vendor lipid names can be brought into one notation from database evidence alone.

The 7 LEOCC rows still without a Species:
- 5 have no SwissLipids entry at any tier: one hydroxylated SM, one d17:1 ceramide, three d18:2 SMs.
- 2 d18:2 ceramides match an entry, but no entry with the same connectivity has a Species parent. SwissLipids places d18:2 sphingoid bases poorly in its hierarchy.

Hierarchy vs hit abbreviation: 97 rows agree, 4 differ. All 4 are glycosylceramides where the hierarchy says HexCer or Hex2Cer and the hit says GlcCer or LacCer with the same chain sums. Comparisons should use the hierarchy's Species name, which is written at one level throughout.

Group-name rows: 29 of the LEOCC complex lipid rows report several molecular species in one name, mostly sphingomyelins. 20 of them map to one Species through their representative InChIKey. Not part of the criteria. Whether group rows can be matched at Species level is an open design question, because it relies on all members of a group sharing one chain sum.

### Pair rules on groups outside the review (173 pairs)

Source: classify_pairs() on results/rows.csv, groups listed in notes/review.csv excluded.

| rule | pairs | relationship |
|---|---|---|
| R03 | 136 | exact |
| R01 | 19 | unresolved |
| R02 | 6 | distinct |
| R05 | 6 | broader_narrower |
| R99 | 4 | unresolved |
| R08 | 2 | unresolved |

- R02: MetLinkR grouped three internal standards with their analytes (glycocholate-d4, inosine-15N4, thymine-d4), not only inosine.
- R05: all six are a Metabolon structure against a Broad Species of the same lipid. None came out exact, which confirms the lipid rules must run before R03.
- R99: four lysophospholipid pairs such as "1-palmitoyl-GPC (16:0)" against "LPC 16:0/0:0". Both have the same Species at structure level and the Broad row has no InChIKey. They probably match at sn-position level. Comparing sn-position names could resolve them, but when that is safe enough to call exact needs more thought. Left for after the evaluation.
- R08: two SM group rows against a Species-level SM. Both groups probably contain that Species. Deciding it needs the group-row check listed under Open questions in design.md.

## Evaluation of v0.1 against the manual review

Labels committed before the audit output was opened (commit: 33ececf). Script: `scripts/21_evaluate_v0.1.R`. Output: `results/evaluation_v0.1/`.

| outcome | groups |
|---|---|
| agree | 28 |
| safe_miss (tool says unresolved, review has an answer) | 14 |
| wrong (tool makes a different definite claim) | 1 |
| E groups not called exact | 10 of 10 |

Of the 42 groups where the review gives a definite relationship (B, C, D), the tool agrees on 28 (67%) and leaves 14 unresolved (33%). It makes no wrong definite claim on them. The one wrong result is xylose (F), where R12 calls the ring and open-chain forms distinct. This was listed under Known limits in design.md before the evaluation and has its own test.

Of the 13 wrong merges by MetLinkR (C), the tool finds 11 as distinct through R12 and leaves 2 unresolved (group-name rows, R08). None of the 10 groups with source errors (E) is called exact.

Where the 14 unresolved results come from:
- R11, 8 groups (6 D, 2 B): both keys carry stereo, and the key alone cannot tell partly specified stereo from different stereo (cholesterol, cholic acid, androsterone vs etiocholanolone). Needs the InChI stereo layers or SMILES.
- R08, 6 groups (2 C, 4 D): group-name rows. Needs a check of what a group contains.

Both gaps were listed in design.md before the evaluation. They are the two clearest targets for the next version.

## v0.2 feasibility: lipid match proposals

Source: `scripts/30_v02_feasibility.R`, reviewed by hand in `notes/v02_proposal_review.csv`. Criteria set in `notes/decisions.md` (2026-10-01) before the review.

| | proposals | correct | wrong |
|---|---|---|---|
| not linked by MetLinkR | 43 | 43 | 0 |
| also linked by MetLinkR | 10 | 10 | 0 |

The criterion (at most 4 wrong among the 43 new proposals) is met with no wrong proposal. Each proposal was checked from the names alone, independently of the tool: same lipid class, and the Metabolon chains add up to the Broad Species.

Harder cases handled correctly:
- Plasmalogens: Metabolon "P-16:0/18:1" against Broad "PC O-34:2". The plasmalogen double bond counts as one extra in O notation, so P-34:1 equals O-34:2. The relationship stays broader_narrower, since O-34:2 also covers other ether lipids.
- One-to-many: Broad "PC 34:4" matches both 14:0/20:4 and 16:1/18:3 from Metabolon, "PC 36:4" matches 18:2/18:2 and 16:0/20:4, and each DG Species matches both the 1,2- and the 1,3-isomer. These are reported as broader_narrower and one_to_many, never exact. For a meta-analysis, a Broad Species value should not be compared with each Metabolon row on its own.

On this data, v0.2 would add 43 lipid links to the 10 that MetLinkR made.

### Open questions

- Resolved: the one invalid InChIKey in Broad_2022Aug is row 354, "Inosine-15N4 [iSTD]", where the InChIKey field holds the literal text "inchikey". MetLinkR gave this row the standardized name "Inosine", the same as natural inosine (UGQMRVRMYYASKQ-KQYNXXCUSA-N in both files). The early exploratory count read "inchikey" as a separate connectivity block and reported 17 different_connectivity groups. The format check in 01_build_rows.R excludes it, which gives 16.
- Two points for later. First, a placeholder value sits in a structure column of the source file, which fits category E. Second, MetLinkR grouped an isotope-labeled internal standard with the natural metabolite. Internal standards are not biological measurements and should not be merged with the analyte in a meta-analysis. The rule set needs a way to catch these. Name tokens such as "[iSTD]" or isotope labels such as "15N4" are candidates, but that has to be checked against more files first.
- "(N(1) + N(8))-acetylspermidine" (LEOCC row 1028) got the MetLinkR name "Asparagine". Its source IDs are comma-separated lists for two isomers and none of them refer to asparagine. Test whether list-valued ID cells cause wrong mappings in MetLinkR.
- Some source names report a group ("(2 or 3)-decenoate") while the InChIKey column holds one specific member. When a name contains "or" or "/", its InChIKey should be treated as one candidate, not as the confirmed structure.
- A group can hold more than one problem at once. In the pregnenolone sulfate group, two rows are a correct match, one row has an InChIKey that does not fit its own name (E), and that same row was wrongly merged (C). The review sheet allows one category per group. Decide before phase 2 whether the tool reports relationships per row pair rather than per group.
- Xylose appears as both its ring form and its open-chain form. These are the same compound, yet the InChIKeys differ in block1. A different connectivity block therefore does not always mean a different compound. The tool needs a way to recognize ring-chain forms and tautomers of the same compound. Otherwise it will raise false alarms on correct merges. If more cases like this turn up, consider a separate category after review.
- Several source InChIKeys carry stereo or position detail that the row name does not (2-hydroxyglutarate as 2S, "(2 or 3)-decenoate" as the 2-isomer). The tool should compare the detail level of the name with the detail level of the key, and not trust the key's extra detail by default.
- Stereo specified against unspecified (R09) is labeled D, but for many natural metabolites (AMP, amino acids, sugars) only one stereo form occurs in biology, so the merge is right in practice. Consider whether the report should say so, for example by noting that R09 pairs are usually safe to merge for downstream analysis.
- In the bilirubin group, the LEOCC keys for "bilirubin (E,E)*" and "bilirubin (E,Z or Z,E)*" are swapped. Found only by opening each key's PubChem record separately: the order of results on the PubChem page did not match the order of the keys, so lookups in review must always be done one key at a time.
- Several pairs have one key with full stereo and one with only part of it (chenodeoxycholic acid, cholesterol). R09 catches only keys with no stereo at all (UHFFFAOYSA), so these reach R11 and come out unresolved. Telling "partly specified" from "differently specified" needs the InChI stereo layers, not just the key. Worth a look after the evaluation.
- Dehydroepiandrosterone sulfate: the LEOCC SMILES and InChIKey describe a stereoisomer with three inverted ring-junction centres, not DHEA-S. Found by recomputing the key from the source SMILES with RDKit and comparing InChI stereo layers. The same recomputation is what check C04 (v0.3) would do automatically. An R route for it (for example rcdk) needs to be tested before v0.3.
- Broad row "Deoxycholic acid" (the one with HMDB0000733, key DGABKXLVXPYZII-SIBKNCMHSA-N) is structurally hyodeoxycholic acid; its name is wrong in the source file. MetLinkR grouped it correctly with hyodeoxycholic acid. Broad also reports deoxycholic acid twice with the same key, which adds to the same-file duplicate question.
- Three LEOCC steroid structures have wrong stereo in the ring system: dehydroepiandrosterone sulfate (three ring-junction centres), glycochenodeoxycholate (C14 inverted, C10 unspecified), glycolithocholate (C13 and C17 inverted). In each case the key was missing from PubChem, and the error was found by recomputing the key from the source SMILES with RDKit and comparing stereocentres with a reference structure. A key that PubChem does not know, for a common metabolite, is a warning sign worth a row check once PubChem lookups exist (v0.3).
- Quinic acid: the Broad InChIKey differs from the key computed from the Broad SMILES in the same row. A second example, after xylose, of a source row whose SMILES and InChIKey disagree. Check C04 (v0.3) would flag both.
- PubChem synonyms can be wrong: CID 168408 is listed as tauro-alpha-muricholate, but its structure has the beta configuration at C6/C7. Check C05 (name vs PubChem synonyms, v0.3) must treat synonyms as weak evidence and should not flag a row on synonyms alone when the structure agrees with the row name.




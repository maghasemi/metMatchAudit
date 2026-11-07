# Findings

Numbers here come from the scripts, not from memory. When a number is copied in, note which results file it came from.

On 2026-10-07 the sections were put back in the order the work was done: the v0.2 to v0.4 results had been added at the top, under the phase 1a heading. No section's text was changed. A short v0.6 section was added, taken from the v0.6 entries in decisions.md, since that version had none here. The git history of this file shows the earlier layout.

## Phase 1a: InChIKey check on MetLinkR example output

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

## v0.4 results

Identifier resolution on the development data: the evaluation of the 53 reviewed groups is unchanged (34 agree, 8 safe_miss, 1 wrong), as expected, since every reviewed row already had a key.

Outside the review, 8 pairs that were unresolved for lack of a key (R01) now have a decision: unresolved pairs 24 to 16, exact 136 to 142, broader_narrower 7 to 9.

Three new lipid proposals, all correct: 1-palmitoyl-GPC (16:0) against LPC 16:0/0:0 (broader_narrower, R09), 1-palmitoyl-GPE (16:0) against LPE 16:0/0:0 and 1-stearoyl-GPE (18:0) against LPE 18:0/0:0 (exact, R03). These were three of the four lysophospholipid pairs left at R99 in v0.2. Proposals: 63 of 63 correct.

A first version resolved identifiers on Species-level lipid names and on internal standards, raising row flags from 15 to 93 and costing 8 rows their lipid profile. Both fixed before the version was frozen (decisions.md, 2026-10-06); row flags are now 17, two of them real disagreements between a vendor's identifier and its own structure (LPC 19:0/0:0, 2-aminobutyric acid).

## 2026-10-06: COMETS external test (blind)

### Setup

metMatchAudit v0.4.1 (tag `v0.4.1`) ran on MetLinkR output for five COMETS datasets: MetLinkR commit a6aa145, RaMP DB 3.0.7. The reference is `ManualMappings.csv` from the MetLinkR supplement. Criteria are in `decisions.md` (2026-10-06, COMETS test) and were committed before the first run. The first three attempts stopped before any comparison because of how expert labels were joined to rows (see the join fix entries in `decisions.md`). The fourth ran to the end. Its state is tagged `comets-test-1`. Script: `scripts/60_comets_test.R`. Output: `results/comets_test/`.

### Result

- 1429 rows. 401 cross-file pairs inside MetLinkR groups. 401 cross-file pairs the expert groups.
- Expert groups spanning two or more files that MetLinkR found whole: 183 of 214 (85.5%). The MetLinkR paper reports 85.3%, so our MetLinkR run reproduces theirs.
- C1, definite claims agree with the expert: 121 of 125 (96.8%). Met.
- C2, precision and F1 not lower with the audit: precision 0.895 in both, F1 0.893 against 0.895. Not met. The audit called two pairs distinct and the expert groups both.
- Counting only exact and broader_narrower as matches: precision 0.978, recall 0.339.
- 739 of 1429 rows got no structure key, because three files carry names only. 259 of 401 pairs ended unresolved.
- New lipid matches: 4, all in the same expert group.

| tool relationship | expert same | expert different |
|---|---|---|
| exact | 121 | 2 |
| broader_narrower | 15 | 1 |
| stereo_variant | 0 | 1 |
| distinct | 2 | 0 |
| unresolved | 221 | 38 |

### The four wrong definite claims

| pair | tool | cause |
|---|---|---|
| scyllo-inositol / myo-inositol | exact | Both identifiers resolved to the inositol key without stereo. Identical stereo-free keys were taken as exact. A gap in the rules. |
| DPA n3 / DPA n6 | exact | Identical keys for two compounds with different double bond positions. An identifier in one source points to the wrong isomer. |
| fructose / fructose | distinct | Ring form against open chain. Same known limit as xylose, now seen twice. |
| 4-vinylphenol sulfate / 4-vinylphenol sulfate | distinct | One identifier points to 4-vinylphenol, not its sulfate. |

Both wrong distinct calls, the ones that failed C2, are pairs with identical names. Identical names with different connectivity mean the identifiers conflict, and the tool cannot tell which one is wrong.

### Post hoc review of the 42 pairs the expert does not group

Done after the result was sealed, so it explains the result and does not replace it. Details in `notes/comets_review.csv`.

| review class | pairs |
|---|---|
| different_compound | 7 |
| stereo_variant | 4 |
| broader_narrower | 12 |
| aggregate (Nightingale class totals and fatty acids) | 5 |
| same_compound_platform_split | 9 |
| expert_inconsistent | 5 |

- No Biocrates acylcarnitine, nor Biocrates ornithine or SM C18:0, is linked to another platform. Each sits alone in its expert group. Biocrates lysoPC a C18:2 is linked to Metabolon, so this is not a rule applied to the whole file.
- The expert groups Nightingale Ala with beta-alanine, not alanine. The Nightingale biomarker list gives Ala as alanine.
- Two pairs with identical names (glutamate, N-acetylcarnosine) are split across groups. Their Metabolon COMP_IDs differ, but other pairs whose COMP_IDs differ between the same two files share one group.
- Only 7 of the 42 are different compounds. The audit called 2 of them exact and left 5 unresolved.

### What this means

- On COMETS the audit did not improve on MetLinkR. Its definite claims are mostly right, but it makes few of them, and it removed none of MetLinkR's wrong merges.
- The main cause is missing evidence: rows from name-only files have no key, so most pairs stop at R01.
- The expert answers "same measurement", the audit answers "same compound". Platform splits and class totals account for 14 of the 42 pairs.
- Changes this test points to: identical names with different connectivity become unresolved with a row flag; identical stereo-free keys from resolved identifiers are not enough for exact when names state different isomers; names that state different omega positions conflict; names joined by a slash count as group names; name to RefMet (v0.5); an aggregate relationship (v0.6).
- Any COMETS number after these changes counts as after adaptation. A fresh dataset with expert mappings is needed for a second blind test.

## 2026-10-07: v0.6 results (after adaptation on COMETS)

The plans and the acceptance check are in `decisions.md` (2026-10-07, v0.6 plan and LIPID MAPS plan). Script: `scripts/10_v0.6/04_acceptance.R`. Output: `results/v06_check/`.

- Development data: no pair changed rule. The evaluation is the same as for v0.4.2: 35 agree, 8 safe miss, 0 wrong, 10 E not exact. LMSD gave 233 rows a lipid profile.
- Lipid match proposals on the development data went from 63 to 110. All 47 new ones were checked by hand and are correct (`notes/v06_dev_lipid_review.csv`).
- COMETS: 23 pairs moved from R01 to a definite answer, 22 of them across files. Two became distinct, valerylcarnitine against glutarylcarnitine and C14:1 against C14:2, both MetLinkR merges of different compounds. Across files, definite claims agree with the expert in 122 of 122 (recounted with `scripts/08_comets/04_recount_all_versions.R`).

## 2026-10-07: second blind test (blind)

### Setup

metMatchAudit v0.6.0 (tag `v0.6.0`) with SwissLipids and LMSD ran on MetLinkR output for two files. MetLinkR was at commit a6aa145 with RaMP DB 3.0.7.

- A is the author's own Metabolon annotation of a mouse dataset, 494 rows. It is private.
- B is the annotation from Ding et al. 2021 (Nat. Commun., DOI 10.1038/s41467-021-26310-y), 1547 rows of mouse brain.

The tool had never seen B's naming style. Plan and criteria are in `decisions.md` (2026-10-07, second blind test). Label md5s and the group label order are in `notes/blind2_labels.md`. All of these were committed before the run. Script: `scripts/93_blind2_test.R`. Counts: `results/blind2/`. Row-level details stay in `data-raw/own/blind2/results/` because A is private.

### Result

| criterion | result | |
|---|---|---|
| C1: definite claims (exact, distinct) agree with the label, at least 90% | 48 of 65 (73.8%) | not met |
| C2: groups labeled A to D with a wrong definite answer, at most 5% | 12 of 58 (20.7%) | not met |
| C3: lipid proposals involving a sampled B row are in the labels, at least 90% | 1 of 1 | met, on one proposal |
| C3 recall | 1 of 5 | no threshold |

The audit has 628 pairs in 395 groups: 216 exact, 205 distinct, 157 unresolved, 47 broader_narrower and 3 stereo_variant. It proposed 57 lipid matches across files.

| label | tool: distinct | unresolved | broader_narrower | exact |
|---|---|---|---|---|
| A | 2 | 8 | 4 | 47 |
| C | 1 | 0 | 0 | 11 |
| D | 1 | 2 | 2 | 2 |
| F | 0 | 1 | 4 | 1 |

### Causes of the wrong answers

The 17 wrong definite claims and the 12 wrong groups come from four causes.

| cause | groups | wrong definite pairs |
|---|---|---|
| 1. B "isomer" features carry the parent's identifiers | 4 | 11 |
| 2. B's name gives a sum, its key a single structure | 3 | 3 |
| 3. Stereo in only one row's key, same compound by name | 4 | 0 |
| 4. Ring against open-chain key, names differ | 1 | 2 |
| F label: threonate against isothreonic acid | 1 | 1 |

1. B reports adenosine, carnosine, glucose 6-phosphate and N-acetylaspartate as features. It also reports further features named "Adenosine isomer", "Carnosine isomer A" and "B", "Hexose phosphate isomer D" and "N-Acetylaspartic acid isomer B". Each of these carries the CID and InChIKey of the parent. The identical keys led to `exact` (R03). R03c reads stereo terms and omega positions in names, but not the word "isomer". The parent's key on these rows is a placeholder for an unknown isomer. The tool cannot see that from the key.
2. B names acylcarnitines as sums ("ACar 17:0;"). Its InChIKeys are single structures, usually the straight chain. The tool follows the key, so straight-chain margaroyl and arachidoyl carnitine in A came out `exact`. For palmitoleoyl the B key is a different 16:1 isomer, so it came out `distinct`. The labels follow the name and say broader or narrower. The same thing happened with myristoylcarnitine in the lipid labels: the tool said exact, the label said D.
3. ADP-ribose, N-acetylmethionine and isovalerylcarnitine each have one key with stereo and one without. The names say the same compound, so the label is A. The tool says broader_narrower (R09), since one record is less specific than the other. Sphingosine 1-phosphate against erythro-sphingosine 1-phosphate got broader_narrower from the lipid rules (R05). These are not definite claims under C1, but they count as wrong in C2. This is the record-against-compound difference already seen with COMETS.
4. B:1 "(-)-N-Acetylneuraminic acid" has the open-chain key. Its connectivity differs from the ring form, so the tool calls it `distinct` from both N-acetylneuraminate rows. With fructose and ribose the same case stayed unresolved through R12n, because the names matched. Here the "(-)-" prefix keeps the names apart. This is the ring-form limit in Known limits, met for the third time.

### Lipid recall

B's lipid names did not parse. Of the 960 B rows matching the lipid class pattern, 485 got a profile from SwissLipids by key, 97 from LMSD, 2 from the SwissLipids fallback and 1 from the name. 375 got none. MS-DIAL-style names such as "PE 40:6;" and "PE 38:6; PE 18:2_20:4;" are not read by the name step.

So B rows got their profile from their key, at structure level, though the name gives only a sum or the chains.

- B:953 PE 38:6 against A:33 (18:2/20:4): both are structure level with different keys. The pair was not kept as a match.
- B:970 PE 40:6 against A:79 (18:0/22:6): same reason.
- B:575 LPE O-18:2 against A:12 "1-(1-enyl-oleoyl)-GPE (P-18:1)*": A:12 got no profile. It has no identifier, and its name was not parsed.
- B:1037 "PE O-38:5; PE O-18:1_20:4; or PE P-38:4; PE P-18:0_20:4;" got no profile. That is right for a row that names two alternatives.

### What this means

- On data with a naming style it had not seen, the tool's definite claims were wrong a quarter of the time. That is far from the 96.8% on COMETS and the 100% after adaptation.
- Most of the error comes from trusting identifiers that the source itself contradicts in the name: placeholder keys on "isomer" features, and structure keys on sum names. The tool checks names against keys for stereo terms and omega positions only.
- The labels answer from the name wherever name and key disagree. Read by the key, several of these pairs would be called the way the tool called them. The labels were fixed first and are not changed here.
- Changes this points to, all of them after adaptation if made:
  - a name stating "isomer" (or "isomer A", "isomer B") blocks `exact` from identical keys, as R03c does;
  - a lipid name that parses below structure level caps the level a key can give;
  - read MS-DIAL names: split on ";", take the most specific annotation, and treat "or" as an alternative;
  - strip optical rotation prefixes such as "(-)-" in name normalization.

### Diagnostic: v0.4.1 on the same files

Done after the blind result was committed, so this is not a blind result. The question is whether the changes made after v0.4.1 (v0.4.2 and v0.6.0) made things worse on data that neither version was built on.

- Script: `scripts/94_blind2_v041_diagnostic.R`.
- Counts: `results/blind2_v041_diagnostic/`.
- v0.4.1 got the PubChem cache of the v0.6.0 run, so both versions saw the same identifiers.
- The script recomputes the v0.6.0 scores and checks that they match `scripts/93_blind2_test.R`.

| measure | v0.4.1 | v0.6.0 |
|---|---|---|
| C1: definite claims that agree | 48 of 67 (71.6%) | 48 of 65 (73.8%) |
| C2: groups labeled A to D with a wrong definite answer | 13 of 58 | 12 of 58 |
| C3: lipid proposals involving a sampled B row, in the labels | 0 of 0 | 1 of 1 |
| C3 recall | 0 of 5 | 1 of 5 |
| lipid matches across files, all | 16 | 57 |

Three labeled pairs changed rule:

- Fructose and ribose against their furanose keys went from `distinct` (R12) to unresolved (R12n). R12n was added in v0.4.2 after COMETS. Here it removed two wrong claims on data it had not seen.
- Sphingosine 1-phosphate against erythro-sphingosine 1-phosphate went from unresolved (R99) to broader_narrower (R05). In v0.6.0 the B row gets a lipid profile, and the pair is decided by the lipid rules. The label is A, so this group went from safe miss to wrong.

That label is not consistent with another label in the same sheet. Sphinganine against erythro-sphinganine is labeled D, on the grounds that only one name states the stereo. The S1P pair has the same difference in its names, but it was labeled A because both rows give the same CID. Under the rule used for sphinganine, the S1P pair would be D, and the tool's answer would agree. The committed label stays A, and C2 stays at 12.

So the versions after v0.4.1 did not make things worse here. They removed two wrong claims, added one wrong group that is open to doubt, and turned up 41 more lipid matches. The drop from COMETS to this test comes from the data, not from the later versions. Every cause listed above was already present in v0.4.1.

## 2026-10-07: v0.7 acceptance (after adaptation)

Script: `scripts/96_v07_check.R`. Log: `results/v07_check_log.txt`. Every run used a copy of the PubChem cache of the v0.6.0 run it was compared with.

| set | result |
|---|---|
| development | 0 rows with a marker, 0 pairs changed. Evaluation 35 agree, 8 safe miss, 0 wrong, 10 E not exact, as before. 110 lipid proposals, none lost or gained. |
| COMETS | 0 rows with a marker, 0 pairs changed. Definite claims across files agree with the expert in 122 of 122 for both versions. |
| second blind test files | 155 rows with a marker, all in B: 55 by the word "isomer", 100 by a lone letter (49 A, 43 B, 8 C). 193 pairs changed rule, all to R03i except one to R12n. |

COMETS is 122 of 122 here and 123 of 123 in the README. Both versions were counted the same way in this script, so v0.7 changes nothing there. The difference came from how the earlier numbers were worked out. It was since checked by scripts/97_comets_recount.R: 122 of 122 is right (decisions.md, 2026-10-07, COMETS numbers for v0.6 corrected).

Every lone-letter name was checked by eye. All are MS-DIAL lipid features ("PC 36:2; A", "Cer-NS d42:2; Cer-NS d18:1/24:1; B").

Second blind test files, scored as in `scripts/93_blind2_test.R`:

| measure | v0.6.0 (blind) | v0.7.0 (after adaptation) |
|---|---|---|
| C1 | 48 of 65 (73.8%) | 48 of 53 (90.6%) |
| C2 | 12 of 58 | 8 of 58 |
| C3 | 1 of 1 | 1 of 1 |
| C3 recall | 1 of 5 | 1 of 5 |

How the labeled pairs changed:
- 11 C pairs, all isomer features, went from `exact` to unresolved.
- 5 F pairs, the 2-hydroxyglutarate isomers, became unresolved.
- B:1 against B:627 (Neu5Ac) went to R12n.

This matches the expectations written in the plan, with the Neu5Ac correction from the feasibility entry. The acylcarnitine pairs did not change, as expected once MS-DIAL name reading was dropped.

The 8 groups still wrong:
- four with stereo in one record only (cause 3);
- three acylcarnitines (cause 2);
- N-acetylneuraminate in A against B:1 (cause 4, names still differ);
- threonate, which is labeled F.

The 193 changed pairs include 62 that were `exact` and 43 that were `distinct` under R06. Most are not labeled. They move to unresolved, which is the safe direction but makes the tool answer less often on files like B.

v0.7.0 meets its acceptance and is tagged.

## 2026-10-07: third blind test (blind)

### Setup

metMatchAudit v0.7.0 (tag `v0.7.0`) with SwissLipids and LMSD ran on MetLinkR output for three files:
- A, the private Metabolon file from the second test, 494 rows;
- C1, the fecal list of Lai et al. 2021 (Nat. Commun., DOI 10.1038/s41467-021-26209-8), 533 rows;
- C2, the serum list of the same paper, 231 rows.

MetLinkR was at commit a6aa145 with RaMP DB 3.0.7. Its mapping rates were 97.4% for A, 78.2% for C1 and 89.6% for C2.

Plan and criteria: `decisions.md` (2026-10-07, third blind test). Label md5s and the stereo rule: `notes/blind3_labels.md`. All of these were committed before the run.

The labels: 60 groups drawn from 175 that span at least two files, 94 pairs (91 A, 3 D), and 40 lipid and carnitine rows. Script: `scripts/103_blind3_test.R`. Counts: `results/blind3/`.

### Result

| criterion | result | |
|---|---|---|
| C1: definite claims (exact, distinct) agree with the label, at least 90% | 78 of 78 (100%) | met |
| C2: groups labeled A to D with a wrong definite answer, at most 5% | 7 of 60 (11.7%) | not met |
| C3: lipid proposals involving a sampled row are in the labels, at least 90% | 13 of 13 | met |
| C3 recall | 13 of 20 (65%) | no threshold |
| C1 on pairs without A | 30 of 30 | reported |

| label | exact | broader_narrower | unresolved |
|---|---|---|---|
| A (91) | 78 (R03) | 11 (R09) | 2 (R11) |
| D (3) | 0 | 0 | 3 (R01) |

There was no wrong `exact` and no `distinct` call at all. The causes that drove the second test's failure, placeholder keys on isomer features and sum names with single-structure keys, did not occur in these files.

### Why C2 failed

All 7 wrong groups are labeled A, and all were called broader_narrower by R09: same connectivity, one key without stereo. That makes 11 pairs:

| group | pairs | names |
|---|---|---|
| 3-hydroxybutyrate | 1 | neither name states stereo |
| CAR 12:0 | 2 | neither name states stereo ("laurylcarnitine", "dodecanoylcarnitine") |
| CAR 18:0 | 1 | neither name states stereo |
| glycerate | 1 | neither name states stereo |
| N-acetylglutamate | 2 | neither name states stereo |
| CAR 16:0 | 2 | one name says L, the natural form |
| malate | 2 | one name says L, the natural form |

In five of the seven groups, neither name says anything about stereo. The keys differ only because one source's PubChem CID is a record without stereo and the other's is a record with it. So the R09 call reflects which record each lab picked, not the compound measured.

This also means the result does not depend on the stereo rule committed for this test. Under the rule used for sphinganine in the second test, where one name stating stereo makes the pair D, only CAR 16:0 and malate would change. C2 would then be 5 of 60, still not met.

This is the same cause as cause 3 in the second test (ADP-ribose, N-acetylmethionine, isovalerylcarnitine), now seen on a second, independent pair of files.

### The other results

- The two unresolved A pairs (malonylcarnitine, N-acetylneuraminate) went to R11: keys with different stereo, and no name states any. These are safe misses.
- The three D pairs involve A rows with no identifier ("aconitate [cis or trans]", "oleate/vaccenate (18:1)"), so R01 applied. These are safe misses too.
- All 13 lipid proposals were right pairs. One, A's 1-palmitoyl-2-oleoyl-GPS against C2's PS(16:0/18:1) with the same CID, came out broader_narrower (R05, structure against sn position), while the label says A. C2's name does not give the double bond position, so the tool's reading can be defended. C3 counts pairs, not relationships, so this does not affect it.
- The 7 labeled lipid pairs the tool missed:
  - three are free carnitine, which gets no lipid profile;
  - three have one row without a profile ("PE(P-16:0e/0:0)", "dihomo-linoleoylcarnitine (C20:2)*" without a CID, trans-2-tetradecenoylcarnitine);
  - one is malonylcarnitine, which R11 left unresolved.

### What this means

- On a naming style the tool had not seen, every definite claim was right.
- C2 failed on one thing: broader_narrower from R09 where the two names agree and only the chosen PubChem record differs in stereo. Across the two blind tests, this is the most common way the tool disagrees with a curator.
- A change this points to, after adaptation if made: when neither name states stereo, a stereo-only key difference should not give broader_narrower. Unresolved is the conservative answer. Reading both blind tests back under such a rule is a diagnostic, not a blind result, and it needs its own plan first.

## 2026-11-07: v0.9 feasibility

Source: `results/v09_feasibility/`, from `scripts/15_v0.9/01_feasibility.R` on v0.7.1. The package did not change.

| check | result | needed |
|---|---|---|
| empty `Metabolite` rows of Broad whose `name` parses | 119 of 119 | at least 80% |
| same Species, where both columns are filled and both parse | 188 of 189 (99.5%) | at least 95% |

The 119 names parsed at three levels: 105 species, 11 sn position, 3 structure.

392 Broad rows have both columns filled, and 183 of them have identical text. The second check covers the 189 where both columns parse. The one disagreement is `Metabolite` "PC 33:1" against `name` "PC 30:1". That row keeps its `Metabolite` under the plan, so the change does not touch it, but it shows that `name` is not a plain copy of the other column.

### Where the blanks are

Broad is the only file with blank names: 119 of 511 rows. The other 11 entries in `blank_names.csv` have none, and that covers LEOCC, the five COMETS files and all the blind test files.

So this change cannot affect COMETS or the blind test files. The 234 R01 pairs on COMETS have another cause: rows without a key, not rows without a name.

### The 13 development pairs

All 13 are R01 today, and each has one Broad row with an empty name. From the rules and the v0.7 rows, I expect these results after the change. This is a reading of the code, not a run.

| pairs | Metabolon side | Broad side after the change | expected |
|---|---|---|---|
| 3 LPC | structure | sn position | R05, broader_narrower |
| 7 SM | group | species | R05, broader_narrower |
| 3 SM (`LEOCC_Metabolon:862`, `:988`, `:990`) | a key, no lipid profile | a profile, no key | R08, unresolved |

That would take development unresolved from 39 to 29. The last three leave R01 but stay unresolved.

None of the 13 groups is in `notes/review.csv`, so the evaluation against the review cannot show this gain. Only the hand check in the acceptance step can.

### A detail for the build

Empty cells in these files come through `fread()` as empty text, not as `NA`, and the same goes for `pubChemId` and `hmdbId` on these rows. The reader has to treat both as blank. The 119 rows also have no key and no identifier, so no PubChem lookup can follow from filling their names.

## 2026-11-07: v0.9 acceptance, automatic part (after adaptation)

Source: `results/v09_check/` and `results/v09_check_log.txt`, from `scripts/15_v0.9/02_acceptance.R` on branch `exp/name-fallback`, version 0.9.0. The unit tests passed first: 345 passed, 0 failed.

| check | result | needed |
|---|---|---|
| development pairs without the fallback | 245 of 245 identical | identical |
| COMETS pairs without the fallback | 414 of 414 identical | identical |
| wrong in the evaluation | 0 | 0 |
| source-error groups called exact | 0 | 0 |
| agree in the evaluation | 35 | at least 35 |
| changed pairs that hold a filled row | 13 of 13 | all |
| lipid proposals lost | 0 | 0 |

All 119 filled rows got a lipid profile from their name. None got a key.

Development pairs, before and after:

| | exact | broader_narrower | stereo_variant | distinct | unresolved |
|---|---|---|---|---|---|
| before | 146 | 35 | 3 | 22 | 39 |
| after | 146 | 44 | 3 | 22 | 30 |

Lipid proposals went from 110 to 129. Nine of the 19 new ones are the nine pairs that moved to R05, and ten are pairs MetLinkR did not group.

### One prediction was wrong

I expected 10 pairs to reach R05 and 29 unresolved. The result was 9 and 30. "1-stearoyl-GPC (18:0)" against "LPC 18:0/0:0" went to R99 and stayed unresolved.

The other two LPC names read at sn position level, but "LPC 18:0/0:0" reads at structure level. With no double bond, the name leaves nothing open. "LPC 20:0/0:0" reads the same way. So both rows of the pair are at structure level, and R05 needs two different levels. Only the Metabolon row has a key and only that row has no molecular annotation, so the key rules and R06 do not apply either.

The pair is a safe miss. The older R99 pair, "LPC 24:0/0:0", has the same shape. Closing this gap needs a new pair rule, which v0.9 does not add, so it stays a separate idea.

### The sheet check

Source: `results/v09_check/hand_check_evidence.csv`, from `scripts/15_v0.9/03_hand_check_evidence.R`. The assistant did this check, not me (see `notes/decisions.md`).

| check, on the 23 pairs | result |
|---|---|
| every alternative in the Metabolon name adds up to the chain sum in the Broad name | 23 of 23 |
| formula of the Metabolon SMILES equals the formula worked out from the Broad name | 23 of 23 |
| claims where the Metabolon row has more detail | 19 of 19 |

Under the counting rule, 19 of 19 claims are right. The 4 unresolved pairs are safe misses. All four fit the same lipid by chain sum and formula, so a stronger rule could answer them.

Three claims are structure against sn position (LPC 16:1, LPC 20:4, LPC O-16:1). Broad gives no double bond position, so broader_narrower can be defended, but a curator might call these the same compound. The third blind test had a similar case with PS(16:0/18:1).

Broad has no m/z and no formula on these rows. The check shows that the names and structures fit together. It does not show that Broad measured the same compound.

## 2026-11-07: v0.10 feasibility

Source: `results/v10_feasibility/` and `results/v10_feasibility_log.txt`, from `scripts/16_v0.10/01_feasibility.R` on version 0.9.0. The package did not change. Only development and COMETS names were sent to PubChem.

| condition | result | needed | |
|---|---|---|---|
| one key for rows that already have a key | 1179 of 1305 (90.3%) | at least 50% | pass |
| of those, same connectivity as the row's key | 1123 of 1179 (95.3%) | at least 98% | **fail** |
| COMETS R01 pairs that would become exact | 142 of 234 | at least 80 | pass |
| of those, agree with the expert | 140 of 142 (98.6%) | at least 98% | pass |

### Three runs

The script was run three times. The numbers above are from the third run, the only one in which no request failed.

- First run: 21 name requests failed. Same connectivity was 1110 of 1165 (95.3%), and 141 COMETS pairs became exact.
- Second run: no name request failed, but 198 rows got a CID and no key. A request for parents or keys that fails is not cached, but within the run it becomes a missing key, or the CID is used as its own parent, and nothing reported it. The script now counts every failed request and gives no go-ahead while any are left (commit `beae16b`).
- Third run: 0 failed requests. The result is close to the first run.

The same silence is in `resolve_identifiers()` in the released package. A rerun repairs it, since failed answers are not cached, but the user is not told to rerun.

### The 56 rows with different connectivity

The check compares PubChem's answer with the row's own key, and the row's key is not always right. These groups come from the names alone. I did not check the structures.

- 19 sulfates, most of them Metabolon names.
- 7 choline esters and choline itself, and 4 acylcarnitines. The name key often ends in `-O` or `-P`, so PubChem picked a charged record.
- 6 sugars (xylose, fructose three times, erythrose, maltotriose), probably the open chain against the ring form.
- 3 bile acids (glycocholic, deoxycholic and taurodeoxycholic acid), where the file key may be the wrong one.
- 17 others, among them 4-hydroxyphenylacetate, where PubChem seems to read the name as the ester.

By set: 713 of 759 on the development data (93.9%) and 410 of 420 on COMETS (97.6%).

### The COMETS R01 pairs

| outcome | expert same | expert different | no expert label |
|---|---|---|---|
| exact | 140 | 2 | 0 |
| both keys from names (would be R01n) | 34 | 1 | 0 |
| same connectivity, different stereo | 15 | 0 | 0 |
| different connectivity, same name | 1 | 0 | 0 |
| different connectivity, different names | 2 | 1 | 0 |
| no key on a side | 10 | 28 | 0 |

Rows without an expert label count as not covered. None of these 234 pairs lacks one.

The two exact pairs the expert separates are "Ornithine ( Orn-PTC )" against "ornithine" and "N-acetylcarnosine" against "N-acetylcarnosine". In both, the names point to the same compound.

Of the three pairs with different connectivity and different names, the expert calls two the same compound: "indolepropionate" and "Indole propionate", once in each direction. PubChem returns two different structures for two spellings of one name. Distinct would have been allowed only at 98% or more, and this is 1 of 3, so those pairs would stay unresolved.

## 2026-11-07: v0.9.1 acceptance

Source: `results/v091_check/` and `results/v091_check_log.txt`, from `scripts/17_v0.9.1/01_acceptance.R` on branch `release-prep/v0.9.1`, version 0.9.1. The unit tests passed first: 367 passed, 0 failed.

| check | development | COMETS |
|---|---|---|
| failed PubChem requests | 0 | 0 |
| keys and their sources identical to 0.9.0 | yes | yes |
| pairs identical to 0.9.0 (245 and 414) | yes | yes |

The runs used copies of the v0.9 caches, so no request had to go to PubChem. This shows that nothing changes when no request fails. The behaviour when requests fail is covered by the unit tests only, with a fake PubChem that answers `PUGREST.ServerBusy`.

## 2026-11-07: v0.11 feasibility

Source: `results/v11_feasibility/` and `results/v11_feasibility_log.txt`, from `scripts/18_v0.11/01_feasibility.R` on the v0.9.1 audits. The package did not change, and nothing was sent to PubChem.

| condition | result | needed |
|---|---|---|
| COMETS R01 pairs that would get `same_name` | 121 | at least 80 |
| of those, in one expert group | 120 of 121 (99.2%) | at least 98% |
| pairs with two keys, different files and the same name, whose keys do not conflict | 186 of 189 (98.4%) | at least 98% |

### The R01 pairs

| set | R01s | different names | same name in one file | same name, flagged name |
|---|---|---|---|---|
| COMETS | 121 | 100 | 8 | 5 |
| development | 0 | 3 | 0 | 0 |

The diagnostic count from before the plan was 126. The difference is the five pairs with a flagged name, which stay R01: three tentative Metabolon names marked with `*` ("bilirubin (E,E)*", "pyroglutamine*", "2-palmitoyl-GPC (16:0)*") and two names that list alternatives ("dihomo-linolenate (20:3n3 or n6)", "dimethylarginine (SDMA + ADMA)"). The expert puts all five in one group. Leaving them out costs five pairs and is the cautious side.

The script's "same file" count is 8, but none of those 8 has the same name on both sides.

112 of the 121 R01s pairs have a key on one side. The other row has no key, so the pair cannot be checked by structure.

The one pair the expert separates is "N-acetylcarnosine" against "N-acetylcarnosine", one side with a key. The same pair was the only exact disagreement in the v0.10 feasibility. The expert puts the two rows in different groups, so either the expert separates two measurements of one compound, or one row is not what its name says. The data cannot tell which.

### The pairs with two keys

| rule | COMETS | development |
|---|---|---|
| R03 (exact) | 119 | 37 |
| R09 | 16 | 6 |
| R09b | 0 | 1 |
| R11 | 6 | 0 |
| R12n | 2 | 1 |
| R99 | 1 | 0 |

The three conflicts are R12n pairs: xylose against Xylose, fructose against fructose and 4-vinylphenol sulfate against itself. The expert puts the two COMETS pairs in one group. All three look like two records for one compound, such as a sugar's ring and open chain, rather than two compounds under one name. None of the 189 pairs is distinct. On COMETS the expert puts 143 of the 144 in one group.

This check passed with little room: 98.4% against 98%. One more conflict among these 189 pairs would have failed it.

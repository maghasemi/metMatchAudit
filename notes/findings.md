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

## v0.4 results

Identifier resolution on the development data: the evaluation of the 53 reviewed groups is unchanged (34 agree, 8 safe_miss, 1 wrong), as expected, since every reviewed row already had a key.

Outside the review, 8 pairs that were unresolved for lack of a key (R01) now have a decision: unresolved pairs 24 to 16, exact 136 to 142, broader_narrower 7 to 9.

Three new lipid proposals, all correct: 1-palmitoyl-GPC (16:0) against LPC 16:0/0:0 (broader_narrower, R09), 1-palmitoyl-GPE (16:0) against LPE 16:0/0:0 and 1-stearoyl-GPE (18:0) against LPE 18:0/0:0 (exact, R03). These were three of the four lysophospholipid pairs left at R99 in v0.2. Proposals: 63 of 63 correct.

A first version resolved identifiers on Species-level lipid names and on internal standards, raising row flags from 15 to 93 and costing 8 rows their lipid profile. Both fixed before the version was frozen (decisions.md, 2026-10-06); row flags are now 17, two of them real disagreements between a vendor's identifier and its own structure (LPC 19:0/0:0, 2-aminobutyric acid).

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

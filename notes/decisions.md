# Decisions

One entry per decision, with the reason. Add new entries at the bottom and never rewrite old ones. If a decision changes, add a new entry that says so.

On 2026-10-07, before the repository was made public, the entries were put back in the order they were written. A few paragraphs that had been pasted under the wrong heading were moved to their own entry, and the first v0.5 plan, which had sat in a separate file, was moved in. No entry's text was changed. The git history of this file shows the earlier layout.

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

## 2026-09-30: Pair rules adjusted while writing their tests

Each rule was tested on a real pair from the review. Five changes came out of that:

1. R02 now runs before R01. "Inosine-15N4 [iSTD]" has the text "inchikey" as its key, so R01 would have called the pair unresolved, while the name alone shows it is an internal standard.
2. R08 now gives unresolved instead of group_overlap. The SM group summing to SM 38:3 would have been called overlapping with SM 36:1, which it is not (manual label C). v0.1 has no way to check what a group contains, so it makes no overlap claim. group_overlap is no longer used.
3. New rule R06b: two rows that report the same Species (or the same molecular species) at the same level are exact. Without it, "SM 36:1;O2" from two vendors ended up unresolved.
4. R10 and R11 now require block2 to differ. Otherwise two keys that differ only in the protonation letter would be called stereo variants. Such pairs now reach R99.
5. The headgroup case is its own rule, R04h, instead of a silent skip inside R04, so it shows up in the output with an ID.

## 2026-09-30: Row check C02 dropped

C02 would flag a row whose InChIKey has stereo while its name states none. Most amino acids, sugars and nucleosides fit that description ("inosine", "glucose"), so it would flag a large share of all rows and carry no signal. The case where it matters, two keys with different stereo and no stereo in either name, is already handled by pair rule R11.

C01 flags only values that fail the format check. A missing key is not flagged, since many rows have none.

## 2026-09-30: v0.2 feasibility and criteria, set before reviewing the proposals

Feasibility run (scripts/30_v02_feasibility.R): 62 cross-file lipid pairs share a Species. MetLinkR linked 15 of them. Of the 47 it did not link, 43 are proposed as broader_narrower (R05), 2 are distinct (R12), 2 unresolved (R99). Proposed matches: 43 new (29 one-to-one, 14 one-to-many) against 10 already linked by MetLinkR.

Group-name lipid rows: all 28 Metabolon rows with a chain list have members that share one chain sum. For 18 of them SwissLipids gives a Species through the representative InChIKey, and in all 18 it matches the sum computed from the name. v0.2 gives a group row a lipid profile only when both routes agree. The other 10 stay without one, since their class would have to be read from the name text.

Precision check: all 43 new proposals, plus the 2 distinct and 2 unresolved unlinked pairs, are reviewed by hand before any v0.2 code is written. Each gets one label: correct (the relationship holds), wrong (it does not), or unclear.

Criteria:
- v0.2 goes ahead if at most 4 of the 43 proposals (about 10%) are wrong.
- Any proposal of the form "exact" that is wrong stops the plan until the cause is found. None of the 43 are exact, so this applies to later runs.
- Every wrong proposal gets its cause written down, and each cause becomes a test.

Recall on MetLinkR's own lipid links is not a real test here, since v0.2 finds them with the same profiles v0.1 used to confirm them. Instead, the report states how many lipid rows have no profile and so cannot be matched at all.

## 2026-10-01: v0.2 code, group rows and R05g

Group-name lipid rows now get a profile at a new level, "group", when the chain sum shared by all listed members matches the Species SwissLipids gives for the row's key. Members written with ether, plasmalogen or hydroxy notation are not summed, since those counts follow different rules.

New rule R05g: a group row against a row at molecular level or finer with the same Species is unresolved, because v0.2 does not check whether that structure is one of the listed members. Against a Species-level row, R05 gives broader_narrower as before.

Side effect on the v0.1 evaluation: the SM 36:1 review group, labeled C, now comes out distinct through R04 instead of unresolved through R08.

propose_lipid_matches() pairs lipid rows from different files that share a Species and keeps exact and broader_narrower pairs. It does not depend on MetLinkR.

Correction (2026-10-01): the first count of group rows with a SwissLipids Species was 18. A recount of results/v02_feasibility/group_rows.csv gives 19, which matches the package output.

## 2026-10-01: v0.3 feasibility criteria, set before running

Question: can PubChem PUG REST turn a SMILES into InChI and InChIKey reliably, including for structures PubChem has no record of?

Reference: scripts/40_rdkit_reference.csv, InChIKeys computed from the same SMILES with RDKit. On the 60 sample rows, the RDKit key agrees with the key stored in the source file on the first 25 characters in 59 cases; the exception is the Broad quinic acid row, whose stored key does not match its own SMILES (known from the review).

Criteria, comparing PubChem with RDKit on the first 25 characters (the protonation letter is ignored):
- All 5 named rows must agree. Four of them are structures PubChem has no record for. If PubChem cannot return an InChI for them, it cannot serve the main purpose of v0.3.
- At least 57 of the 60 sample rows must agree, and every disagreement must be explained.

If either criterion fails, the PubChem route is dropped and a local toolkit (rinchi or rcdk) is tested next.

Request time and coverage (rows with a SMILES or InChI, R11 pairs with a structure on both sides) are recorded but are not criteria.

## 2026-10-01: Local toolkit for v0.3

ChemmineOB (Bioconductor) did not build: it needs the Open Babel system library, which every user would have to install separately. rinchi is no longer on CRAN.

Proposed: RDKit through reticulate, as an optional dependency. Only the v0.3 structure checks need it; v0.1 and v0.2 work without it. RDKit builds InChI with the IUPAC InChI library.

The reference keys were made with RDKit, so testing RDKit in R against them only checks the R-to-Python bridge. That RDKit itself gives correct keys rests on independent sources: on the 60 sample rows its keys agree with the keys stored by Metabolon and Broad in 59 cases (the exception is the Broad quinic acid row, whose stored key does not fit its own SMILES), and with PubChem on every structure PubChem had a record of, except one inositol stereo block (1-stearoyl-2-linoleoyl-GPI) that PubChem most likely changed during standardisation.

## 2026-10-01: v0.3 stereo comparison, design and RDKit feasibility

RDKit through reticulate works on the development machine (py_require("rdkit")). It reproduced the DHEA-S key that PubChem could not compute.

Stereo comparison (scripts/41_inchi_stereo_functions.R) reads the /t, /m, /s and /b layers of two InChIs with the same main layer. /m1 flips every sign in /t. Undefined centres ("?") are ignored. Result: same, subset, conflict, or relative (any /s other than 1).

Self-test on 317 LEOCC structures with stereo, run with RDKit: removing some stereo marks gave subset in 314 cases and conflict in 3 (symmetric molecules such as inositol, where fewer defined centres change InChI's canonical numbering). Inverting one centre gave conflict in all 317 and subset in none.

Design consequence: only "subset" creates a new claim, through new rule R09b (same connectivity, one stereo set is a subset of the other, broader_narrower). "conflict" creates no new claim: those pairs continue to R10 and R11 as before. The known false conflicts on symmetric molecules can therefore only leave a pair unresolved.

A row's InChI is used only when its key agrees with the row's stored key on the first 25 characters. Otherwise the row is a C04 case and its structure is not trusted.

Transparency: this script was already run on the example data in a separate environment, without notes/review.csv, before these criteria were written. What it showed: 65 of 65 reference keys reproduced; among the 16 R11 pairs, every reviewed D pair with a usable structure came out subset, androsterone against etiocholanolone came out conflict; 21 Broad rows have an InChI that disagrees with their stored key, 2 of them in connectivity (xylose, L-urobilin). The run on the development machine is a confirmation, not a blind test.

Criteria for the confirmation run:
- 65 of 65 reference keys reproduced.
- No R11 pair from a group labeled B is classified as subset.
- Broad xylose is flagged with a connectivity difference between its InChI and its stored key.

## 2026-10-06: Identifier resolution, criteria set before running

Reason: the five COMETS datasets with expert mappings carry no InChIKey or SMILES, only names and in part PubChem or HMDB IDs. Many real users are in the same position. The tool therefore needs a way to get an InChIKey from an identifier.

Design (scripts/50_id_functions.R):
- PubChem CID: InChIKey of that record from PubChem.
- HMDB ID: CIDs that PubChem lists for the HMDB registry ID, then their InChIKeys. Resolved only if all those CIDs give the same key.
- Cells listing several IDs are not resolved.
- HMDB IDs are normalised to seven digits ("HMDB00641" becomes "HMDB0000641").
- Every answer is cached on disk.
- A derived key describes the database record, not necessarily what was measured. Rules such as R11 already handle keys that claim more stereo than the name.

The feature is built and tested on the development files (LEOCC, Broad), where rows carry identifiers and an InChIKey, so the stored key is the reference. The COMETS data play no part in building it.

Criteria:
- CID route: at least 95% of resolved rows agree with the stored key on the first 14 characters.
- HMDB route (random sample of 300, seed 20261006): at least 90%.
- Every disagreement on the first 14 characters is listed and explained.

Agreement on the first 25 characters is reported but is not a criterion, since database records often carry more stereo than the source files.

## 2026-10-06: Identifier resolution, first run

HMDB route: 282 of 297 resolved rows agree with the stored key on the first 14 characters (94.9%). Criterion 90%: met.

CID route: 1098 of 1158 (94.8%). Criterion 95%: not met, by three rows. The criterion stands as written.

Reading the 60 CID disagreements by name (not yet checked record by record): 3 are group-name rows, which the tool should never have resolved (R08 already distrusts their keys); at least 2 are source errors found in the review (pregnenediol sulfate, andro steroid monosulfate); 3 are sugars in ring against open-chain form; many look like salts or charged forms (choline, thiamin, heme, acylcholines, several drugs); the rest look like rows where the CID and the InChIKey given by the vendor describe different compounds. In most of these cases PubChem returned the record for the CID correctly; the CID and the key in the same Metabolon row disagree. About 5% of LEOCC rows with both identifiers carry inconsistent ones.

Decision pending: drop the CID route, or revise the design (skip group and alternative names; compare parent compounds for salts) and test the revised design on data not used so far.

## 2026-10-06: Identifier resolution retest, criteria set before running

Test data: annotations.csv from the MetLinkR SI (2574 rows, 1535 with both a PubChem CID and an InChIKey). Never used in development, and not one of the five files with expert mappings. Only identifier resolution is tested; no grouping is looked at.

Rows are kept only if neither their CID nor the first 14 characters of their key occur in LEOCC or Broad, because the revised design was made after reading disagreements there.

Design as revised (2026-10-06): group and alternative names not resolved; key taken from PubChem's parent compound.

Criteria, unchanged from the first run:
- CID route: at least 95% of resolved rows agree with the stored key on the first 14 characters.
- HMDB route (random sample of 300, seed 20261006): at least 90%.
- Every disagreement listed.

If the CID route fails again, it is dropped from v0.4 and only the HMDB route goes in.

## 2026-10-06: Identifier resolution retest, result

On 783 fresh rows of annotations.csv (no CID or key block shared with LEOCC or Broad):
- CID route: 779 of 779 resolved rows agree on the first 14 characters (100%); 3 skipped for group or alternative names, 1 unresolved. Criterion met.
- HMDB route: 288 of 297 (97%); 2 skipped for names. Criterion met.

Both routes go into v0.4.

Caveat: the Broad-format file probably took its InChIKeys from PubChem using the same CIDs, so the CID result shows that the lookup is reproduced correctly, not that the vendor's CIDs are right. In LEOCC, about 5% of rows carry a CID and an InChIKey that describe different compounds; identifier resolution passes such errors on. Every derived key therefore keeps its source in the output.

The 9 HMDB disagreements include a dipeptide with the residues reversed (Ala-Phe against what looks like Phe-Ala), several lipids whose HMDB ID points to another species or isomer, and two small molecules whose HMDB ID points to a related compound. HMDB IDs in vendor files are less reliable than CIDs.

## 2026-10-06: Identifier resolution skips lipid names that are not one structure

First full run of v0.4 on the development data: the evaluation of the 53 reviewed groups did not change (34, 8, 1), but row flags rose from 15 to 93. Most new flags were on Broad rows such as "TG 42:0", "LPE 18:1" or "SM 34:2;O2". Broad gives no InChIKey for these Species-level names but does give an HMDB ID, and an HMDB record is always one specific structure. The derived key then claimed a structure the row never reported. Effects: 8 rows lost their lipid profile through lipid_source_conflict, and about 70 C04 flags compared two arbitrary representatives of the same Species.

Fix: resolve_identifiers() skips rows whose name rgoslin reads at Species, molecular species or sn-position level. Names that fix one structure, such as "LPE 16:0/0:0", are still resolved.

Found on development data, before v0.4 was frozen, so the COMETS test is not affected.

Second fix, same run: internal standards are also skipped. Broad "Inosine-15N4 [iSTD]" carries the HMDB ID of unlabelled inosine and got that compound's key. No wrong claim followed (R02 runs first), but the key described the wrong molecule.

## 2026-10-06: COMETS external test, criteria set before running

Data: the five COMETS files with expert mappings (ManualMappings.csv, MetLinkR SI), none used in development. MetLinkR run with default settings, RaMP database 3.0.7, metLinkR commit a6aa145 (metlinkr-run/run_metlinkr_comets.R). Global mapping rate 69.4%.

Tool: metMatchAudit v0.4.0 as tagged. scripts/60_comets_test.R stops if the package code differs from the tag. Identifiers: HMDB for inputs_1539354673, PubChem CID for inputs_1555609921; no file has InChIKey or SMILES.

Expert labels are matched to rows by file and metabolite ID. If the match is not one to one, the script stops; only the join may then be fixed, and the fix is logged before rerunning.

Only pairs of rows from different files count. Expert positive: same Manual_Metabolite_Group_Label.

Criteria:
- C1: among cross-file pairs the tool calls exact or distinct, at least 90% agree with the expert (exact with the same expert group, distinct with different groups).
- C2: removing the pairs the tool calls distinct from MetLinkR's links does not lower precision and does not lower F1, compared with MetLinkR alone.

Reported without a criterion: recall; the expert split of broader_narrower and unresolved pairs; the share of new lipid matches in one expert group; the share of expert cross-file groups found whole by MetLinkR (the paper reports 85.3%); Nightingale pairs on their own; key sources.

Run once. Every wrong definite claim is listed and reviewed afterwards, with the reason, which can be a tool error or an expert choice. Results go into findings.md and the README whatever they are.

## 2026-10-06: COMETS test runs on v0.4.1

The first attempt to run scripts/60_comets_test.R stopped at its version check, before reading any data: R/identifiers.R differed from tag v0.4.0. The difference is the internal standard fix (decisions.md, 2026-10-06, second fix), made on the development data and used in the last development evaluation, but copied in after the tag was set. The tag is left as it is. The fix is released as v0.4.1, and the test runs on that version. Criteria are unchanged.

## 2026-10-06: COMETS test, join fix

The first attempts of scripts/60_comets_test.R stopped before any comparison. ManualMappings.csv lists Nightingale Omega-6 six times and LA four times, each with a group, and its row order differs from the input files.

Change, made before any result was seen: rows are matched to expert labels by file and metabolite ID. A row can belong to several expert groups. Two rows count as the same for the expert when they share at least one group, and a pair in two shared groups counts once. The script also checks that audit rows are in input order before using IDs by position. It stops if any input row has no expert line, any expert line matches no row, or an ID repeats within a file.

Criteria and definitions unchanged.

The next attempt stopped too: 740 rows had no expert line. ManualMappings.csv writes IDs in upper case and, in inputs_1555609921, without spaces inside them. A check that printed only IDs and column names, not group labels, confirmed this. IDs on both sides are now compared in upper case with all spaces removed. A clash this creates within one file stops the test.

## 2026-10-07: comets-test-1 tag corrected

The tag was first created on 7aec708 by mistake and pushed. A minute later it was deleted and recreated on 9389f61, the commit the COMETS test ran on. Nothing referred to the wrong tag in between.

## 2026-10-07: v0.4.2 plan, written before any code

The COMETS test (`findings.md`, 2026-10-06) found four wrong definite claims. Three come from rules that trust identifiers when the names say otherwise. v0.4.2 changes only that. Name to RefMet stays in v0.5 and the aggregate relationship in v0.6.

The principle is the same as before: when the evidence does not settle a pair, the answer is unresolved. None of these changes adds a new way to call a pair distinct.

### Changes

1. **Names that contradict identical keys (R03c).** Runs just before R03. If the two InChIKeys are identical but the names state different isomers, the pair is unresolved, not exact. "State different isomers" means one of:
   - both names carry stereo descriptors and the sets differ. Descriptors are the existing ones (R/S, E/Z, cis/trans, alpha/beta, D/L) plus the cyclitol and sugar prefixes myo, scyllo, epi, chiro, neo, muco, allo, threo, erythro;
   - both names state an omega position (n3, n-3, omega-3, ω-3) and the positions differ.
   Evidence text says the identifier in one source may point to the wrong isomer. Fixes scyllo/myo-inositol and DPA n3/n6.

2. **Same name, different connectivity (R12n).** Runs just before R12. If the names are identical after normalization (lower case, spaces collapsed, a trailing `*` dropped) but block1 differs, the pair is unresolved. The keys conflict and the tool cannot tell which one is wrong. Isomer labels such as `(1)` and `(2)` stay part of the name, so DG 18:1_18:2 [1] and [2] still reach R12. Fixes fructose and 4-vinylphenol sulfate.

3. **Slash names as group names.** `is_group_name()` also matches two words of three or more letters joined by a slash, such as "mannitol/sorbitol" and "oleate/vaccenate (18:1)". Digits on either side of the slash do not count, so lipid names like "PC 16:0/18:1" are not affected.

### Acceptance, checked before the release

- Development evaluation (`scripts/21_evaluate_v0.1.R` on the v0.4.2 audit): at least 34 agree and at most 1 wrong, as in v0.4.
- Lipid match proposals: the same 63.
- Every row flag that changes is listed and explained.
- Unit tests for each change, including cases that must stay as they are: glycine against Gly stays exact, "PC 16:0/18:1" is not a group name, DG [1] against [2] stays distinct.
- COMETS is rerun and reported as after adaptation. No criterion is claimed as met from that rerun.

Correction, same day, before any code: change 3 is dropped. is_group_name() already treats two words joined by a slash as a group name. In COMETS, mannitol/sorbitol and oleate/vaccenate stopped at R01 because they had no key. Also, the existing stereo descriptors are R/S, E/Z, cis/trans, D/L and optical signs. alpha and beta are left out on purpose, as states_stereo() says. Expected side effect of R12n: the xylose pair in the development data becomes unresolved instead of distinct.

v0.4.2 acceptance, checked 2026-10-07: 286 tests pass. Development evaluation 35 agree, 8 safe_miss, 0 wrong, 10 E_not_exact (v0.4.1: 34, 8, 1, 10). Lipid proposals 63, row flags unchanged. Two pairs changed rule. xylose against Xylose went from R12 distinct to R12n unresolved, which agrees with its review category F. (R)- against (S)-3-hydroxybutyrylcarnitine went from R03 exact to R03c unresolved. That pair is outside the reviewed groups, so the wrong exact call in v0.4.1 had not been seen before.

## 2026-10-07: COMETS rerun with v0.4.2 (after adaptation)

Script `scripts/61_comets_after_v042.R`, output `results/comets_v0.4.2/`, same PubChem cache as the blind run. Not a blind test: v0.4.2 was written to fix the four wrong claims of the blind run.

- Exactly four pairs changed, the four targeted ones: fructose and 4-vinylphenol sulfate to R12n, scyllo/myo-inositol and DPA n3/n6 to R03c. All four are now unresolved.
- Definite claims: 121 of 121 agree with the expert.
- With the audit, precision and F1 equal MetLinkR alone (0.895). The audit now calls no COMETS pair distinct, so it removes nothing. C2 holds only in this trivial sense.
- 263 of 401 pairs are unresolved. None of the seven MetLinkR merges of different compounds is separated.
- On the development data the same changes also caught a wrong exact call that no review had seen ((R)- against (S)-3-hydroxybutyrylcarnitine). That is the only evidence so far that the new rules carry beyond the pairs they were written for. A second blind dataset is still needed.

## 2026-10-07: v0.5 feasibility plan, name to RefMet, written before any code

### Question

739 of 1429 COMETS rows have no structure key, because three files carry names only. Can a row's own name give it a key that is accurate enough to audit with?

### Source

RefMet from the Metabolomics Workbench. The REST service returns the whole table at `https://www.metabolomicsworkbench.org/rest/refmet/all/xx/all/txt`, with name, regno, pubchem_cid, inchi_key, exactmass, formula, sys_name, main_class, sub_class and synonyms. The table is downloaded once, its md5, download date and row count are logged here, and lookups run offline against that file, the same way SwissLipids is used. The `match` service is not used: it changes on the server without a version, so results could not be repeated.

### Independence

MetLinkR also uses RefMet. If MetLinkR built its standardized names from the same RefMet lookup, a key from RefMet would agree with MetLinkR's grouping by construction and would audit nothing. Before the feasibility run, read the MetLinkR code at commit a6aa145 and write down how it uses RefMet. The `Standardized name` column is never used as input to the lookup.

### Lookup, fixed before the run

1. Normalize the row's name with `normalize_name()`.
2. Look it up against RefMet names and synonyms, normalized the same way. Exact matches only.
3. Collect the InChIKeys of all hits. If they share one block1, the row gets that key, and the stereo block only if all hits agree on it (otherwise the key is kept at block1 and marked stereo unknown). If block1 differs between hits, the row gets no key and is counted as ambiguous.
4. Rows skipped by `resolve_identifiers()` are skipped here too: group names, alternative names, Species and molecular lipid names, internal standards.

### Data

- Accuracy set: every row in the development data (LEOCC and Broad) and in COMETS that already has a key, from the file or from identifier resolution. COMETS expert labels are not read.
- Coverage set: the 739 COMETS rows without a key. Only the count of rows that get a key is reported, no comparison.

### Criteria

- Accuracy: among accuracy-set rows that get a RefMet key, block1 agrees with the existing key in at least 95%. Same bar as identifier resolution.
- Coverage: at least 30% of the 739 keyless rows get a key. Below that, the route is not worth building and v0.5 goes another way.
- Both must hold. Every block1 disagreement is listed and read.

If either fails, the criteria stay as written and the result is logged as a failed feasibility test.

## 2026-10-07: v0.5 plan, rows checked against their RefMet entry, written before any code

### Why the first v0.5 plan was dropped

The plan to give keyless rows a key from their RefMet name was dropped before it was logged. MetLinkR's `Standardized name` is a RefMet name (`workflow.R` at commit a6aa145 queries RefMet name conversion first, then RaMP synonyms, then RefMet again). A key taken from RefMet would agree with MetLinkR's grouping by construction. The downloaded table also has no synonyms, so exact name lookup would have found very little.

### What v0.5 checks instead

MetLinkR claims, for each row, "this row is RefMet entry X". v0.5 tests that claim against the row's own evidence and reports conflicts as row flags. It does not change any pair result. A conflict on one row says nothing yet about the other row of a pair, so the pair rules stay as they are.

### Source

`data-raw/refmet/refmet_all.txt`, from `https://www.metabolomicsworkbench.org/rest/refmet/all/xx/all/txt`, downloaded 2026-10-07. md5 `63079d6af6aff5fec9ddc9d3f952afa2`, 208170 entries, 35618 with a valid InChIKey, no name used twice. A different file stops the run until the change is logged. Group names missing from the table are reported and left without a check.

### The check, fixed before the run

For each row whose group name is a RefMet entry, the first case that applies decides:

1. **Entry is a lipid at species level** (its name parses with rgoslin at species level, such as `CAR 5:0` or `MG 18:2`). Its InChIKey, when present, stands for one member of a sum composition and is not used.
   - If the row has a lipid profile, compare Species. Headgroup naming differences (HexCer against GlcCer) are not a conflict, as in R04h.
   - Else, if the row has a key, compare molecular formulas: the formula of the row's parent compound against the entry's formula.
2. **Entry is a lipid at molecular level or finer and the row has a lipid profile:** compare Species.
3. **Entry has an InChIKey and the row has a key:** compare block1 of the two parent compounds.
4. Otherwise there is nothing to compare.

Parent compounds come from PubChem, through the client and cache built for identifier resolution: the row's key and the entry's PubChem CID are each mapped to their parent compound (`cids_type=parent`), which removes salts and counterions. If several parents come back and they disagree on block1, or the lookup fails, the row is reported as unchecked, not flagged.

When the comparison finds a difference:

- If the row's name and the entry's name are the same after `normalize_name()`, the row gets **C07 `refmet_identifier_conflict`**. The names agree and the structures do not, so one of the identifiers is wrong or the two describe different forms (ring and open chain).
- Otherwise the row gets **C06 `refmet_assignment_conflict`**: the row's own evidence says it is not the compound MetLinkR assigned it to.

Group names, alternative names and internal standards are skipped, as in identifier resolution.

### Known limits, stated now

- block1 ignores stereo, so a row assigned to the wrong stereoisomer is not caught.
- Rows with no key and no lipid profile (most of the Biocrates, LC/MS and Nightingale files in COMETS) are not checked.
- Ring and open-chain forms with different names are flagged C06 even when they are the same sugar.

### Acceptance, on the development data

- Every C06 flag is reviewed by hand and labeled real or not real. Real means the row is not the entry's compound, or for a species-level entry, not within its sum composition. At least 80% must be real.
- These six rows, which the manual review placed in C or E groups, must get C06: `andro steroid monosulfate C19H28O6S (1)*`, `(15:0)-anacardic acid`, `docosapentaenoate (n3 DPA; 22:5n3)`, `N-acetyl-1-methylhistidine*`, `octadecadienedioate (C18:2-DC)*`, `pregnenediol sulfate (C21H34O5S)*`.
- Pair results and the development evaluation stay identical to v0.4.2, since no pair rule changes. All tests pass.
- COMETS is run and reported as after adaptation, with the number of C06 and C07 flags and the rows behind them. No expert label is used to tune anything.

If an acceptance check fails, it is logged as failed. The check is not tuned on the same data to make it pass.

## 2026-10-07: v0.5 acceptance failed

Run with scripts/71_v05_check.R. Pair results identical to v0.4.2 on both datasets, all 308 tests pass, and all six required rows got C06. The 80% check failed: of 40 C06 flags on the development data, 17 are rows MetLinkR assigned to the wrong compound, 22 are not, 1 is unclear (42.5%).

Why the 22 are not assignment errors:
- 15 come from the formula comparison. PubChem's parent compound for acylcarnitines and sphingomyelins is the cation (C17H32NO4+), while RefMet stores the neutral inner salt (C17H31NO4). Same compound, formulas one H and one charge apart.
- 2 are errors in RefMet itself: CAR 20:4 is stored as C31H54NO4 (exact mass 504.405, true formula C27H45NO4), and 5-hydroxy-6-methoxyduloxetine sulfate as C19H21NO7S, one sulfur short.
- 2 are probably wrong keys in the source file (oxindolylalanine, S-carboxyethylcysteine).
- 3 are one compound with two records or two tautomers (3-hydroxyhippurate sulfate, imidazole lactate, imidazole propionate).

The criteria stay as written and v0.5 is not released on this result. The check finds real conflicts between a row and its RefMet entry, but the flag name blames the assignment, and the cause can be the assignment, the source file or RefMet. Any change from here (neutral formulas, a flag that names the conflict but not its cause) is after adaptation on the development data and needs data not used so far for a real test.

## 2026-10-07: RefMet row check shelved

After the failed acceptance, the RefMet row check is set aside, not fixed. It only adds row flags and changes no pair result, so even a working version would not help the main COMETS problem: rows from name-only files carry no evidence. The code and the acceptance script are kept on branch shelved/refmet-check. main stays at v0.4.2. The two RefMet record errors found on the way (CAR 20:4 stored as C31H54NO4, 5-hydroxy-6-methoxyduloxetine sulfate one sulfur short) are kept in this log.

## 2026-10-07: v0.6 plan, Biocrates p180 names from the vendor's own table, written before any package code

### Why

In COMETS, the Biocrates file (inputs_1539353822) carries names only. Its rows have no key and no lipid profile, so every pair with them stops at R01. Most of the remaining COMETS problems sit in these pairs: isobaric sums such as C18:1 or SM C18:0, and wrong merges such as C14:2 with C14:1 or C5 with glutarylcarnitine.

### Source

"biocrates-p180 old vs. new nomenclature (v1-2023)", last modified 2023-02-24, `https://biocrates.com/wp-content/uploads/2023/12/p180-nomenclature.pdf`, md5 `4a1f84e43f5ad2f9dfe5aaa093fc0195`. Biocrates states the new names follow Liebisch et al. 2020 (LIPID MAPS shorthand). `scripts/80_biocrates_table.R` extracts the table from the PDF. Nothing is typed by hand. The script stops unless every section has the count the PDF states. It found 21 amino acids, 21 biogenic amines, 1 monosaccharide, 40 acylcarnitines, 90 glycerophospholipids and 15 sphingolipids, 188 in all. The PDF is not committed; the extracted table is, with its source.

### Mapping, fixed now

- Glycerophospholipids and sphingolipids: the new name from the table, as printed. This includes the vendor's own choices, such as `lysoPC a C28:0` to `PC O-28:0` and `SM (OH) C22:1` to `SM 41:1`.
- Acylcarnitines: the table keeps the old abbreviation and gives the chemical name, which fixes the chain. `Cx`, `Cx:y` and `Cx:y-OH` with x of at least 2 and no alternative become `CAR x:y`, plus `;OH` for the hydroxy forms. C0 (carnitine), the -DC and -M forms, and the four analytes listed with an alternative (C4-OH (C3-DC), C5-DC (C6-OH), C5-OH (C3-DC-M), C6 (C4:1-DC)) get nothing.
- Amino acids, biogenic amines and hexoses get nothing in v0.6. Their table entries give a name, not a structure.
- A row is matched by the abbreviation in the last pair of parentheses of its name, spaces removed, against the table abbreviation (with its alternative in parentheses, if any). Only analytes in the table are mapped. No rule extends the table to analytes it does not list.

### In the package, after the feasibility check passes

- A new optional `vendor` column in the files table. `"biocrates_p180"` is the only value accepted; any other value stops the run.
- Rows of such a file get a lipid profile from the mapped name through rgoslin, with `lipid_source` set to `"vendor_table"`. If the row already has a profile from another source and the Species differ, it gets C03 as today.
- No pair rule changes. The existing lipid rules (R04, R05, R05g, R06, R06b) then apply to these rows.

### Feasibility, `scripts/81_biocrates_feasibility.R`, counts only

- F1: every mapped name parses with an accepted rgoslin grammar.
- F2: at least 90% of the COMETS Biocrates acylcarnitine and lipid rows are found in the table. While testing the matching code, the names alone (no labels) gave 145 of 148. The three not found are C8:1, lysoPC a C6:0 and PC ae C40:0, which this release of the table does not list.
- F3: at least 20 cross-file pairs in which the Biocrates row would get a profile and the other row already has a lipid profile. Below that, the lipid rules would rarely fire and the work would change little. In that case the reason (likely Metabolon acylcarnitines without a lipid profile) is logged and v0.6 is reconsidered.

v0.6 feasibility, run 2026-10-07 with scripts/81_biocrates_feasibility.R. F1 met: 136 of 136 mapped names parse. F2 met: 145 of 148 COMETS Biocrates acylcarnitine and lipid rows found in the table (98.0%). F3 not met: 26 cross-file pairs have a mapped Biocrates row, and in only 2 does the other row have a lipid profile. The other 24 are mostly Metabolon acylcarnitines, whose names rgoslin does not read and whose keys SwissLipids does not resolve to a profile. The package change is not built on this result.

Also found: rgoslin reads one-chain names such as "CAR 14:2" at structure level and LPC names at molecular level. Biocrates measures by FIA-MS and cannot separate isomers, so a profile from the vendor table should be at species level whatever rgoslin reports. To be fixed in any next plan.

### Acceptance, after the package change

- All tests pass. Development evaluation identical to v0.4.2, since the development data has no Biocrates rows.
- COMETS rerun reported as after adaptation, with every pair whose rule changed listed and read.
- A blind test still needs another dataset with Biocrates p180 rows and expert or manual mappings.

scripts/80_biocrates_table.R could not run on the Mac: the pdftools binary failed to load (symbol _libintl_newlocale not found). It was run instead on Linux with pdftools 3.4.0 and poppler 24.02.0 on the same PDF (md5 4a1f84e43f5ad2f9dfe5aaa093fc0195). All six section counts matched the PDF, 188 analytes. Output data-raw/biocrates/p180_nomenclature.csv, md5 60efe0cbc17e40ff4370d4af2943cbd6, committed.

## 2026-10-07: plan, LIPID MAPS as a second lipid source, written before any package code

### Why

The v0.6 feasibility check failed on F3: in 24 of 26 COMETS pairs with a mapped Biocrates row, the other row had no lipid profile. Those rows are mostly Metabolon acylcarnitines. SwissLipids does not give them a profile. The LIPID MAPS Structure Database (LMSD) lists fatty acyl carnitines with InChIKeys and LIPID MAPS abbreviations, for example `CAR 14:1` with key `NNCBVXBBLABOCB-SEYXRHQNSA-N`.

### Source

LMSD standard SDF, dated 2026-03-06 on `https://www.lipidmaps.org/databases/lmsd/download`, CC BY 4.0. Zip md5 `8ee98567e23ffb7741ec5baeddb3ed84`. Unzipped `structures.sdf`, md5 `134c2b5c6d76e9883d46339a6289eefd`, 50586 records. Neither file is committed. Attribution to LIPID MAPS goes into the README and the package documentation.

### Lookup, fixed now

Same as SwissLipids: full key, then the first 25 characters, then the first 14, stopping at the first tier with a hit. The hits' abbreviations are parsed with rgoslin (accepted grammars only). The row gets a Species only if every hit gives the same Species. Group names, alternative names and internal standards are skipped.

### In the package, if the checks pass

- SwissLipids stays the first source. LMSD is used only for rows that get no profile from SwissLipids or from the name, so no result from v0.4.2 changes through LMSD alone. `lipid_source` gets the value `"lmsd"`.
- The Biocrates vendor table from the v0.6 plan comes with it. Profiles from the vendor table are set to species level whatever rgoslin reports, because Biocrates measures by FIA-MS and cannot separate isomers.

### Feasibility, `scripts/82_lmsd_feasibility.R`, counts only

- L1 (reported, no threshold): rows with a key that get a Species from LMSD and have none today, per dataset and file.
- L2: at least 20 COMETS cross-file pairs in which the Biocrates row gets a profile from the vendor table and the other row has a Species from SwissLipids or LMSD. Same bar as F3.
- L3: among rows where both SwissLipids (or the name) and LMSD give a Species, at least 95% agree, counting headgroup naming differences (HexCer against GlcCer) as agreement. At least 50 such rows are needed for the check to count. Every disagreement is listed and read.

If L2 or L3 fails, it is logged and neither LMSD nor the vendor table goes into the package on this plan.

LMSD feasibility, run 2026-10-07 with scripts/82_lmsd_feasibility.R. L1: 230 development rows and 154 COMETS rows get a Species from LMSD and have none today. L2 met: 22 of 26 COMETS pairs with a mapped Biocrates row now have a Species on the other row. L3 met: where SwissLipids or the name and LMSD both give a Species, they agree in 149 of 149 development rows and 66 of 66 COMETS rows.

Correction to the plan, before building: LMSD changes no row's existing profile, but rows that gain a profile can change pair results, since lipid rules run before R03. The statement that no v0.4.2 result changes was wrong. Acceptance is therefore: development evaluation not worse than v0.4.2 (35 agree, 0 wrong), every pair whose rule changes is listed and read, and every new lipid match proposal on the development data is reviewed by hand, as the 63 of v0.2 were.

v0.6 acceptance, checked 2026-10-07 with scripts/83_v06_check.R. Development data: no pair changed rule, evaluation identical to v0.4.2 (35 agree, 8 safe_miss, 0 wrong, 10 E_not_exact). LMSD gave 233 rows a profile. Lipid match proposals went from 63 to 110; all 47 new ones were reviewed by hand and are correct (notes/v06_dev_lipid_review.csv). 46 of them were already grouped by MetLinkR and only entered the lipid matches because the rows now have a profile; one (nonadecanoate against Nonadecylic acid) is a match MetLinkR missed.

COMETS, after adaptation: 23 pairs moved from R01 to a definite answer. Two became distinct (valerylcarnitine C5 against glutarylcarnitine, C14:1 against C14:2), both MetLinkR merges of different compounds that the expert also separates. 21 became broader_narrower (a Biocrates sum against one Metabolon structure); the expert groups 12 of them. 36 new lipid matches across files, 18 in the same expert group. The other 18 are all Biocrates sums matched one-to-many to single structures (C4, C5, sn-2 LPC, diacyl PC), chemically consistent but against the expert's one-to-one practice.

## 2026-10-07: second blind test, plan, written before MetLinkR is run on the data

### Data

- A: Metabolon annotation of the author's own mouse dataset, 494 rows. Columns copied by hand from the original Metabolon file into a new workbook, no values and no edits. Eight CAS numbers were turned into dates by Excel on the way; CAS is used neither by MetLinkR nor by this tool. Private: the file stays in `data-raw/own/` and is never committed.
- B: annotation from Ding J. et al. A metabolome atlas of the aging mouse brain. Nat. Commun. 2021, 12. DOI 10.1038/s41467-021-26310-y. 1547 rows with name, InChIKey, PubChem CID and SMILES, copied the same way from the paper's supplementary data. Source file and sheet: https://pubmed.ncbi.nlm.nih.gov/34654818/.
- md5 of both workbooks: MD5 (data-raw/own/2ndData.xlsx) = 1b4876e548f55567bfcdd33f494c978a
MD5 (data-raw/own/MyData.xlsx) = b2d29cb02ffc65cc53fd5b08e4b85805.

The tool has never seen B's naming style. B's names follow MS-DIAL-like notation ("PE 38:1; PE 18:0_20:1;", "ACar 12:0;"), 67 of its PubChem CIDs carry invisible characters, and 15 InChIKeys start with "null_". None of this is cleaned. The test is meant to show how the tool does on data as it comes.

### Frozen code

metMatchAudit at tag `v0.6.0`, with LMSD (md5 `134c2b5c6d76e9883d46339a6289eefd`) and SwissLipids (md5 `b6ddd5bc3003967dce95e905e5748953`). MetLinkR at commit a6aa145 with RaMP DB 3.0.7, as for COMETS. The test script stops if the package code differs from the tag.

### Order, each step committed before the next

1. This plan, with the two fill-ins done.
2. `scripts/91_blind2_prepare.R` reads both workbooks with every column as text and writes CSV copies, logging rows and md5.
3. MetLinkR on the two CSV files. Only MetLinkR's own output is looked at: it is what is being audited.
4. `scripts/92_blind2_sample.R`, seed 20261007:
   - 60 groups drawn at random from MetLinkR groups that hold rows from both files (all of them if there are fewer);
   - 40 rows drawn at random from B whose name starts with a lipid class: `^(L?P[ACEGIS]|SM|Cer|HexCer|Hex2Cer|GlcCer|TAG|TG|DAG|DG|MAG|MG|CE|ACar|CAR|FA)[ -]`.
5. Labels by the author, with Claude's help, before the tool is run on these files:
   - every pair in the 60 groups gets a category from `notes/review_guide.md` (A exact, B stereo, C different compounds, D broader or narrower, E source error, F unclear) and the record behind the decision;
   - every one of the 40 lipid rows gets the rows of A that are the same compound or a broader or narrower annotation of it, or "none".
   The labels answer "same compound", the question the tool asks, not "same measurement". Labels are committed before step 6.
6. `scripts/93_blind2_test.R` runs the audit once and compares.

### Criteria

- C1, pairs: of the tool's definite claims (exact, distinct) on labeled pairs, at least 90% agree with the label. exact agrees with A; distinct agrees with C or E; exact on an E pair is wrong.
- C2, groups: each group's most severe relationship against its most severe label, as in `scripts/21_evaluate_v0.1.R`. Among groups labeled A to D, at most 5% get a wrong definite answer. Agree and unresolved are reported.
- C3, lipid matches: of the tool's lipid match proposals that involve a sampled B row, at least 90% are pairs in the labels. Recall against the labels is reported, with no threshold.

All three are reported whatever the outcome. No rule changes between the labels and the run, and none are made on this data afterwards without being marked as after adaptation.

### Known weakness

The labeler wrote the tool. To limit that: labels come before any tool output on these files, every decision cites its record, and if a colleague can label 20 of the groups independently, agreement between the two is reported.

Two changes to the second blind test plan, before any step ran: B's column "Pubchem CID" is renamed to PUBCHEM in the CSV copy so MetLinkR reads it without a space; no value changes. The labeling sheets hold names from the private file A, so they stay in data-raw/own/blind2/ and their md5 is committed instead, once they are filled in. This fixes the labels in git before the audit runs without publishing them.

## 2026-10-07: second blind test, result and what follows

C1 and C2 were not met. C3 was met, but on one proposal. Details and causes are in `findings.md` (2026-10-07, second blind test). The result is reported as it came out, in the README too.

- No rule is changed on the basis of this run without the label "after adaptation". A rerun of these files after any change is not a blind result.
- The labels stay as committed. The labels and the tool disagree in two places about what a pair is: name against key (causes 2 and 3), and the parent's key on "isomer" features (cause 1). Both are noted in the findings and not resolved by relabeling.
- Next version (v0.7): candidates are listed at the end of the findings entry. Before building, each gets a written plan and an acceptance check on the development data, COMETS and these files, all three reported as after adaptation.
- Any further claim of blind performance needs a third test on data not used so far, with criteria committed first.

## 2026-10-07: local copy of the sealed development audit overwritten

`git status` before the second blind test showed `results/audit_v0.1_sealed.rds` modified. The files in `results/audit_v0.1/` and `results/evaluation_v0.1/` were also modified. The local copy of the rds is a v0.4.2 run from 2026-10-07 09:28. The committed copy is a v0.4.0 run from 2026-10-06, and git history keeps every committed version. Most likely `scripts/20_audit_example.R` and `scripts/21_evaluate_v0.1.R` were run again with v0.4.2. Both write to these fixed paths.

`scripts/83_v06_check.R` (run 12:45) read the local copy as its v0.4.2 baseline. Its message says so, and the copy is indeed v0.4.2. So the v0.6 acceptance comparison stands.

Fix:
- the local v0.4.2 outputs move to `results/dev_v0.4.2/`;
- the committed files are restored;
- `scripts/83_v06_check.R` now reads its baseline from `results/dev_v0.4.2/audit.rds`.

Nothing is rerun. From now on, scripts that rerun an old step write to a folder named after the version.

## 2026-10-07: v0.7 plan, written before any code

v0.7 addresses three of the four causes found in the second blind test (`findings.md`, 2026-10-07). Every result on the development data, COMETS and the second blind test files after these changes is after adaptation. A blind claim needs a third test.

### Changes

1. **Isomer marker in names (new rule R03i).** A name that contains the word "isomer" or "isomers", alone or followed by a letter or number ("Carnosine isomer A", "Hexose phosphate isomer D"), marks a feature whose identity is not known. Its identifiers may belong to the parent compound, as in B. When either row of a pair carries the marker, no claim is made from keys or from R06b, and the pair is unresolved. The rule sits right after R03c and before R06b, R03 and all later key rules. "stereoisomer" and other words that only contain "isomer" do not count. The lipid rules R05 and R06 still apply, because they rest on the name, not the key.
2. **MS-DIAL lipid names.** Before rgoslin reads a name, the name is cleaned:
   - a trailing ";" is dropped;
   - a name with several annotations separated by "; " ("PE 38:6; PE 18:2_20:4") is read from its most specific annotation;
   - a name with " or " between annotations is treated as alternatives and gets no profile, as alternative names already do;
   - "ACar" is read as "CAR" if rgoslin does not know it.

   No new level rule is needed. `lipid_evidence()` already lets the name set the level when name and key agree on the Species. The B rows got a structure level only because their names did not parse.
3. **Optical rotation prefixes.** `normalize_name()` drops a leading "(+)-", "(-)-" or "(±)-". Two rows whose names then match and whose keys differ in connectivity go to R12n (unresolved) instead of R12 (distinct). This is the N-acetylneuraminic acid case.

Not in v0.7: cause 3, stereo given in only one record. The tool compares records, and calling these pairs broader_narrower is by design. The labels compared compounds by name. This stays a documented difference. The inconsistent S1P label is noted in findings and is not changed.

### Order, each step committed before the next

1. This plan.
2. Feasibility, `scripts/95_v07_feasibility.R`, without changing the package:
   - count names with the isomer marker in the development data, COMETS, A and B, and list them;
   - run the name cleaning as a stand-alone function on B's 960 lipid rows, and report how many get a profile from the name and which do not;
   - count names with an optical rotation prefix in all four sets.

   If fewer than 80% of B's lipid rows get a name profile after cleaning, change 2 is dropped from v0.7 and the reason is logged.
3. Build, with unit tests for the marker, the cleaning and the prefix.
4. Acceptance, `scripts/96_v07_check.R`.

### Acceptance

- **Development data:**
  - the evaluation against `notes/review.csv` is no worse than v0.6.0: 35 agree, 8 safe miss, 0 wrong, 10 E not exact;
  - every pair whose rule changes is listed and checked by hand;
  - any of the 110 lipid match proposals that is lost is listed with its reason;
  - new proposals are checked by hand, and all must be correct.
- **COMETS, after adaptation:**
  - no definite claim that agreed with the expert may change to one that does not;
  - changed pairs are listed with the expert's view;
  - precision and F1 are reported.
- **Second blind test files, after adaptation:**
  - scored with the same code as `scripts/93_blind2_test.R`;
  - reported, with no threshold, since the changes come from these files.

  Expected beforehand:
  - the 11 wrong `exact` calls on isomer features become unresolved;
  - the three acylcarnitine pairs become broader_narrower through R05, if B's names parse;
  - the two Neu5Ac pairs become unresolved.

  Anything else that changes is listed.

If any development or COMETS condition fails, v0.7 is not tagged until the cause is found and logged.

### After v0.7

A third blind test on data not used so far, with criteria C1 to C3 as in the second test, written and committed first. Ideally a file from a vendor or a naming style the tool has not seen, labeled before the run.

## 2026-10-07: v0.7 feasibility result and changes to the plan

The feasibility script is `scripts/95_v07_feasibility.R`, with its log in `results/v07_feasibility_log.txt`. Its results:

- **Isomer marker.** 55 rows of B carry the word "isomer", 53 of them with a key. There are none in the development data, COMETS or A, so the rule cannot change any result there.
- **MS-DIAL names.** After cleaning, 598 of B's 960 lipid rows parse from their name, which is 74.9% of the 798 that could. That is below the 80% set in the plan, so change 2 is dropped from v0.7, as the plan says.

  The names that fail fall into three groups:
  - MS-DIAL class names that rgoslin does not read (Cer-NS, Cer-NDS, HexCer-NS and the like);
  - fatty acids written with a common name after a semicolon ("FA 16:1; (palmitoleic acid)");
  - names whose last part is a lone letter ("CE 20:5; A").

  Reading these would need a class map for MS-DIAL notation. That is a plan of its own for a later version, not a fix to make here.

  Of the names that parse and also have a profile from their key, 8 disagree on the Species:
  - five are GlcCer against HexCer, a naming difference;
  - LPE 21:0 against PE 21:0;
  - PC 35:3 against PE 38:3, the same formula with a different head group;
  - PC 36:3 against PC O-36:3.

  The last three are cases where B's key and B's name describe different lipids. Name against key is a check the tool could make for lipids, as C04 does for structures.
- **Rotation prefix.** One row, B:1 "(-)-N-Acetylneuraminic acid". The plan expected both Neu5Ac pairs to become unresolved. That holds only for B:1 against B:627 "N-Acetylneuraminic acid". Against A's "N-acetylneuraminate" the names still differ after the prefix is dropped (acid against -ate), so that pair stays distinct. The expectation in the plan stays as written, and this is what the check should show.

Two changes to the plan, made before the build:

1. **The marker also covers MS-DIAL's lone letter.** A name whose last ";"-separated part is a single capital letter ("Cer-NS d42:2; Cer-NS d18:1/24:1; A") marks the same thing as "isomer A". B uses both forms, and the feasibility run found the letter form among the names that fail to parse.
2. **R03i sits right after R01, before the lipid rules, not after R03c.** The plan said the lipid rules rest on the name and can still apply. In B that is not true: the marked lipid rows got their profile from the key, the same key the marker says not to trust. So a marked row decides nothing, whether through its key or through a profile found by that key. R02 (internal standard) and R01 (no evidence) still come first.

v0.7 is now R03i and the rotation prefix. Acceptance is unchanged. The word marker and the prefix touch no row in the development data or COMETS. The letter form was not counted there, so `scripts/96_v07_check.R` reports the marker count for every set. If those counts are zero, both sets should come out unchanged, and any difference is a bug.

## 2026-10-07: COMETS numbers for v0.6 corrected

The README gave 123 of 123 definite claims for v0.6, with precision 0.900 and F1 0.897. `scripts/96_v07_check.R` counted 122 of 122. `scripts/97_comets_recount.R` reruns the evaluation code of `scripts/60_comets_test.R` on the stored audits of every version. It reproduces the blind v0.4.1 result exactly, 121 of 125 with F1 0.893 against 0.895.

| | v0.4.1 | v0.4.2 | v0.6.0 | v0.7.0 |
|---|---|---|---|---|
| definite claims right | 121 of 125 | 121 of 121 | 122 of 122 | 122 of 122 |
| precision with the audit | 0.895 | 0.895 | 0.897 | 0.897 |
| F1 with the audit | 0.893 | 0.895 | 0.896 | 0.896 |

MetLinkR alone has precision 0.895 and F1 0.895 throughout.

The earlier numbers were worked out by hand from the list of changed pairs and never counted. That list covered all pairs, and one of the two new `distinct` calls is a pair within one file, which the COMETS evaluation does not count. The README now gives the recounted numbers. From now on, every number in the README comes from a script.

## 2026-10-07: third blind test, plan, written before MetLinkR is run on the data

### Data

- **A**: the same private Metabolon file as in the second test, `data-raw/own/A_metabolon.csv`, 494 rows. Step 2 logs its md5 next to the md5 from the second test, and the test script checks it. A's naming is no longer new to the tool. It is here because it is a real vendor file the new files can be matched against.
- **C1** and **C2**: Supplementary Data 2 (feces, 533 rows) and 3 (serum, 231 rows) of Lai Y. et al. High-coverage metabolomics uncovers microbiota-driven biochemical landscape of interorgan transport and gut-brain communication in mice. Nat. Commun. 2021, 12, 6000. DOI 10.1038/s41467-021-26209-8. Open access, CC BY 4.0. Workbook md5: C1 `24b6e7dcf7e99befe9a49aa222dca921`, C2 `22ca253e0e86ffdcb6fb6046f3794667`.
  - Every row has a name and a PubChem CID.
  - The names were annotated by the authors at MSI level 1 or 2 and curated by hand. They follow HMDB and LIPID MAPS style ("PI(16:0/18:1(9Z))", "LysoPE(0:0/14:1(9Z))", "Isovaleryl-L-carnitine", "D-(+)-Pantothenate", "Sucrose (or Cellobiose)"), a style the tool has not been tested on.
  - The lists hold only metabolites that differ between germ-free and conventional mice. That is a selection, but it does not bear on how names match.
- Only the name and the PubChem CID are passed on. Name becomes `BIOCHEMICAL` and PubChem_CID becomes `PUBCHEM`, as for B. The workbooks and the CSV copies are public and go in `data-raw/lai2021/`.

What has been seen of C1 and C2 before this plan: the column headers, about 60 names printed while choosing the data, and the 40 lipid and carnitine names listed in the step below. No tool output on these files exists.

### Frozen code

metMatchAudit at tag `v0.7.0`, with SwissLipids (md5 `b6ddd5bc3003967dce95e905e5748953`) and LMSD (md5 `134c2b5c6d76e9883d46339a6289eefd`). MetLinkR at commit a6aa145 with RaMP DB 3.0.7. The test script stops if the package code differs from the tag.

### Order, each step committed before the next

1. This plan.
2. `scripts/101_blind3_prepare.R`: CSV copies of C1 and C2, every column as text, with row counts and md5 logged.
3. MetLinkR on A, C1 and C2 together. The output goes to `data-raw/own/metlinkr_output_blind3/`, so the second test's output is not touched.
4. `scripts/102_blind3_sample.R`, seed 20261008:
   - 60 groups drawn at random from MetLinkR groups that hold rows from at least two of the three files (all of them if there are fewer);
   - every row of C1 and C2 whose name matches `^((Lyso)?(PC|PE|PS|PI|PG|PA)|SM|Cer|CerP|GlcCer|GalCer|LacCer|TG|DG|MG|CE)\(|carnitine$` (case ignored). That is 40 rows, so all are taken.
5. Labels, made as in the second test, before the tool is run on these files:
   - every pair in the 60 groups gets a category (A exact, B stereo, C different compounds, D broader or narrower, E source error, F unclear) and the record behind it;
   - every one of the 40 rows gets the rows of the other two files that are the same compound or a broader or narrower annotation of it, or "none".

   The md5s of the label sheets go into `notes/blind3_labels.md`, which is committed before step 6.
6. `scripts/103_blind3_test.R` runs the audit once and compares.

### Criteria

The same as in the second test:

- **C1:** of the tool's definite claims (exact, distinct) on labeled pairs, at least 90% agree with the label. exact agrees with A; distinct agrees with C or E.
- **C2:** among groups labeled A to D, at most 5% get a wrong definite answer. A group's label is its most severe pair label, in the order E > C > F > D > B > A.
- **C3:** of the tool's lipid match proposals that involve a sampled row, at least 90% are pairs in the labels. Recall is reported, with no threshold.

All three are reported whatever the outcome.

If fewer than 30 groups hold rows from at least two files, the test still runs, and the report says the sample is too small to carry the thresholds.

### Known weaknesses

- The labeler wrote the tool, as before. No independent labeler is available.
- A was used in the second test, so its naming is known. The new part is C1 and C2, and every result is also reported for pairs that do not involve A.

## 2026-10-07: v0.8 plan, written before any code

All 7 wrong groups in the third blind test, and part of cause 3 in the second, come from one situation:
- two rows name the same compound;
- neither name says anything about stereo;
- one lab's PubChem record has stereo and the other's does not.

R09 then calls the pair broader_narrower. The difference lies in which record each lab picked, not in what was measured. Every result on data already used after this change is after adaptation.

### Change

R09 and R09b make a claim only when at least one of the two names states stereo (`name_states_stereo`). When neither does, the pair goes to a new rule **R09n**, unresolved, placed right before R09. This mirrors how R10 and R11 already treat two specified keys: a claim when a name states stereo, unresolved when only the keys differ.

The change can only turn a broader_narrower answer into unresolved. It cannot add a wrong definite answer.

### Order, each step committed before the next

1. This plan.
2. Feasibility, `scripts/104_v08_feasibility.R`, without changing the package. On the v0.7.0 audits of all four sets (development, COMETS, second and third blind test files), it lists every R09 and R09b pair where neither name states stereo, with its label:
   - development: the review category of the group;
   - COMETS: whether the expert groups the pair;
   - blind tests: the pair label.

   **Go-ahead:** among the touched pairs that carry a pair or group label (development, second and third test), at least two thirds must have a label other than D. A label other than D means broader_narrower was not the curator's answer there. If fewer do, the change is dropped and the reason is logged.
3. Build, with unit tests.
4. Acceptance, `scripts/105_v08_check.R`.

### Acceptance

- **Development:**
  - wrong stays at 0;
  - every changed pair is listed;
  - agree turning into safe miss is reported and allowed;
  - lipid match proposals lost are listed with the reason.
- **COMETS:**
  - definite claims unchanged (122 of 122);
  - precision and F1 recomputed with `scripts/97_comets_recount.R`. They should not move, since only `distinct` pairs are removed from MetLinkR's matches.
- **Second and third test files:** scored with the same code as their test scripts and reported as after adaptation, with no threshold.

  Expected beforehand:
  - third test: C2 goes from 7 to 2 wrong groups (CAR 16:0 and malate stay, since one name there says L);
  - second test: the groups of cause 3 where neither name states stereo leave the wrong column.

  Anything else that changes is listed.

Lipid proposals: the change also removes proposals that rest on R09 alone. Third-test C3 recall will fall. That is reported and not counted against acceptance, since a pair the tool cannot settle should not be proposed as a match.

## 2026-10-07: v0.8 feasibility result, R09n dropped

`scripts/104_v08_feasibility.R` (log in `results/v08_feasibility_log.txt`) counted the R09 and R09b pairs where neither name states stereo:

| set | touched | labeled | label other than D |
|---|---|---|---|
| development | 27 | 25 | 1 (an E pair) |
| COMETS | 16 | | 15 grouped by the expert |
| second test | 25 | 5 | 3 |
| third test | 23 | 9 | 9 |

Only 13 of the 39 labeled pairs, 33.3%, have a label other than D. The go-ahead needed two thirds, so **R09n is dropped**, as the plan says.

The count shows something more important than the rule. My own reviews used two different conventions for the same situation:
- The development review (2026-10, v0.1) labeled pairs such as lactate against lactic acid, malate against malic acid and cholate against cholic acid as D, broader or narrower. One record has stereo and the other does not, and the review followed the records.
- The blind labels of the second and third tests followed the names instead. When neither name states stereo, the label is A.

The tool was built and checked against the first convention. Its R09 calls agree with it in 24 of 24 development pairs. The third test's C2 failure therefore comes from the labeling convention fixed for that test, which differs from the one the tool was evaluated against. This does not change the result: C2 was not met and stays not met. It explains the result.

Decision:
- No rule change. The two conventions answer different questions. Which record is more specific is the development review's question. Whether a curator merging by name would call two rows the same compound is the blind labels' question.
- The README will say plainly that broader_narrower from R09 means the two records differ only in whether stereo is given, and that a user who matches by name can treat such a pair as the same compound.
- Which convention the tool should follow is a design decision. It has to be made explicitly, before any new rule. The go-ahead check above mixed labels from both conventions, which is why it could not settle the question. It is not rerun with another threshold.
- Any further test fixes in advance which convention its labels follow.

## 2026-10-08: README shortened

The README had grown to a full report of every test. It now gives one summary table and the main limits. The sections on each test, what did not work and the open problems moved unchanged to `notes/evaluation.md`. The limits left out of the README were added to `design.md`. No number changed.

## 2026-11-07: README graphics

The README now opens with a header image in the style of my GitHub profile, and has a diagram of how a pair gets its answer and a chart of the blind test results. The images are in `man/figures/`. The chart shows the blind results only, with criteria that were not met. The tables stay, so every number can still be read as text. No number changed.

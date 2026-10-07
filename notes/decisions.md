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

Correction, same day, before any code: change 3 is dropped. is_group_name() already treats two words joined by a slash as a group name. In COMETS, mannitol/sorbitol and oleate/vaccenate stopped at R01 because they had no key. Also, the existing stereo descriptors are R/S, E/Z, cis/trans, D/L and optical signs. alpha and beta are left out on purpose, as states_stereo() says. Expected side effect of R12n: the xylose pair in the development data becomes unresolved instead of distinct.

v0.4.2 acceptance, checked 2026-10-07: 286 tests pass. Development evaluation 35 agree, 8 safe_miss, 0 wrong, 10 E_not_exact (v0.4.1: 34, 8, 1, 10). Lipid proposals 63, row flags unchanged. Two pairs changed rule. xylose against Xylose went from R12 distinct to R12n unresolved, which agrees with its review category F. (R)- against (S)-3-hydroxybutyrylcarnitine went from R03 exact to R03c unresolved. That pair is outside the reviewed groups, so the wrong exact call in v0.4.1 had not been seen before.


## 2026-10-07: v0.4.2 plan, written before any code

The COMETS test (`findings.md`, 2026-10-06) found four wrong definite claims. Three come from rules that trust identifiers when the names say otherwise. v0.4.2 changes only that. Name to RefMet stays in v0.5 and the aggregate relationship in v0.6.

The principle is the same as before: when the evidence does not settle a pair, the answer is unresolved. None of these changes adds a new way to call a pair distinct.

## 2026-10-07: COMETS rerun with v0.4.2 (after adaptation)

Script `scripts/61_comets_after_v042.R`, output `results/comets_v0.4.2/`, same PubChem cache as the blind run. Not a blind test: v0.4.2 was written to fix the four wrong claims of the blind run.

- Exactly four pairs changed, the four targeted ones: fructose and 4-vinylphenol sulfate to R12n, scyllo/myo-inositol and DPA n3/n6 to R03c. All four are now unresolved.
- Definite claims: 121 of 121 agree with the expert.
- With the audit, precision and F1 equal MetLinkR alone (0.895). The audit now calls no COMETS pair distinct, so it removes nothing. C2 holds only in this trivial sense.
- 263 of 401 pairs are unresolved. None of the seven MetLinkR merges of different compounds is separated.
- On the development data the same changes also caught a wrong exact call that no review had seen ((R)- against (S)-3-hydroxybutyrylcarnitine). That is the only evidence so far that the new rules carry beyond the pairs they were written for. A second blind dataset is still needed.

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

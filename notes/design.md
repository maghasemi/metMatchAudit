# Design

This file says what the tool does before the code does it. If the code needs to behave differently, change this file first and log the change in `notes/decisions.md`.

Status: revised 2026-09-29, after the lipid feasibility test (`notes/findings.md`).

One principle runs through every rule below. When the evidence does not settle a question, the answer is `unresolved`. A rule that fixes one problem by guessing, and so creates a new one, does not go in.

## What the tool answers

Given rows from two or more datasets that some matching step put into the same group, the tool says how each pair of rows relates. MetLinkR is the first matching step we support, through an adapter. The core does not depend on it.

Answers are per pair of rows, not per group. One group can hold a correct match, a wrong merge and a source error at the same time (pregnenolone sulfate and SM 36:1 groups in the review).

A second capability comes in v0.2: for lipid rows that the matching step left out, propose matches at the finest level both rows support. This is the problem that started the project, validating Metabolon data against another vendor's lipid names.

## Input

One line per row:

| column | required | notes |
|---|---|---|
| row_id | yes | unique, e.g. `LEOCC_Metabolon:132` |
| file | yes | source dataset |
| original_name | yes | as written in the source file |
| group_id | yes | the group the matching step assigned; for MetLinkR, `Standardized name` |
| inchikey | no | from the source file, never edited |
| smiles | no | from the source file, used only for a consistency check |

The MetLinkR adapter builds this table from the `*_metLinkR.csv` files. `scripts/01_build_rows.R` is its first version.

## Evidence per row

Nothing is typed in by hand. Each value comes from the input file, from a fixed parsing rule, or from a versioned database file.

### InChIKey

- valid or not
- block1 (connectivity), block2, protonation letter
- stereo unspecified: block2 is `UHFFFAOYSA`

### Name checks

Regular expressions, one function each, each with tests.

- group name: " or ", "and or", "and/or", " + ", or a list inside parentheses such as "(d18:1/20:2, d18:2/20:1)"
- tentative: trailing `*`
- isomer label: `(1)`, `(2)`, `[1]`, `[2]`
- internal standard: `[iSTD]`, or an isotope label such as `15N4`, `13C6`, `d3`
- stereo stated: R/S, E/Z, cis/trans, alpha/beta, D-/L- tokens

### Lipid evidence

Two sources, both tested in the feasibility run.

**SwissLipids, by InChIKey.** For rows with a valid key, no group name and no alternative names.

1. Look the key up in three tiers: full key, first 25 characters (SwissLipids stores pH 7.3 forms, so the protonation letter often differs), first 14 characters. Stop at the first tier that finds anything. Record the tier.
2. Walk up the parent links from each hit and collect the Species and Molecular subspecies ancestors.
3. If no hit reaches a Species, collect every entry with the same first 14 key characters and use their Species, but only if they all agree on exactly one. Same connectivity means same chain composition, so a disagreement here points to a database problem, and the row stays without a Species.
4. Parse the collected SwissLipids abbreviations with rgoslin, so they come out in the same notation as shorthand names from other vendors.

The SwissLipids file is identified by its md5 (`b6ddd5bc3003967dce95e905e5748953` for the release used so far). A different file stops the run until the change is logged.

**rgoslin, by name.** For every row. A parse counts only if rgoslin reports a grammar other than `NOT_PARSEABLE` and a level finer than class. Each name is parsed on its own, and rgoslin's text value "NA" is read as missing.

From both sources the row gets a lipid profile:

- `lipid_species`: e.g. `SM 36:1;O2`
- `lipid_molecular`: e.g. `SM 18:1;O2/18:0`, if known
- `lipid_level`: the finest level the evidence supports
- `lipid_source`: `swisslipids`, `swisslipids_fallback`, or `rgoslin_name`

When both sources give a Species for the same row and they differ, the row gets a `lipid_source_conflict` flag and no lipid profile.

## Pair rules

Rules run in this order for each pair in a group. The first rule that applies decides. Every result keeps its rule ID.

| ID | condition | relationship |
|---|---|---|
| R01 | either row has an invalid InChIKey and no lipid profile | unresolved |
| R02 | either row is an internal standard | distinct |
| R03 | both InChIKeys identical | exact |
| R04 | both rows have a lipid profile and `lipid_species` differs | distinct |
| R05 | same `lipid_species`, rows reported at different levels | broader_narrower |
| R06 | same `lipid_species`, both at molecular level or finer, `lipid_molecular` differs | distinct |
| R07 | both formulas known and different (v0.3) | distinct |
| R08 | either name is a group name | group_overlap |
| R09 | same block1, exactly one key has unspecified stereo | broader_narrower |
| R10 | same block1, both keys give stereo, at least one name states stereo | stereo_variant |
| R11 | same block1, both keys give stereo, neither name states stereo | unresolved |
| R12 | different block1 | distinct |
| R99 | nothing above applies | unresolved |

Where the rules come from:

- R02: inosine-15N4 [iSTD] merged with inosine.
- R04: SM 36:1 merged with an SM 38:3 group.
- R05: Broad "SM 36:1;O2" against Metabolon "stearoyl sphingomyelin (d18:1/18:0)". Both resolve to Species `SM 36:1;O2`.
- R06 and R12: DG 18:1_18:2 [1] and [2] share Species `DG 36:3` but differ in structure (1,2- against 1,3-DG). R06 does not separate them because their molecular species is the same, so R12 does, through the InChIKey.
- R07 before R08: androsterone sulfate against "andro steroid monosulfate C19H28O6S". A vague name does not make a pair broader/narrower when the formulas differ.
- R08: "(2 or 3)-decenoate", "2-hydroxyvaleric and or 2-hydroxyisovaleric acid". The key of a group row is one member, not a confirmed structure.
| R08b | either name lists alternative names at its top level ("X or Y") | unresolved |
- R11: 2-hydroxystearate as 2S in one file and 2R in the other, with no stereo in either name.

A note on headgroups: SwissLipids writes glycosylceramides at Species level as HexCer or Hex2Cer, while a name may say GlcCer or LacCer. v0.1 does not treat HexCer as a parent of GlcCer. Such pairs differ in `lipid_species`, which would make R04 call them distinct. To avoid that false result, R04 skips a pair when the two Species strings differ only in the headgroup and the chain part is identical. The pair falls through, usually to `unresolved`.

## Row checks

| ID | condition | flag |
|---|---|---|
| C01 | InChIKey fails the format check | invalid_key |
| C02 | the key has stereo but the name states none | key_more_specific_than_name |
| C03 | SwissLipids and rgoslin give different Species for the row | lipid_source_conflict |
| C04 | the row's SMILES and InChIKey describe different structures (v0.3) | smiles_key_mismatch |
| C05 | the row's name is not among the PubChem synonyms of its key (v0.3) | name_key_mismatch |

C04 catches xylose in the Broad file, where the SMILES is the ring form and the InChIKey the open-chain form. C05 catches pregnenediol sulfate and the 2-hydroxycaproic acid row. Both need a chemistry toolkit or PubChem, so they wait.

## Output

1. `pairs`: row_a, row_b, group_id, relationship, rule_id, evidence text
2. `row_flags`: row_id, flag
3. `row_evidence`: row_id with its InChIKey parts, name checks and lipid profile
4. `groups`: group_id, number of rows, most severe relationship

Severity, most severe first: source_conflict, distinct, unresolved, group_overlap, broader_narrower, stereo_variant, exact.

`merge_level` sets which relationships count as "same" downstream. Default: `exact` only.

## Known limits

- Ring and open-chain forms of sugars (xylose) have different block1 values. R12 calls them distinct. C04 will flag the source row in v0.3, but the pair result stays wrong. Reported in the evaluation.
- Some lipids have no usable SwissLipids entry: d18:2 and d17:1 sphingoid bases and one hydroxylated SM in the test data. They get no lipid profile unless their name parses with rgoslin.
- Glycosylceramide pairs that differ only in headgroup naming stay unresolved (see the note under the pair rules).
- Group-name rows get no lipid profile in v0.1, even though 20 of 29 in the test data map to one Species through their representative key. Using that needs proof that all members of a group share one chain sum.
- Name regexes are tuned on Metabolon and Broad naming.

## Versions

- v0.1: audit of existing groups. InChIKey, name checks, SwissLipids file, rgoslin. No network access during a run. Rules R01 to R06 and R08 to R12, row checks C01 to C03.
- v0.2: lipid match proposals for rows the matching step left out, using the same lipid profiles. Output: candidate pairs with the level they agree at.
- v0.3: PubChem lookup with cache and a chemistry toolkit. Adds R07, C04, C05.

## Evaluation

Only after the manual review is finished and committed. Review categories map to relationships:

| review category | expected relationship |
|---|---|
| A | exact |
| B | stereo_variant |
| C | distinct |
| D | broader_narrower or group_overlap |
| E | source_conflict or a C04/C05 flag (v0.3) |
| F | unresolved, or read the note |

The review labels whole groups, so the comparison uses each group's most severe relationship. Report agreement per category and list every disagreement with its reason. Also report where the tool is right for a different reason than the reviewer gave, since that shows which rule carried the result.

## Open questions

- Should a species-level name against one specific molecule (R05) ever count as a match for meta-analysis? Left to `merge_level` for now.
- Two Broad rows sit in the 3-methyladipic acid group. Rows from one file in one group may mean double counting downstream. A row check could catch it.
- Group-name rows at Species level, see Known limits.

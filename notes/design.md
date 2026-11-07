# Design

This file describes what the tool does, as of v0.11.0. Changes to behaviour are planned in `notes/decisions.md` first. This file is then updated to match the released code. Earlier versions of it are in the git history.

One principle runs through every rule below. When the evidence does not settle a question, the answer is `unresolved`. A rule that fixes one problem by guessing, and so creates a new one, does not go in. The one answer that is not about structure, `same_name`, says only that two names agree, and it is never counted as a claim.

## What the tool answers

Some matching step has put rows from two or more datasets into the same group. For each pair of rows in a group, the tool says how the two relate. MetLinkR is the first matching step it supports, through an adapter. The core does not depend on it.

Answers are per pair of rows, not per group. One group can hold a correct match, a wrong merge and a source error at the same time, as the pregnenolone sulfate and SM 36:1 groups in the first review did.

The tool also proposes lipid matches across files, including rows the matching step left in different groups, at the finest level both rows support.

"Same compound" is read at the level of the records. When one record gives stereo and the other does not, the pair is broader_narrower, even if both names mean the same natural compound (see Known limits).

## Input

One line per row:

| column | required | notes |
|---|---|---|
| row_id | yes | unique, such as `LEOCC_Metabolon:132` |
| file | yes | source dataset |
| original_name | yes | as written in the source file; an empty cell can be filled from a second column, `name_fallback_col` |
| name_source | no | `name_col`, `fallback`, or empty when neither column has a name |
| group_id | yes | the group the matching step assigned; for MetLinkR, `Standardized name` |
| inchikey | no | from the source file, never edited |
| smiles, inchi | no | from the source file, used for structure checks |
| pubchem, hmdb | no | used to find a key for rows without one |
| vendor | no | per file; `biocrates_p180` is the only value so far |

`read_metlinkr()` builds this table from the `*_metLinkR.csv` files.

## Evidence per row

Nothing is typed in by hand. Each value comes from the input file, from a fixed parsing rule, or from a versioned database file.

### Identifier resolution

A row without a valid InChIKey gets one from PubChem when it carries a single PubChem CID (tried first) or a single HMDB ID. The key is the one of PubChem's parent compound. These are not resolved:
- group names and alternative names;
- lipid names that are not one structure;
- internal standards;
- cells listing several IDs;
- an HMDB ID whose CIDs give different keys.

The file value stays in `inchikey_file`. Every key has an `inchikey_source` (`file`, `pubchem_cid` or `hmdb_via_pubchem`), and pair evidence names the source when a key rule used a derived key. PubChem answers are cached on disk.

A request that fails, through a network error or a PubChem error such as a busy server, is not cached. The rows that needed it get no key in that run, and a row whose parent could not be asked does not fall back to the CID's own key. `resolve_identifiers()` gives a message with the number of failed requests, and the audit keeps it in `info$pubchem_failed` and in `run_info.txt`. Running again asks the failed requests again. "Not found" is a real answer and is cached.

### InChIKey

- valid or not
- block1 (connectivity), block2, protonation letter
- stereo unspecified: block2 is `UHFFFAOYSA`

### Name checks

Regular expressions, one function each, each with tests.

- group name: " or ", "and or", "and/or", " + ", a list inside parentheses such as "(d18:1/20:2, d18:2/20:1)", or two words joined by a slash
- alternative names: "X or Y" at the top level
- tentative: trailing `*`
- isomer label (Metabolon): `(1)`, `(2)`, `[1]`, `[2]`
- isomer marker: the word "isomer", or a lone capital letter after the last semicolon (MS-DIAL)
- internal standard: `[iSTD]`, or an isotope label such as `15N4`, `13C6`, `d3`
- stereo stated: R/S, E/Z, cis/trans, D-/L- before a lowercase letter, optical signs
- stereo terms and omega positions, for comparing two names
- normalized name: lower case, spaces collapsed, trailing `*` and a leading optical rotation sign dropped

### Structure

When a row has an InChI or SMILES, `structure_evidence()` builds its InChI with RDKit and trusts it only when its key matches the row's stored key. Trusted structures allow the stereo comparison in R09b. Without RDKit the audit runs without these checks.

### Lipid evidence

Sources, in order:

1. **SwissLipids, by InChIKey.** For rows with a valid key, no group name and no alternative names.
   - The key is looked up in three tiers: full key, first 25 characters, first 14. The lookup stops at the first tier with a hit.
   - Species and Molecular subspecies ancestors are collected and parsed with rgoslin.
   - If no hit reaches a Species, every entry with the same first 14 characters is used, but only if they all agree on one Species.
2. **rgoslin, by name.** Only the Shorthand2020, Goslin and LipidMaps grammars count, at a level finer than class.
3. **LIPID MAPS (LMSD), by InChIKey.** Used only for rows that get nothing from the first two. Same three tiers. The row gets a Species only if every hit gives the same one.
4. **Vendor table.** Rows from a Biocrates p180 file get a Species-level profile from the vendor's own old-to-new nomenclature table.

When SwissLipids and the name both give a Species and they differ, the row gets no profile and is flagged. When they agree, or only the name parses, the profile follows the name, since the name says what was measured.

A group-name row gets a profile at level `group` only when its members share one chain sum and that sum matches the Species SwissLipids gives for its key.

Database files are identified by md5. A different file stops the run until the change is logged.

Each row's lipid profile has these fields:

- `lipid_species`, such as `SM 36:1;O2`
- `lipid_molecular`, if known
- `lipid_level`: `group`, `species`, `molecular`, `sn_position` or `structure`
- `lipid_source`: `swisslipids`, `swisslipids_fallback`, `rgoslin_name`, `lmsd`, `vendor_table` or `group_chain_sum`

## Pair rules

Rules run in this order. The first rule that applies decides, and every result keeps its rule ID.

| ID | condition | relationship |
|---|---|---|
| R02 | either row is an internal standard | distinct |
| R01s | R01 would apply, the rows are in different files, and the names are the same after normalization; neither name is a group name, a list of alternatives, an isomer label or marker, or a tentative name | same_name |
| R01 | either row has neither a valid InChIKey nor a lipid profile | unresolved |
| R03i | either name carries an isomer marker | unresolved |
| R04h | lipid Species differ only in hexose headgroup naming (HexCer vs GlcCer) | unresolved |
| R04 | lipid Species differ | distinct |
| R05g | same Species, one row is a group-name row and the other is at molecular level or finer | unresolved |
| R05 | same Species, reported at different levels | broader_narrower |
| R06 | same Species, both molecular or finer, molecular species differ | distinct |
| R03c | identical keys, or a pair R06b would call exact, but the names state different isomers | unresolved |
| R06b | same Species and same level, at species or molecular level, same name at that level | exact |
| R03 | identical InChIKeys | exact |
| R08 | either name reports a group of compounds | unresolved |
| R08b | either name lists alternative names | unresolved |
| R09 | same connectivity, exactly one key without stereo | broader_narrower |
| R09b | same connectivity, the stereo of one trusted structure is a subset of the other's | broader_narrower |
| R10 | same connectivity, different stereo, at least one name states stereo | stereo_variant |
| R11 | same connectivity, different stereo, no name states stereo | unresolved |
| R12n | different connectivity, but the same normalized name | unresolved |
| R12 | different connectivity | distinct |
| R99 | none of the above | unresolved |

R04 to R06b apply only when both rows have a lipid profile. R09 to R12 need a valid key on both rows.

Where the rules come from:

- R02: inosine-15N4 [iSTD] merged with inosine.
- R01s: COMETS, where 234 of the 250 unresolved pairs had no evidence on one side and 121 of those had the same name in two files. MetLinkR groups by name, so the label repeats part of its work. It sorts unresolved pairs; it does not check them.
- R03i: "Carnosine isomer A" carrying the key of carnosine (second blind test).
- R04: SM 36:1 merged with an SM 38:3 group.
- R05: Broad "SM 36:1;O2" against Metabolon "stearoyl sphingomyelin (d18:1/18:0)". Both resolve to Species `SM 36:1;O2`.
- R06 and R12: DG 18:1_18:2 [1] and [2] share Species `DG 36:3` but differ in structure. R06 does not separate them because their molecular species is the same, so R12 does, through the key.
- R03c: scyllo- and myo-inositol with one key without stereo, and DPA n3 and n6 with one key (COMETS).
- R08: "(2 or 3)-decenoate". The key of a group row is one member, not a confirmed structure.
- R11: 2-hydroxystearate as 2S in one file and 2R in the other, with no stereo in either name.
- R12n: fructose and 4-vinylphenol sulfate, where one identifier points to a different structure (COMETS).

A note on headgroups: SwissLipids writes glycosylceramides at Species level as HexCer or Hex2Cer, while a name may say GlcCer or LacCer. R04h keeps such pairs from being called distinct when the chain part is identical.

## Row checks

| ID | condition | flag |
|---|---|---|
| C01 | the InChIKey column holds a value that is not a standard InChIKey (a missing key is not flagged) | invalid_key |
| C03 | SwissLipids and the name give different Species | lipid_source_conflict |
| C04 | the row's own structure and its stored key differ in connectivity | structure_key_connectivity |
| C04 | same connectivity, but the stereo or isotope block differs | structure_key_stereo |

## Output

`audit_metlinkr()` returns, and `write_audit()` writes:

1. `pairs`: row_a, row_b, group_id, relationship, rule_id, evidence
2. `groups`: one line per group, with the most severe relationship and its rule
3. `row_flags`: row_id, check, flag, detail
4. `rows`: every row with its evidence
5. `lipid_matches`: lipid pairs across files, with relationship, shape (one_to_one or one_to_many) and whether the matching step had already grouped them

Severity, most severe first: distinct, unresolved, same_name, broader_narrower, stereo_variant, exact.

## Known limits

- broader_narrower from R09 means one record has stereo and the other does not. Often the names agree and the difference only reflects which PubChem record a lab chose. A curator matching by name would call such pairs the same compound. Changing this was tried in v0.8 and dropped (decisions.md, 2026-10-07).
- Sugars in ring and open-chain form get different keys. With the same name they stay unresolved (R12n); with different names they are called distinct.
- A sum name with a single-structure key ("ACar 17:0" with a straight-chain key) gets the key's answer.
- MS-DIAL-style lipid names are not read.
- Rows with only a name and no usable identifier get no evidence. A key looked up from the name was tried in v0.10 and dropped (decisions.md, 2026-11-07). Since v0.11 such a pair gets `same_name` when the names agree across files, which says nothing about structure.
- Two rows at structure level with the same Species stay unresolved (R99) when only one has a key, as with "1-stearoyl-GPC (18:0)" and "LPC 18:0/0:0".
- Glycosylceramide pairs that differ only in headgroup naming stay unresolved.
- Name checks are tuned on Metabolon, Broad, Biocrates, MS-DIAL and HMDB-style names.
- Nightingale class totals (Omega-6, total fatty acids) and other aggregate measures have no relationship of their own.
- The stereo comparison can report a false conflict for symmetric molecules such as inositol. A conflict never leads to a claim on its own.
- A structure whose SMILES and key agree but which does not match its name (the Metabolon steroids in the first review) is not detected.
- The Biocrates table covers the p180 kit as published in 2023; analytes it does not list get no profile.
- Pair checks cover only rows the matching step grouped. Lipid matches cover all lipid rows with a profile.

## Versions

- v0.1: audit of existing groups. InChIKey, name checks, SwissLipids, rgoslin. Rules R01 to R12, row checks C01 and C03.
- v0.2: lipid match proposals across files; group-name rows with a shared chain sum (R05g).
- v0.3: structure checks with RDKit (R09b, C04).
- v0.4: identifier resolution through PubChem. v0.4.1 skips internal standards.
- v0.4.2: R03c and R12n, after the COMETS test.
- v0.5: a check of rows against their RefMet entry. Failed its acceptance, not released, kept on branch `shelved/refmet-check`.
- v0.6: LIPID MAPS as a second lipid source and the Biocrates p180 vendor table.
- v0.7: R03i for isomer markers; optical rotation signs dropped in name normalization.
- v0.8: R09n planned and dropped at feasibility.
- v0.9: a second name column for rows with an empty name. No new pair rule. v0.9.1 reports failed PubChem requests instead of hiding them.
- v0.10: a key from the row's own name through PubChem, planned and dropped at feasibility. Kept on branch `shelved/name-key-lookup`.
- v0.11: R01s and the relationship `same_name`, for pairs without evidence whose names agree across files.

## Evaluation

Each test fixes its criteria and labels in `notes/decisions.md` before the tool runs on its data. Review categories map to relationships:

| review category | expected relationship |
|---|---|
| A | exact |
| B | stereo_variant |
| C | distinct |
| D | broader_narrower |
| E | not exact (source error) |
| F | unresolved |

A group's label is its most severe pair label, in the order E > C > F > D > B > A. The comparison uses each group's most severe relationship. `same_name` counts as unresolved does: a safe miss, never a definite claim.

# metMatchAudit 0.7.1

* No change in behaviour. The omega sign in the pattern behind `omega_positions()` is written as a `\u` escape, so the package code is plain ASCII and `R CMD check` passes without warnings.
* Comments in the package code point to the new script paths (`scripts/README.md`).

# metMatchAudit 0.7.0

* New rule R03i: when either name marks a feature of unknown identity, the pair stays unresolved. A marker is the word "isomer", as in "Carnosine isomer A", or a lone capital letter after the last semicolon, as in MS-DIAL's "CE 20:5; A". Such rows can carry the identifiers of the parent compound. The rule comes right after R01, before the lipid rules, because a lipid profile found through that key is no better than the key. New helper `has_isomer_marker()` and column `name_isomer_marker`.
* `normalize_name()` drops a leading optical rotation sign, so "(-)-N-Acetylneuraminic acid" and "N-Acetylneuraminic acid" meet in R12n.
* Written after the second blind test, where 11 of 17 wrong definite claims were `exact` calls on "isomer" features carrying their parent's key. All results on the test files after this change are after adaptation.
* Reading MS-DIAL lipid names was planned for this version and dropped: 74.9% of the test file's lipid names parsed after cleaning, short of the 80% set beforehand.

# metMatchAudit 0.6.0

* New `load_lmsd()` and `lmsd_lookup()`: the LIPID MAPS Structure Database as a second lipid source, looked up by InChIKey in three tiers like SwissLipids. `lipid_evidence()` uses it only for rows that get no profile from SwissLipids or from the name; `lipid_source` is then `"lmsd"`. Where both sources gave a Species, they agreed in 149 of 149 development rows and 66 of 66 COMETS rows.
* New `biocrates_p180_table()` and `vendor_lipid_profiles()`: rows from a file in Biocrates p180 nomenclature get a lipid profile from the vendor's own old-to-new table, always at species level. The table ships in `inst/extdata`, extracted from the vendor PDF by `scripts/10_v0.6/01_biocrates_table.R`.
* `read_metlinkr()` and `audit_metlinkr()` accept an optional `vendor` column; `audit_metlinkr()` takes `lmsd` and `lmsd_md5`.
* Development data: no pair result changed, evaluation as in 0.4.2. Lipid match proposals went from 63 to 110, all checked by hand and correct.
* COMETS, after adaptation: 23 pairs moved from unresolved to a definite answer, among them two MetLinkR merges of different compounds now called distinct.

There is no 0.5.0 release. The RefMet row check planned for it failed its acceptance test and is kept on branch `shelved/refmet-check`.

# metMatchAudit 0.4.2

* New rule R03c: identical InChIKeys, or a lipid annotation R06b would call exact, are not enough for `exact` when the names state different isomers (stereo descriptors such as myo and scyllo, or omega positions such as n3 and n6). The pair is unresolved.
* New rule R12n: the same name with different connectivity means the keys conflict. The pair is unresolved, not distinct.
* New name helpers `normalize_name()`, `stereo_terms()` and `omega_positions()`.
* Written after the blind COMETS test to fix its four wrong definite claims. On the development data the xylose pair is no longer wrong, and R03c caught a wrong `exact` call outside the reviewed groups. Evaluation: 35 agree, 8 unresolved, 0 wrong.

# metMatchAudit 0.4.1

* `resolve_identifiers()` skips internal standards, whose identifier usually belongs to the unlabelled compound. This fix was made on the development data before the COMETS test but missed the v0.4.0 tag.

# metMatchAudit 0.4.0

* New `resolve_identifiers()`: rows without a valid InChIKey get one from their PubChem CID or HMDB ID, through PubChem's parent compound. Group names, alternative names, ID lists and lipid names that are not one structure are skipped. Needs httr2 and jsonlite; answers are cached on disk.
* Every key has a source (`file`, `pubchem_cid`, `hmdb_via_pubchem`), and pair evidence names it when a derived key was used.
* New helpers `normalize_pubchem_id()` and `normalize_hmdb_id()`.
* `read_metlinkr()` accepts optional `pubchem_col` and `hmdb_col`.
* Tested on 783 rows not used in development: CID-derived keys agree with stored keys in 779 of 779 cases, HMDB-derived keys in 288 of 297.

# metMatchAudit 0.3.0

* New optional structure checks, using RDKit through reticulate. Without RDKit the audit runs as before.
* `structure_evidence()` builds an InChI for each row from its InChI or SMILES column and trusts it only when its key matches the row's stored InChIKey.
* `compare_stereo()` and `inchi_stereo()` compare the stereo layers of two InChIs.
* New rule R09b: when one row's stereo is a subset of the other's, the pair is broader_narrower. Conflicts make no new claim.
* New row check C04: a row's own structure and its stored key differ in connectivity or in stereo.
* `read_metlinkr()` accepts optional `smiles_col` and `inchi_col`.
* Evaluation against the manual review: 34 agree, 8 unresolved, 1 known limit (xylose). Lipid matches: 60 of 60 correct.

# metMatchAudit 0.2.0

* New `propose_lipid_matches()`: pairs lipid rows from different files that share a Species and reports exact and broader/narrower matches, marked one-to-one or one-to-many. Does not need MetLinkR. On the example data it proposes 59 matches, 48 of them not linked by MetLinkR; all 59 were checked by hand and found correct.
* Group-name lipid rows such as "sphingomyelin (d18:1/20:2, d18:2/20:1)*" get a profile at the new level "group" when the chain sum of their members matches the Species SwissLipids gives for their key.
* New rule R05g: a group row against a specific structure with the same Species stays unresolved.
* `audit_metlinkr()` returns `lipid_matches`; `write_audit()` writes `lipid_matches.csv`.
* Evaluation against the manual review: 29 agree, 13 unresolved, 1 known limit.

# metMatchAudit 0.1.0

First release.

* `audit_metlinkr()` checks every pair of rows that MetLinkR put in one group and reports their relationship with the rule behind it.
* Evidence from InChIKeys, name checks, SwissLipids and rgoslin.
* Evaluated on the MetLinkR example data against a blind manual review of 53 groups: 28 agree, 14 unresolved, 1 known limit (xylose), no source-error group called exact.

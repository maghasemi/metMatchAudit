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
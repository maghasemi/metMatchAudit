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
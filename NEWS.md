# metMatchAudit 0.1.0

First release.

* `audit_metlinkr()` checks every pair of rows that MetLinkR put in one group and reports their relationship with the rule behind it.
* Evidence from InChIKeys, name checks, SwissLipids and rgoslin.
* Evaluated on the MetLinkR example data against a blind manual review of 53 groups: 28 agree, 14 unresolved, 1 known limit (xylose), no source-error group called exact.
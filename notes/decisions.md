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

Clarification during review: some groups hold three or more rows with different relationships between pairs. The review sheet takes one category per group, so the group gets the most severe category present, in this order: E > C > D > B > A. The per-pair relationships go into evidence_note. This does not change any category already assigned.

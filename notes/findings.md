# Findings

Numbers here come from the scripts, not from memory. When a number is copied in, note which results file it came from.

## Phase 1a: InChIKey check on MetLinkR example output

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
| A | |
| B | |
| C | |
| D | |
| E | |
| F | |

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

### Open questions

- Resolved: the one invalid InChIKey in Broad_2022Aug is row 354, "Inosine-15N4 [iSTD]", where the InChIKey field holds the literal text "inchikey". MetLinkR gave this row the standardized name "Inosine", the same as natural inosine (UGQMRVRMYYASKQ-KQYNXXCUSA-N in both files). The early exploratory count read "inchikey" as a separate connectivity block and reported 17 different_connectivity groups. The format check in 01_build_rows.R excludes it, which gives 16.
- Two points for later. First, a placeholder value sits in a structure column of the source file, which fits category E. Second, MetLinkR grouped an isotope-labeled internal standard with the natural metabolite. Internal standards are not biological measurements and should not be merged with the analyte in a meta-analysis. The rule set needs a way to catch these. Name tokens such as "[iSTD]" or isotope labels such as "15N4" are candidates, but that has to be checked against more files first.
- "(N(1) + N(8))-acetylspermidine" (LEOCC row 1028) got the MetLinkR name "Asparagine". Its source IDs are comma-separated lists for two isomers and none of them refer to asparagine. Test whether list-valued ID cells cause wrong mappings in MetLinkR.
- Some source names report a group ("(2 or 3)-decenoate") while the InChIKey column holds one specific member. When a name contains "or" or "/", its InChIKey should be treated as one candidate, not as the confirmed structure.
- A group can hold more than one problem at once. In the pregnenolone sulfate group, two rows are a correct match, one row has an InChIKey that does not fit its own name (E), and that same row was wrongly merged (C). The review sheet allows one category per group. Decide before phase 2 whether the tool reports relationships per row pair rather than per group.
- Xylose appears as both its ring form and its open-chain form. These are the same compound, yet the InChIKeys differ in block1. A different connectivity block therefore does not always mean a different compound. The tool needs a way to recognize ring-chain forms and tautomers of the same compound. Otherwise it will raise false alarms on correct merges. If more cases like this turn up, consider a separate category after review.
- Several source InChIKeys carry stereo or position detail that the row name does not (2-hydroxyglutarate as 2S, "(2 or 3)-decenoate" as the 2-isomer). The tool should compare the detail level of the name with the detail level of the key, and not trust the key's extra detail by default.




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

### Open questions

- Resolved: the one invalid InChIKey in Broad_2022Aug is row 354, "Inosine-15N4 [iSTD]", where the InChIKey field holds the literal text "inchikey". MetLinkR gave this row the standardized name "Inosine", the same as natural inosine (UGQMRVRMYYASKQ-KQYNXXCUSA-N in both files). The early exploratory count read "inchikey" as a separate connectivity block and reported 17 different_connectivity groups. The format check in 01_build_rows.R excludes it, which gives 16.
- Two points for later. First, a placeholder value sits in a structure column of the source file, which fits category E. Second, MetLinkR grouped an isotope-labeled internal standard with the natural metabolite. Internal standards are not biological measurements and should not be merged with the analyte in a meta-analysis. The rule set needs a way to catch these. Name tokens such as "[iSTD]" or isotope labels such as "15N4" are candidates, but that has to be checked against more files first.
- "(N(1) + N(8))-acetylspermidine" (LEOCC row 1028) got the MetLinkR name "Asparagine". Its source IDs are comma-separated lists for two isomers and none of them refer to asparagine. Test whether list-valued ID cells cause wrong mappings in MetLinkR.
- Some source names report a group ("(2 or 3)-decenoate") while the InChIKey column holds one specific member. When a name contains "or" or "/", its InChIKey should be treated as one candidate, not as the confirmed structure.
- A group can hold more than one problem at once. In the pregnenolone sulfate group, two rows are a correct match, one row has an InChIKey that does not fit its own name (E), and that same row was wrongly merged (C). The review sheet allows one category per group. Decide before phase 2 whether the tool reports relationships per row pair rather than per group.
- Xylose appears as both its ring form and its open-chain form. These are the same compound, yet the InChIKeys differ in block1. A different connectivity block therefore does not always mean a different compound. The tool needs a way to recognize ring-chain forms and tautomers of the same compound. Otherwise it will raise false alarms on correct merges. If more cases like this turn up, consider a separate category after review.
- Several source InChIKeys carry stereo or position detail that the row name does not (2-hydroxyglutarate as 2S, "(2 or 3)-decenoate" as the 2-isomer). The tool should compare the detail level of the name with the detail level of the key, and not trust the key's extra detail by default.




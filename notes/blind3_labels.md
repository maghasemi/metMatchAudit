# Third blind test: labels

Written 2026-10-07, after labeling and before the audit is run on A, C1 and C2. The plan and criteria are in `notes/decisions.md` (2026-10-07, third blind test).

The sheets hold names from the private file A, so they stay in `data-raw/own/blind3/`. Their md5 is committed here instead, and `scripts/103_blind3_test.R` stops if either file differs.

| file | md5 |
|---|---|
| `data-raw/own/blind3/label_pairs.csv` | `898c093671f12638e61cf3637661ec3d` |
| `data-raw/own/blind3/label_lipids.csv` | `1a769432ed8b2c9d54d5f461dc2eeca2` |

## How the labels were made

A category and a reason for every line were first drafted with an AI assistant (see AI usage in the README). The proposals came from the names, the PubChem CIDs and HMDB IDs in the sheets, and the lipid naming rules. PubChem could not be queried at the time. Where CIDs differ, the reason says so and the label follows the names. The author accepted the proposals without changes. No tool output on these files had been seen.

## Stereo rule

The second test applied stereo inconsistently (sphinganine against S1P, see findings). This test fixes the rule before the run:

- A plain name ("alanine", "malate", "palmitoylcarnitine") means the natural form of the compound, by the convention of Metabolon, HMDB and most curated lists. A name that adds the natural form's descriptor ("L-alanine", "D-gluconate", "palmitoyl-L-carnitine") is the same compound: **A**.
- A descriptor that picks one of several forms found in nature, against a name that leaves the choice open, is broader or narrower: **D**. Examples are R or S on 3-hydroxybutyrate, erythro or threo, and cis against "cis or trans".
- When neither name states stereo and only the identifiers differ, the label is **A**.
- "L-Glycine" is glycine, which has no stereo.

This rule is expected to work against the tool in some pairs. The tool compares records, and a key without stereo against a key with it gives broader_narrower (rule R09). This was cause 3 in the second test. The rule was chosen because it is how a curator merging the two lists would decide, not to help the tool.

## Pairs

60 groups, 94 pairs, all across files. 31 pairs do not involve A.

| category | pairs | groups, by most severe label |
|---|---|---|
| A same compound | 91 | 58 |
| D broader or narrower | 3 | 2 |
| B, C, E, F | 0 | 0 |

The D pairs are:
- A's "aconitate [cis or trans]" against C2's cis-aconitate;
- A's "oleate/vaccenate (18:1)" against oleate in C1 and in C2.

All 60 groups are labeled A to D. Under C2, at most 3 may get a wrong definite answer, since 5% of 60 is 3.

The sample is dominated by A pairs, so C1 mostly tests whether the tool's `exact` calls are right. It says little about whether `distinct` calls are right.

## Lipid and carnitine rows

40 rows: 22 from C1 and 18 from C2. 17 have at least one match in another file, and every match is labeled A (same compound). 23 have none.

Three "none" labels were decided against a close neighbour:
- linoelaidyl (18:2 trans,trans) carnitine against A's linoleoylcarnitine (cis,cis);
- trans-2-hexadecenoylcarnitine against A's palmitoleoylcarnitine (9Z);
- C2's hydroxylated PE(18:0/20:4) against A's PE 18:0/20:4.

In the sheet, `matches_in_other_files` and `relationship` are lists separated by ";", aligned with each other.

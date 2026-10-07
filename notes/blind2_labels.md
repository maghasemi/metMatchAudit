# Second blind test: labels

Written 2026-10-07, after labeling and before the audit is run on A and B. The plan and criteria are in `notes/decisions.md` (2026-10-07, second blind test).

The label sheets hold names from the private file A, so they stay in `data-raw/own/blind2/` and are not committed. Their md5 is committed here instead. `scripts/93_blind2_test.R` stops if either file differs.

| file | md5 |
|---|---|
| `data-raw/own/blind2/label_pairs.csv` | `3f70bedbd9a0c5a1f21204571a37a75d` |
| `data-raw/own/blind2/label_lipids.csv` | `bd170037bd8d528ddfc5d984fc5e9813` |

## How the labels were made

Claude proposed a category and a reason for every line, working from the names and identifiers in the sheets and from PubChem, HMDB and the lipid naming rules. Claude saw no tool output on these files. The author went through the proposals and accepted them without changes. Neither of us has run the audit on A or B.

## Pairs

60 groups and 86 pairs: 72 pairs across files and 14 within one file.

| category | pairs | across files | groups, by most severe label |
|---|---|---|---|
| A same compound | 61 | 57 | 47 |
| C different compounds | 12 | 6 | 5 |
| D broader or narrower | 7 | 6 | 6 |
| F unclear | 6 | 3 | 2 |
| B stereo, E source error | 0 | 0 | 0 |

The C pairs are mostly B features named "isomer" set against their parent compound (adenosine, carnosine, N-acetylaspartate, hexose phosphate isomer D), plus oleoyl ethanolamide against oleoyl ethylamide. The D pairs are acylcarnitines where A names one chain and B gives only the sum, plus two sphinganine pairs. The F pairs are the 2-hydroxyglutarate isomer features and threonate against isothreonic acid.

### Group label

A group's label is its most severe pair label. F is ranked as follows:

E > C > F > D > B > A

A C pair settles the question of whether the group holds different compounds, whatever its other pairs are. An F pair next to D or A pairs leaves the group's true answer open. F groups are reported but are not among the "groups labeled A to D" that C2 counts.

This gives 58 groups labeled A to D. Under C2, at most 2 of them may get a wrong definite answer, since 5% of 58 is 2.9.

The plan did not fix F's place in this order. It is fixed here, before the run.

## Lipid rows

There are 40 lipid rows from B. 35 have no match in A. Five have a match, and all five are D:

| B row | A row | reason |
|---|---|---|
| B:102 ACar 14:0 | A:323 | B gives a sum, A the straight chain |
| B:575 LPE O-18:2 | A:12 | P-18:1 is one form of O-18:2 |
| B:953 PE 38:6 (18:2_20:4) | A:33 | same chains, A also gives positions |
| B:970 PE 40:6 | A:79 | B gives a sum, A gives chains |
| B:1037 PE O-38:5 (P-18:0_20:4) | A:18 | same chains, A also gives positions |

Three decisions are worth recording:

- B:529 LPC 18:2-SN2 is "none". Its acyl chain sits at sn-2, so it is not the same compound as A's 1-linoleoyl-GPC, and it is not a broader or narrower form of it either.
- B:575 is D, not A. O-18:2 also covers an ether lipid with the double bond elsewhere in the chain.
- B:953 is matched only to the A row with the same chains. A PE 38:6 with other chains would be a different compound from B's annotation.

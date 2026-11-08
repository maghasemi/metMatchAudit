# Fourth blind test: labels

Written 2026-11-08, after the user accepted the labels and before the audit is run on F1, F2 and F3. The plan and criteria are in `notes/decisions.md` (2026-11-08, fourth blind test, plan). The draft and its judgment calls are in the entry "step 4, the sample and the draft labels".

The files are CC BY-NC-ND, so the sheets stay in `data-raw/blind4/labels/`, which git ignores. Their md5 is committed here instead, and `scripts/19_blind_test_4/04_test.R` must stop if either file differs.

| file | md5 |
|---|---|
| `data-raw/blind4/labels/label_pairs.csv` | `f9f03def8a5b04546c8839673ea6d811` |
| `data-raw/blind4/labels/label_lipids.csv` | `70e2213390b427e7d58adf24ec81a051` |

## How the labels were made

The assistant drafted a category and a reason for every line from the names, the PubChem CIDs and HMDB IDs in the sheets, and a few PubChem records read for the hard cases (CID 5280818, 5280819, 44457542, 442106, 5283572 and others). The user accepted all of them without changes, and cannot check the chemistry. This is the weakest part of the test, as the plan says. No tool output on these files had been seen. The reason of the 27 same-name pairs is a template that only names the shared name and the identifiers.

## Stereo rule

Unchanged from `notes/blind3_labels.md`: a plain name means the natural form; a name that adds the natural form's descriptor is the same compound (A); a descriptor that picks one of several natural forms against a name that leaves the choice open is broader or narrower (D); when neither name states stereo and only identifiers differ, A.

## Pairs

57 groups (27 same name, 30 different names), 67 pairs, 64 across files.

| category | pairs |
|---|---|
| A same compound | 55 (27 same name, 28 different names) |
| C different compounds | 8 |
| D broader or narrower | 3 |
| F unclear | 1 |
| B, E | 0 |

The F pair is F1:21 against F2:695, a plasmalogen PC whose ether chain is at sn-1 in F1 and, by F2's order of chains, at sn-2. PubChem has no record under either name. Under C2 an F group sits between C and D in the order of severity. The same-name stratum has 27 groups, below the 30 planned, so the report says so.

## Lipid and carnitine rows

40 rows: 4 from F1, 36 from F2, none from F3. 6 have an A match in another file, 6 a D match, 28 none. The F2:726 row is labeled none for the same doubt as the F pair.

## Judgment calls the user accepted

- "Urobilin" against L-urobilin is C, because the formulas in PubChem differ.
- F2:277 PC(18:0) against F1:107 is D (position open against sn-1).
- F2:723 against F1:44 is D (double bond positions open in F1).

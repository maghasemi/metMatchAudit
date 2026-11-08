# Evaluation

Every test of the tool, in the order it was run. Each test fixed its criteria and labels in [decisions.md](decisions.md) before the tool ran on its data. Results after a change are marked as after adaptation. The README gives a summary.

## Manual review of the example data

I first checked the MetLinkR example output by hand. Of 53 groups whose rows carried different InChIKeys:

| what the rows turned out to be | groups |
|---|---|
| one annotation broader than the other (for example SM 36:1;O2 against SM d18:1/18:0) | 27 |
| different compounds merged under one name | 13 |
| an error in a source file | 10 |
| stereoisomers | 2 |
| unclear | 1 |

Some of the wrong merges matter biologically: omega-3 and omega-6 docosapentaenoic acid, the endocannabinoid 2-arachidonoylglycerol and its isomerisation product 1-arachidonoylglycerol, 1,2- and 1,3-diacylglycerols, and asparagine with acetylspermidine. Running the audit on the rest of the example data showed internal standards merged with the metabolites they stand in for (inosine-15N4, glycocholate-d4, thymine-d4).

The source errors were a surprise. Five steroid structures in the Metabolon example file have inverted stereocentres in the ring system, so the file's own InChIKey and SMILES describe a different stereoisomer than the name (dehydroepiandrosterone sulfate, glycochenodeoxycholate, glycolithocholate, taurocholate, taurodeoxycholate). Two rows in the Broad file carry the name of one bile acid and the structure of another.

## Development data against the manual review

The manual review of the 53 groups was finished and committed before the tool's output on those groups was opened. Each version was compared with it group by group:

| outcome | v0.1 | v0.2 | v0.3 | v0.4 | v0.4.2 and v0.6 |
|---|---|---|---|---|---|
| agrees with the review | 28 | 29 | 34 | 34 | 35 |
| unresolved where the review has an answer | 14 | 13 | 8 | 8 | 8 |
| a different definite answer | 1 | 1 | 1 | 1 | 0 |
| source-error groups called exact | 0 of 10 | 0 of 10 | 0 of 10 | 0 of 10 | 0 of 10 |

The one wrong answer up to v0.4 was xylose, whose ring and open-chain forms have different InChIKeys. Since v0.4.2, two rows with the same name and conflicting keys stay unresolved, which is what the review says for that group.

The 8 unresolved cases are 5 groups with a group-name row the tool cannot check yet, and 3 where stereo really differs but no name states it, or where a row's structure does not match its own key. By design these stay unresolved.

Lipid matches on the example data: 110 proposals between the Metabolon and Broad files, all checked by hand and all correct ([v02_proposal_review.csv](v02_proposal_review.csv), [v06_dev_lipid_review.csv](v06_dev_lipid_review.csv)).

## External test on COMETS

To see how the tool does on data it was not built on, I ran MetLinkR on five COMETS datasets and compared the result with the COMETS expert mappings from the MetLinkR supplement. The criteria were written and committed before the run, and the expert labels were not opened until the tool's output was final. Our MetLinkR run found 183 of 214 multi-file expert groups whole (85.5%; the paper reports 85.3%).

Blind result, v0.4.1:

| criterion | result | |
|---|---|---|
| C1: definite claims (exact, distinct) agree with the expert, at least 90% | 121 of 125 (96.8%) | met |
| C2: with the audit, precision and F1 not lower than MetLinkR alone | F1 0.893 against 0.895 | not met |

The audit's definite claims were mostly right, but it made few of them. 259 of 401 pairs stayed unresolved, mostly because three of the five files carry names only. Its two `distinct` calls were both wrong, and it separated none of MetLinkR's real wrong merges.

Reading the 42 pairs the expert does not group, after the result was sealed: 7 are different compounds, 4 stereoisomers, 17 a broader against a narrower annotation or a class total, 9 the same compound measured on a different platform that the expert keeps apart, and 5 look like inconsistencies in the expert labels (Nightingale alanine is grouped with beta-alanine, for example). The expert answers "same measurement"; this tool answers "same compound". Details are in [findings.md](findings.md) and [comets_review.csv](comets_review.csv).

Changes made after that test, and their results on COMETS. These numbers are after adaptation, not blind:

- v0.4.2: identical keys are not enough for `exact` when the names state different isomers (R03c), and the same name with conflicting keys stays unresolved (R12n). The four wrong claims of the blind run became unresolved. On the example data, R03c also caught a wrong `exact` call no review had seen, (R)- against (S)-3-hydroxybutyrylcarnitine.
- v0.6: LMSD as a second lipid source and the Biocrates p180 table. 23 COMETS pairs moved from unresolved to a definite answer, 22 of them across files. Two became `distinct`, valerylcarnitine against glutarylcarnitine and C14:1 against C14:2, both MetLinkR merges of different compounds; one of the two is a pair across files. Across files, definite claims now agree with the expert in 122 of 122, precision goes from 0.895 to 0.897 and F1 from 0.895 to 0.896 (recounted by `../scripts/08_comets/04_recount_all_versions.R`; an earlier version of this README gave 123, 0.900 and 0.897, worked out by hand instead of counted). 36 new lipid matches across files; the expert groups 18 of them, and the other 18 are Biocrates sums matched to single structures, which the expert links one to one or not at all.

## Second blind test

v0.6.0 ran on two files it had not seen. One is my own Metabolon annotation of a mouse dataset (494 rows, private; only counts are published). The other is the annotation from the aging mouse brain atlas of Ding et al. 2021 (1547 rows, MS-DIAL-style names). Before the run I labeled 86 pairs in 60 randomly drawn MetLinkR groups and 40 random B lipid rows, then committed the criteria and the md5 of the label files.

| criterion | result | |
|---|---|---|
| C1: definite claims (exact, distinct) agree with the label, at least 90% | 48 of 65 (73.8%) | not met |
| C2: groups labeled A to D with a wrong definite answer, at most 5% | 12 of 58 (20.7%) | not met |
| C3: lipid proposals involving a sampled B row are in the labels, at least 90% | 1 of 1 | met, on one proposal |
| C3 recall | 1 of 5 | |

11 of the 17 wrong definite claims come from one habit of the B file. Features named "Adenosine isomer", "Carnosine isomer A" and the like carry the parent compound's identifiers, so the tool called them `exact` with the parent. Most of the rest are rows whose name and key disagree: "ACar 17:0" names a sum but carries a straight-chain key. The tool follows the key, and the labels follow the name. Lipid recall was low because the tool could not read MS-DIAL names such as "PE 38:6; PE 18:2_20:4;". 375 of 960 B lipid rows got no profile, and the rest got a structure-level profile from their key that their names do not support.

The tool held up on the earlier test sets but did not carry over to an unseen naming style. Running v0.4.1 on the same files gave 48 of 67 and 13 of 58, so the later versions did not cause the drop. Details: [findings.md](findings.md).

v0.7, after adaptation: a new rule R03i keeps a pair unresolved when a name marks a feature of unknown identity, either by the word "isomer" or, in MS-DIAL lipid names, by a lone letter after the last semicolon ("PC 36:2; A"). On the same files, C1 goes to 48 of 53 (90.6%) and C2 to 8 of 58. The development data and COMETS do not change. These numbers come from the files the rule was written for, so they show the fix works, not that it generalizes. Reading MS-DIAL lipid names was planned for v0.7 and dropped: only 74.9% of the names parsed, short of the 80% set beforehand.

## Third blind test

v0.7.0 ran on my Metabolon file and on two published lists from Lai et al. 2021 (Nat. Commun.; mouse feces, 533 rows, and serum, 231 rows, each with PubChem CIDs and HMDB-style names). As before, I drew 60 MetLinkR groups and 40 lipid and carnitine rows at random, labeled them, and committed the criteria and the label md5s before the run. The rule for stereo in names was also fixed in advance.

| criterion | result | |
|---|---|---|
| C1: definite claims (exact, distinct) agree with the label, at least 90% | 78 of 78 (100%) | met |
| C2: groups labeled A to D with a wrong definite answer, at most 5% | 7 of 60 (11.7%) | not met |
| C3: lipid proposals involving a sampled row are in the labels, at least 90% | 13 of 13 | met |
| C3 recall | 13 of 20 | |

Every `exact` call was right, and no lipid proposal was wrong. All seven wrong groups have one cause. Two rows name the same compound, but one lab's PubChem record has stereo and the other's does not, and rule R09 then calls the pair broader_narrower. In five of the seven, neither name says anything about stereo. The same cause appeared in the second test. The sample was mostly same-compound pairs, so it tests `exact` calls far more than `distinct` ones. Details: [findings.md](findings.md).

## A second name column (v0.9)

In the Broad example file, `Metabolite` is empty on 119 of 511 rows and the name is in column `name`. v0.9 can read that column, and these results are after adaptation.

On the development data, unresolved pairs went from 39 to 30, and lipid proposals went from 110 to 129. The evaluation against the review did not change, because none of the 13 affected groups is in it. COMETS and the blind test files have no blank names, so nothing changed there.

The 19 new claims passed a scripted check of chain sums and formulas. A curator did not review them, and Broad has no m/z on these rows.

## Pairs that only share a name (v0.11)

v0.11 gives a pair the relationship `same_name` when R01 would leave it unresolved, the rows are in different files and their names agree after normalization. No structure is compared, and MetLinkR groups by name, so this sorts unresolved pairs and is not a claim. These results are after adaptation.

On COMETS, 121 of the 234 R01 pairs moved to `same_name`. The expert puts 120 of them in one group. The one it separates is "N-acetylcarnosine" against "N-acetylcarnosine". Definite claims, lipid proposals and the development evaluation did not change. The two blind tests have only six R01 pairs between them, so a test of this label needs new data with name-only files.

## What did not work

Three ideas were tried and set aside. All are logged in [decisions.md](decisions.md).

- Giving name-only rows a key from RefMet. MetLinkR's standardized names are RefMet names, so a key from RefMet would agree with MetLinkR by construction and audit nothing.
- Checking each row against the RefMet entry of its group. The check found real problems, among them a Metabolon row for 2-hydroxycarbamazepine filed under norcodeine and two RefMet records with wrong formulas (CAR 20:4 stored as C31H54NO4). But only 17 of 40 flags were wrong assignments, short of the 80% set beforehand, and the check could not change any pair result. The code is kept on branch `shelved/refmet-check`.
- Giving name-only rows a key from PubChem, looked up by the row's own name. On COMETS, 142 of the 234 R01 pairs would have become exact, and 140 of them agree with the expert. But on rows that already have a key, the name's key had the same connectivity in only 95.3%, short of the 98% set beforehand. The code is kept on branch `shelved/name-key-lookup`.

## Open problems


- A stereo-only key difference between two names that state no stereo should not give broader_narrower. This was the only cause of wrong answers in the third test and one of the causes in the second.
- A further blind test with more pairs of different compounds, since the third test had almost none.
- Reading MS-DIAL lipid names, with a class map for names such as Cer-NS and HexCer-NDS.
- A name against key check for lipids. In the second test, B names PC 35:3 for a row whose key is PE 38:3.
- Nightingale biomarker definitions and an aggregate relationship for class totals.
- Checking which structures a group-name row lists.

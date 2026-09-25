# How to fill in review.csv

`scripts/03_make_review_sheet.R` creates `notes/review.csv`. Each line is one MetLinkR name whose rows carry more than one InChIKey.

For each line:

1. Open every link in `pubchem_links`. Each one searches PubChem for one InChIKey.
2. Compare what PubChem says each key is with the names in `original_names`.
3. Pick one category from `notes/decisions.md` (A to F) and put the letter in `category`.
4. In `evidence_note`, write what decided it in a few words, for example "both CID records are bilirubin, E,E vs Z,Z double bond geometry" or "C15:0 vs C15:2 side chain, different compounds".

If a key returns nothing in PubChem, write that in `evidence_note` and use F unless something else in the row settles it.

Do not change any other column. Do not reorder lines. Save as CSV with UTF-8 encoding. If you edit in Excel, check afterwards that it did not turn InChIKeys or IDs into something else.

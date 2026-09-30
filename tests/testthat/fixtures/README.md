# Test fixtures

`swisslipids_mini.tsv` has the same header as the SwissLipids download (`lipids.tsv`) and 25 made-up entries. Entry IDs start at SLM:000000901 so they cannot be mistaken for real ones.

The InChIKeys attached to the isomeric entries are real keys from the MetLinkR example files, so tests can use real row names:

| key | row in the example data |
|---|---|
| BLZVZPYMHLXLHG-RQOIEFAZSA-N | oleoyl-linoleoyl-glycerol (18:1/18:2) [1], LEOCC |
| GREDRAMJRDQWEJ-LTEAFHAISA-N | oleoyl-linoleoyl-glycerol (18:1/18:2) [2], LEOCC (stored here with protonation letter O, to test tier 2) |
| RZRNAYUHWVFMIP-KTKRTIGZSA-N | 1-oleoylglycerol (18:1), LEOCC |
| LKQLRGMMMAHREN-YJFXYUILSA-N | stearoyl sphingomyelin (d18:1/18:0), LEOCC, and SM 36:1;O2, Broad |
| VJLLLMIZEJJZTE-NNTBDIJYSA-N | glycosyl-N-palmitoyl-sphingosine (d18:1/16:0), LEOCC |
| YVTVVJRJLRHYSJ-COUPUTHFSA-N | sphingomyelin (d18:1/20:2, d18:2/20:1, d16:1/22:2)*, LEOCC (placed under SM(d38:3)) |
| IWHLTJGDYLJQSH-FLPFCPAFSA-N | sphingomyelin (d18:2/16:0, d18:1/16:1)*, LEOCC (placed under SM(d34:2)) |

The hierarchy copies situations seen in the feasibility run: an entry with no Level and no parent that needs the connectivity fallback (1-oleoylglycerol), a connectivity block whose entries point to two different Species (no fallback allowed), a glycosylceramide under a HexCer Species, and a loop in the parent links.

## metlinkr/

Two small files in the layout MetLinkR writes (`<name>_metLinkR.csv`, with a `Standardized name` column). Rows and InChIKeys are copied from the MetLinkR example output (LEOCC Metabolon and Broad 2022 files). The standardized names follow what MetLinkR assigned there. "Glycine" has "-" as its standardized name, which MetLinkR uses internally for unmapped rows.
